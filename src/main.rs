#![allow(dead_code)]

use std::f32::consts::{PI};

use image::{ RgbImage, Rgb32FImage, Rgb};
use core::ops::Index;

mod sampling;
mod vec3;

use vec3::{Vec3, dot, cross};
use sampling::*;

// const K_1_PI : f32 = 0.318309886;
const LUT_RES :u32 = 32;

#[derive(Copy, Clone)]
struct Point([f32; 2]);


impl Index<usize> for Point {
    type Output = f32;
    #[inline(always)]
    fn index(&self, _index: usize) -> &f32 {
        &self.0[_index]
    }
}

fn hammerslay(i : usize, n_bits : usize) -> Point
{
    let k = 1 << n_bits;
    let mut x = 0.0;
    let y = (i as f32 + 0.5) / k as f32;
    {
        let mut n = i;
        let mut base = 1.0/ 2.0;
        while n != 0
        {
            if n & 1 != 0
            {
                x += base;
            }
            n /= 2;
            base /= 2.0;
        }
    }
    return Point([x, y]);
}

fn v3_from_sample(s : &Sample3D) -> Vec3
{
    Vec3::new(s[0], s[1], s[2])
}

struct TangentSpace
{
    t : Vec3,
    b : Vec3,
    n : Vec3,
}

impl TangentSpace
{
    fn from_normal(n: &Vec3) -> TangentSpace
    {
        let up = if n.z.abs() < n.x.abs() { Vec3::new(0.0, 0.0, 1.0) } else { Vec3::new(0.0, 1.0, 0.0) };
        let tx = cross(&up, &n).normalize();
        let ty = cross(&n, &tx);
        TangentSpace{t: tx, b:ty, n: *n}
    }

    fn transform(&self, v : &Vec3) -> Vec3
    {
        return self.t.scale(v.x) + self.b.scale(v.y) + self.n.scale(v.z);
    }
}

fn fresnel_schlick(f0 : f32, f90 : f32, u : f32) -> f32
{
    f0 + (f90 - f0) * (1.0 - u).powf(5.0)
}

fn diffuse_lambert() -> f32
{
    1.0 / PI
}

fn diffuse_burley(n_v : f32, n_l : f32, l_h : f32, a : f32) -> f32
{
    let fd90 = 0.5 + 2.0 * l_h * l_h * a;
    let l_scatter = fresnel_schlick(1.0, fd90, n_l);
    let v_scatter = fresnel_schlick(1.0, fd90, n_v);
    l_scatter * v_scatter / PI
}

fn d_ggx(n_h : f32, a2 : f32) -> f32
{
    let f = n_h * n_h * (a2 - 1.0) + 1.0;
    a2 / (f * f) / PI
}

fn g_smith_1(n_x : f32, a2: f32) -> f32
{
    n_x / (n_x * (1.0 - a2) + a2)
}

// Smith uncorrelated
fn g_smith_u(n_l : f32, n_v : f32, a2 : f32) -> f32
{
    g_smith_1(n_l, a2) * g_smith_1(n_v, a2)
}

// Visibility function simplify with the 4 (n_v)(n_l) denominator
fn v_smith_1(n_x : f32, a2: f32) -> f32
{
    let n_x2 = n_x*n_x;
    1.0 / (n_x + (a2 + (1.0 - a2) * n_x2).sqrt())
}

// Smith uncorrelated
fn v_smith_u(n_l : f32, n_v : f32, a2 : f32) -> f32
{
    v_smith_1(n_l, a2) * v_smith_1(n_v, a2)
}

// Smith correlated
fn v_smith(n_l : f32, n_v : f32, a2 : f32) -> f32
{
    let s = 1.0 - a2;
    let v = n_v * (n_l * n_l * s + a2).sqrt();
    let l = n_l * (n_v * n_v * s + a2).sqrt();
    return 0.5 / (v + l);
}

fn brdf_d(n_v : f32, n_l : f32, v_h : f32, a : f32) -> f32
{
    // diffuse_lambert()
    diffuse_burley(n_v, n_l, v_h, a)
}

fn brdf_s(n_v : f32, n_l : f32, n_h : f32, v_h : f32, a : f32, f0 : f32) -> f32
{
    let a2 = a*a;
    let d = d_ggx(n_h, a2);
    let f = fresnel_schlick(f0, 1.0, v_h);
    let g = v_smith_u(n_l, n_v, a2);
    d * f * g // (4.0 * n_h * n_v) // included in v_smith
}

fn encode_srgb(v : f32) -> u8
{
    // let t = if v < 0.0031308 {
    //     v * 12.92
    // } else {
    //     1.055 * v.powf(1.0/2.4) - 0.055
    // };
    let t = v;
    return (t.clamp(0.0, 1.0) * 255.0) as u8;
}

fn tonemap(c : f32) -> f32
{
    (c*(2.51*c+0.03))/(c*(2.43*c+0.59)+0.14)
}

fn furnace_test(mu: f32, a2: f32, samples: &[Vec3]) -> f32
{
    assert!(mu >= 0.0 && mu <= 1.0);
    let mut acc = 0.0;
    let n = Vec3::new(0.0, 0.0, 1.0);
    let view_dir = Vec3::new((1.0-mu*mu).sqrt(), 0.0, mu);
    for xi in samples
    {
        let h = xi;
        let light_dir = 2.0 * dot(&view_dir, h) * *h - view_dir;
        let n_l = dot(&n, &light_dir);
        if n_l > 0.0
        {
            let n_v = mu;
            let v_h = dot(&view_dir, &h).clamp(0.0, 1.0);
            let n_h = dot(&n, &h).clamp(0.0, 1.0);
            let d = 1.0;//d_ggx(n_h, a2); // simplify with the pdf
            let f = 1.0;//fresnel_schlick(f0, 1.0, v_h); // Always 1.0
            let g = v_smith(n_l, n_v, a2);
            let s = d * f * g; // (4.0 * n_h * n_v) simplify with V
            // let s = brdf_s(n_v, n_l, n_h, v_h, a, f0);
            let pdf = 4.0 * v_h / n_h;
            // let pdf = 1.0 / TAU; // pdf is d_ggx() cos(theta) sin(theta)
            // let pdf = 1.0; // pdf is d_ggx() cos(theta) sin(theta)
                           // but it simplify with s that should be multiplied by cos(theta)
            acc += s * n_l * pdf;// / (4.0 * n_v);
        }
    }
    acc / samples.len() as f32
}

fn furnace_test_with_lut(mu: f32, rough: f32, samples: &[Vec3], lut: &Rgb32FImage) -> f32
{
    assert!(mu >= 0.0 && mu <= 1.0);
    let a = rough * rough;
    let a2 = a * a;
    let mut acc = 0.0;
    let n = Vec3::new(0.0, 0.0, 1.0);
    let view_dir = Vec3::new((1.0-mu*mu).sqrt(), 0.0, mu);
    for xi in samples
    {
        let h = xi;
        let light_dir = 2.0 * dot(&view_dir, h) * *h - view_dir;
        let n_l = dot(&n, &light_dir);
        if n_l > 0.0
        {
            let n_v = mu;
            let v_h = dot(&view_dir, &h).clamp(0.0, 1.0);
            let n_h = dot(&n, &h).clamp(0.0, 1.0);
            let d = 1.0;//d_ggx(n_h, a2); // simplify with the pdf
            let f = 1.0;//fresnel_schlick(f0, 1.0, v_h); // Always 1.0
            let g = v_smith(n_l, n_v, a2);
            let s = d * f * g; // (4.0 * n_h * n_v) simplify with V
            let pdf = 4.0 * v_h / n_h;
            acc += (s + f_multi_scatter(n_v, n_l, rough, lut)) * n_l * pdf;
        }
    }
    acc / samples.len() as f32
}


fn lerp(a : Rgb<f32>, b : Rgb<f32>, t : f32) -> Rgb<f32>
{
    let r = a[0] + (b[0] - a[0]) * t;
    let g = a[1] + (b[1] - a[1]) * t;
    let b = a[2] + (b[2] - a[2]) * t;
    Rgb([r, g, b])
}

fn sample_lut(u : f32, v : f32, lut: &Rgb32FImage) -> Rgb<f32>
{
    let px = (u as u32);
    let px_next = (px + 1).min(LUT_RES - 1);
    let py = (v as u32);
    let py_next = (py + 1).min(LUT_RES - 1);

    let a = *lut.get_pixel(px, py);
    let b = *lut.get_pixel(px_next, py);
    let c = *lut.get_pixel(px, py_next);
    let d = *lut.get_pixel(px_next, py_next);

    let ku = u - u.floor();
    let kv = v - v.floor();

    let e = lerp(a, b, ku);
    let f = lerp(c, d, ku);
    lerp(e, f, kv)
}

fn f_multi_scatter(n_v : f32, n_l : f32, rough : f32, lut: &Rgb32FImage) -> f32
{
    let mu_v = n_v * LUT_RES as f32;
    let mu_l = n_l * LUT_RES as f32;
    let a = rough * LUT_RES as f32;
    let Ev = sample_lut(mu_v, a, lut);
    let El = sample_lut(mu_l, a, lut);
    let E_avg = PI * Ev[1];
    Ev[0]*El[0]/E_avg
}

fn render_sphere( img: &mut RgbImage, rough: f32, f0: f32, lut: &Rgb32FImage)
{
    let dx = 1.0 / img.width() as f32;
    let dy = 1.0 / img.height() as f32;
    let g = encode_srgb(0.5);
    let black = Rgb([0, 0, 0]);
    let gray = Rgb([g, g, g]);

    let view_dir = Vec3::new(0.0, 0.0, -15.0).normalize();
    let light_dir = Vec3::new(-4.0, 4.0, -2.0).normalize();

    let a = rough * rough;
    let a2 = a * a;

    for j in 0..img.height()
    {
        let v = (0.5 + j as f32) * dy;
        let y = v * -2.0 + 1.0;
        for i in 0..img.width() {
            let u = (0.5 + i as f32) * dx;
            let x = u * 2.0 - 1.0;
            let mut out_color = if y > 0.0 { black } else { gray };
            let delta = x*x + y*y;
            if delta < 1.0
            {   // inside the sphere
                let z = -(1.0 - delta).sqrt();
                let n = Vec3::new(x, y, z);
                // let furnace = furnace_test(-z, a2, &ggx_samples);
                let mut result = 0.0;

                let h = (view_dir + light_dir).normalize();
                let n_l = dot(&n, &light_dir);
                if n_l > 0.0
                {
                    let n_v = dot(&n, &view_dir).clamp(0.0, 1.0);
                    let v_h = dot(&view_dir, &h).clamp(0.0, 1.0);
                    let n_h = dot(&n, &h).clamp(0.0, 1.0);
                    let d = d_ggx(n_h, a2);
                    let f = fresnel_schlick(f0, 1.0, v_h); // Always 1.0
                    let g = v_smith(n_l, n_v, a2);
                    let s = d * f * g; // (4.0 * n_h * n_v) simplify with V
                    // result = (s) * n_l;
                    result = (s + f_multi_scatter(n_v, n_l, rough, lut)) * n_l;
                }

                let t = tonemap(result);
                let r_ = encode_srgb(t);
                let g_ = r_;//encode_srgb(y*0.5+0.5);
                let b_ = r_;//encode_srgb(g*0.5+0.5);
                out_color = Rgb([r_, g_, b_]);
            }
            img.put_pixel(i, j, out_color);
        }
    }
}

fn main() {
    const BIT_DEPTH :usize = 12;
    const SAMPLE_COUNT :usize = 1 << BIT_DEPTH;
    let mut samples : Vec<Point> = Vec::new();
    for k in 0..SAMPLE_COUNT {
        samples.push(hammerslay(k, BIT_DEPTH));
    }
    let mut lut = Rgb32FImage::new(LUT_RES, LUT_RES);
    let lut_scale = 1.0 / LUT_RES as f32;
    for j in 0..LUT_RES {
        let rough = (0.5 + j as f32) * lut_scale;
        let a = rough * rough;
        let a2 = a * a;
        let mut ggx_samples : Vec<Vec3> = Vec::new();
        for xi in &samples {
            ggx_samples.push(v3_from_sample(&ggx_sampling(xi[0], xi[1], a2)));
        }
        let mut e_agv = 0.0;
        for i in 0..LUT_RES {
            let mu = (0.5 + i as f32) * lut_scale;
            e_agv += furnace_test(mu, a2, &ggx_samples);
        }
        e_agv /= LUT_RES as f32;

        for i in 0..LUT_RES {
            let mu = (0.5 + i as f32) * lut_scale;
            let f = 1.0 - furnace_test(mu, a2, &ggx_samples);
            lut.put_pixel(i, j, Rgb([f, e_agv, f]));
        }
    }
    lut.save("lut.exr").expect("Saving img failed");

    for i in 0..LUT_RES {
        let mu = (0.5 + i as f32) * lut_scale;
        let rough = 0.5;
        let a = rough * rough;
        let a2 = a * a;
        let mut ggx_samples : Vec<Vec3> = Vec::new();
        for xi in &samples {
            ggx_samples.push(v3_from_sample(&ggx_sampling(xi[0], xi[1], a2)));
        }
        let f = furnace_test_with_lut(mu, rough, &ggx_samples, &lut);

        println!("mu {} -> {}", mu, f);
    }

    // const OUT_RES :u32 = 512;

    // for i in 0..10 {
    //     let mut out_img = RgbImage::new(OUT_RES, OUT_RES);
    //     let rough = (i as f32 + 0.5) / 10.0;
    //     render_sphere(&mut out_img, rough, 1.0, &lut);
    //     let img_name = format!("img/sphere{:02}.png", i);
    //     println!("Saving {}", img_name);
    //     out_img.save(img_name).expect("Saving img failed");
    // }
}
