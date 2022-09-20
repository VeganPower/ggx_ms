use std::f32::consts::{PI, TAU};
use core::ops::Index;

pub struct Sample3D([f32; 3]);

impl Index<usize> for Sample3D {
    type Output = f32;
    #[inline(always)]
    fn index(&self, _index: usize) -> &f32 {
        &self.0[_index]
    }
}

pub fn from_polar(c_theta: f32, phi: f32) -> Sample3D
{
    let c_theta2 = c_theta * c_theta;
    let s_theta = (1.0 - c_theta2).sqrt();
    Sample3D([s_theta * phi.cos(), s_theta * phi.sin(), c_theta])
}

// PDF = 1.0 / (pi * pi)
pub fn uniform_sampling(u: f32, v: f32) -> Sample3D
{
    // theta = v * pi / 2
    let phi = u * TAU;
    let theta = v * PI * 0.5;
    let c_theta = theta.cos();
    from_polar(c_theta, phi)
}

// PDF = sin(theta) / (2.0 pi)
pub fn hemisphere_sampling(u: f32, v: f32) -> Sample3D
{
    // theta = acos(1-v)
    let phi = u * TAU;
    let c_theta = 1.0 - v;
    from_polar(c_theta, phi)
}

// PDF = sin(theta) * cos(theta) / pi
pub fn cos_weighted_sampling(u: f32, v: f32) -> Sample3D
{
    // theta = asin(sqrt(v))
    let phi = u * TAU;
    let s_theta = v.sqrt();
    let c_theta = (1.0 - s_theta * s_theta).sqrt();
    from_polar(c_theta, phi)
}

// PDF = (n + 1) / (2 * pi) * cos(theta)^n sin(theta)
pub fn phong_sampling(u: f32, v: f32, n : f32) -> Sample3D
{
    // theta = acos((1-v)^(1/n+1))
    let phi = u * TAU;
    let c_theta = (1.0 - v).powf(1.0 / (n + 1.0));
    from_polar(c_theta, phi)
}

// PDF = D_ggx(theta) * cos(theta) * sin(theta)
pub fn ggx_sampling(u: f32, v: f32, a2 : f32) -> Sample3D
{
    // theta = acos(sqrt((1-v)/(v(a2-1)+1)))
    let phi = u * TAU;
    let c_theta = ((1.0 - v) / (v * (a2 - 1.0) + 1.0)).sqrt();
    from_polar(c_theta, phi)
}