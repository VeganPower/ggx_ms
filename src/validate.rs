use std::f32::consts::{PI, TAU};
mod sampling;

fn main()
{
    const BIT_DEPTH :usize = 10;
    const SAMPLE_COUNT :usize = 1 << BIT_DEPTH;

    // Integrate an hemisphere in the standard way
    // Expected result 2 PI
    let mut acc = 0.0;
    let d_theta = (PI * 0.5) / (SAMPLE_COUNT as f32);
    for j in 0..SAMPLE_COUNT
    {
        let theta = (j as f32 + 0.5) * d_theta;
        acc += theta.sin() * d_theta;
    }
    acc *= TAU;
    println!("Expected {}; Result  {}", TAU, acc);

    const K :usize = 64;
    let d_t = 1.0 / K as f32;
    let d_p = 1.0 / K as f32;
    let dx = 1.0 / (K * K) as f32;

    // This compute int_omega sin(theta)*dw
    acc = 0.0;
    for i in 0..K
    {
        let x = (i as f32 + 0.5) * d_p;
        for j in 0..K
        {
            let y = (j as f32 + 0.5) * d_t;
            let s = sampling::uniform_sampling(x, y);
            let cos_theta = s[2];
            let sin_theta = (1.0 - cos_theta*cos_theta).sqrt();
            let pdf = 1.0 / (PI*PI);
            acc += sin_theta / pdf;
        }
    }
    println!("Expected {}; Result  {}", TAU, acc * dx);

    // This compute int_omega sin(theta)*dw
    acc = 0.0;
    for i in 0..K
    {
        let x = (i as f32 + 0.5) * d_p;
        for j in 0..K
        {
            let y = (j as f32 + 0.5) * d_t;
            let s = sampling::hemisphere_sampling(x, y);
            let cos_theta = s[2];
            let _sin_theta = (1.0 - cos_theta*cos_theta).sqrt();
            let pdf = 1.0 / TAU;
            acc += 1.0 / pdf;
        }
    }
    println!("Expected {}; Result  {}", TAU, acc * dx);

    // This compute int_omega cos(theta)*sin(theta)*dw
    acc = 0.0;
    for i in 0..K
    {
        let x = (i as f32 + 0.5) * d_p;
        for j in 0..K
        {
            let y = (j as f32 + 0.5) * d_t;
            let s = sampling::cos_weighted_sampling(x, y);
            let cos_theta = s[2];
            let _sin_theta = (1.0 - cos_theta*cos_theta).sqrt();
            let pdf = 1.0 / PI;
            acc += 1.0 / pdf;
        }
    }
    println!("Expected {}; Result  {}", PI, acc * dx);
}