use core::ops::{Add, Sub, Mul};

#[derive(Copy, Clone)]
pub struct Vec3 {
    pub x : f32,
    pub y : f32,
    pub z : f32,
}

impl Vec3
{
    pub fn new(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3 { x, y, z}
    }

    pub fn scale(&self, k : f32) -> Vec3 {
        Vec3 { x : self.x * k, y: self.y * k, z: self.z * k}
    }

    pub fn normalize(&self) -> Vec3
    {
        let l = dot(self, self).sqrt();
        self.scale(1.0 / l)
    }
}

pub fn cross(a : &Vec3, b : &Vec3) -> Vec3
{
    let cx = a.y * b.z - a.z * b.y;
    let cy = a.z * b.x - a.x * b.z;
    let cz = a.x * b.y - a.y * b.x;
    Vec3::new(cx, cy, cz)
}

pub fn dot(a : &Vec3, b : &Vec3) -> f32
{
    a.x * b.x + a.y * b.y + a.z * b.z
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {x: self.x + other.x, y: self.y + other.y, z: self.z + other.z}
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {x: self.x - other.x, y: self.y - other.y, z: self.z - other.z}
    }
}

impl Mul<f32> for Vec3 {
    type Output = Self;

    fn mul(self, k: f32) -> Self {
        Self {x: self.x * k, y: self.y * k, z: self.z * k}
    }
}

impl Mul<Vec3> for f32 {
    type Output = Vec3;

    fn mul(self, k: Vec3) -> Vec3 {
        Vec3 {x: k.x * self, y: k.y * self, z: k.z * self}
    }
}