use crate::prelude::*;

pub fn cross(a: &F3, b: &F3) -> F3 {
    a.cross(b)
}

pub fn dot(a: &F3, b: &F3) -> F {
    a.dot(b)
}

pub fn normalize(a: & F3) -> F3 {
    a.normalize()
}

pub fn length(a: & F3) -> F {
    a.length()
}

pub fn mix(a: &F3, b: &F3, v: &F) -> F3 {
    F3::new(
        (1.0 - v) * a.x + b.x * v,
        (1.0 - v) * a.y + b.y * v,
        (1.0 - v) * a.z + b.z * v
    )
}

pub fn pow(a: &F3, exp: &F3) -> F3 {
    F3::new(
        a.x.powf(exp.x),
        a.y.powf(exp.y),
        a.z.powf(exp.z),
    )
}