use crate::prelude::*;

#[inline(always)]
pub fn cross(a: &F3, b: &F3) -> F3 {
    a.cross(b)
}

#[inline(always)]
pub fn dot(a: &F3, b: &F3) -> F {
    a.dot(b)
}

#[inline(always)]
pub fn normalize(a: & F3) -> F3 {
    a.normalize()
}

#[inline(always)]
pub fn floor(a: & F3) -> F3 {
    a.floor()
}

#[inline(always)]
pub fn fract(a: & F3) -> F3 {
    a.fract()
}

#[inline(always)]
pub fn length(a: & F3) -> F {
    a.length()
}

#[inline(always)]
pub fn mix(a: &F3, b: &F3, v: &F) -> F3 {
    F3::new(
        (1.0 - v) * a.x + b.x * v,
        (1.0 - v) * a.y + b.y * v,
        (1.0 - v) * a.z + b.z * v
    )
}

#[inline(always)]
pub fn mix_f(a: &F, b: &F, v: &F) -> F {
    (1.0 - v) * a + b * v
}

#[inline(always)]
pub fn smoothstep(e0: F, e1: F, x: F) -> F {
    let t = ((x - e0) / (e1 - e0)). clamp(0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

#[inline(always)]
pub fn pow(a: &F3, exp: &F3) -> F3 {
    F3::new(
        a.x.powf(exp.x),
        a.y.powf(exp.y),
        a.z.powf(exp.z),
    )
}

#[inline(always)]
pub fn min(x: F, y: F) -> F {
    if x < y {
        x
    } else {
        y
    }
}

#[inline(always)]
pub fn max(x: F, y: F) -> F {
    if x > y {
        x
    } else {
        y
    }
}