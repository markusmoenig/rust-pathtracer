pub mod pinhole;

use crate::prelude::*;

#[allow(unused)]
pub trait Camera3D : Sync + Send {

    fn new() -> Self where Self: Sized;
    fn gen_ray(&self, p: Vector2<PTF>, fov: PTF, offset: PTF2, width: PTF, height: PTF) -> Ray;
}
