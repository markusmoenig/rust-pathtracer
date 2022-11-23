pub mod pinhole;

use crate::prelude::*;

/// Trait for abstracting cameras.
#[allow(unused)]
pub trait Camera3D : Sync + Send {

    fn new() -> Self where Self: Sized;

    /// Set the origin and center of the camera.
    fn set(&mut self, origin: PTF3, center: PTF3);
    /// Set the fov of the camera.
    fn set_fov(&mut self, fov: PTF);

    /// Generate a ray.
    fn gen_ray(&self, p: Vector2<PTF>, offset: PTF2, width: PTF, height: PTF) -> Ray;
}
