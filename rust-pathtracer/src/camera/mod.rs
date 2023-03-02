pub mod pinhole;

use crate::prelude::*;

/// Trait for abstracting cameras.
#[allow(unused)]
pub trait Camera3D : Sync + Send {

    fn new() -> Self where Self: Sized;

    /// Set the origin and center of the camera.
    fn set(&mut self, origin: F3, center: F3);
    /// Set the fov of the camera.
    fn set_fov(&mut self, fov: F);

    /// Generate a ray.
    fn gen_ray(&self, p: F2, offset: F2, width: F, height: F) -> Ray;
}
