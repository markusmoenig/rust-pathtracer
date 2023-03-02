//! This is a port of the excellent [GLSL_Pathtracer](https://github.com/knightcrawler25/GLSL-PathTracer) to Rust utilizing an abstracted, trait based backend.
//! Perfect for rendering procedural content.
//!

pub type I = i64;
pub type F = f64;

const PI : F = std::f64::consts::PI;
const INV_PI : F = 1.0 / std::f64::consts::PI;
const TWO_PI : F = std::f64::consts::PI * 2.0;

/// Rays are stored in an array [origin, direction]
pub type Ray = [crate::fx::F3; 2];

pub mod fx;
pub mod math;
pub mod tracer;
pub mod camera;
pub mod scene;
pub mod material;
pub mod globals;
pub mod buffer;
pub mod light;

pub mod prelude {

    pub use crate::I;
    pub use crate::F;

    pub use crate::fx::F2;
    pub use crate::fx::F3;
    pub use crate::fx::B3;
    pub use crate::math::*;

    pub use crate::Ray;

    pub use crate::camera::Camera3D;
    pub use crate::camera::pinhole::Pinhole;

    pub use crate::light::AnalyticalLight;

    pub use crate::material::*;
    pub use crate::globals::*;

    pub use crate::buffer::ColorBuffer;

    pub use crate::scene::Scene;
    pub use crate::tracer::Tracer;
}
