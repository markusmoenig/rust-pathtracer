extern crate nalgebra_glm as glm;

pub type PTF = f64;
pub type PTF2 = glm::DVec2;
pub type PTF3 = glm::DVec3;
pub type PTF4 = glm::DVec4;
const PI : PTF = std::f64::consts::PI;
const INV_PI : PTF = 1.0 / std::f64::consts::PI;
const TWO_PI : PTF = std::f64::consts::PI * 2.0;

// pub type PTF = f32;
// pub type PTF2 = glm::Vec2;
// pub type PTF3 = glm::Vec3;
// pub type PTF4 = glm::Vec4;
// const PI : PTF = std::f32::consts::PI;
// const INV_PI : PTF = 1.0 / std::f32::consts::PI;
// const TWO_PI : PTF = std::f32::consts::PI * 2.0;

pub type Color = [PTF; 4];
pub type Ray = [PTF3; 2];

pub mod tracer;
pub mod camera;
pub mod scene;
pub mod material;
pub mod globals;
pub mod buffer;
pub mod light;

pub mod prelude {
    pub use crate::PTF;
    pub use crate::PTF2;
    pub use crate::PTF3;
    pub use crate::PTF4;
    pub use crate::Color;
    pub use crate::Ray;

    pub use nalgebra::*;

    pub use crate::camera::Camera3D;
    pub use crate::camera::pinhole::Pinhole;

    pub use crate::light::AnalyticalLight;

    pub use crate::material::*;
    pub use crate::globals::*;

    pub use crate::buffer::ColorBuffer;

    pub use crate::scene::Scene;
    pub use crate::tracer::Tracer;
}
