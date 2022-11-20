extern crate nalgebra_glm as glm;

pub type PTF = f64;
pub type PTF2 = glm::DVec2;
pub type PTF3 = glm::DVec3;
pub type PTF4 = glm::DVec4;

// pub type PTF = f32;
// pub type PTF2 = glm::Vec2;
// pub type PTF3 = glm::Vec3;
// pub type PTF4 = glm::Vec4;

pub type Color = [PTF; 4];
pub type Ray = [PTF3; 2];

pub mod tracer;
pub mod camera;
pub mod scene;
pub mod material;
pub mod globals;

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

    pub use crate::material::*;
    pub use crate::globals::*;

    pub use crate::scene::Scene;
    pub use crate::tracer::Tracer;
}
