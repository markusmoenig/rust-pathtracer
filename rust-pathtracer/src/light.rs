
use crate::prelude::*;

#[derive(PartialEq, Clone, Debug)]
/// An analytical light.
pub struct AnalyticalLight {
    pub light           : Light,
}

impl AnalyticalLight {

    /// Creates a spherical light based on its position, radius and emission.
    pub fn spherical(position: F3, radius: F, emission: F3) -> Self {

        let light = Light {
            light_type          : LightType::Spherical,
            position            : position,
            emission            : emission,
            u                   : F3::zeros(),
            v                   : F3::zeros(),
            radius              : radius,
            area                : 4.0 * crate::PI * radius * radius
        };

        Self {
            light,
        }
    }
}