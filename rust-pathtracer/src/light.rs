
use crate::prelude::*;

/// An analytical light.
pub struct AnalyticalLight {
    pub light           : Light,
}

impl AnalyticalLight {

    /// Creates a spherical light based on its position, radius and emission.
    pub fn spherical(position: PTF3, radius: PTF, emission: PTF3) -> Self {

        let light = Light {
            light_type          : LightType::Spherical,
            position            : position,
            emission            : emission,
            u                   : PTF3::zeros(),
            v                   : PTF3::zeros(),
            radius              : radius,
            area                : 4.0 * crate::PI * radius * radius
        };

        Self {
            light,
        }
    }
}