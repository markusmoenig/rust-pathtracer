
use crate::prelude::*;

pub struct AnalyticalLight {
    pub light           : Light,
}

impl AnalyticalLight {

    pub fn spherical(position: PTF3, radius: PTF, emission: PTF3) -> Self {

        let light = Light {
            light_type          : LightType::Spherical,
            position            : position,
            emission            : emission,
            u                   : PTF3::zeros(),
            v                   : PTF3::zeros(),
            radius              : radius,
            area                : 10.0,
        };

        Self {
            light,
        }
    }
}