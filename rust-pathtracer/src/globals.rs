use crate::prelude::*;

// State

#[derive(PartialEq, Clone, Debug)]
pub struct State {
    pub depth               : u16,
    pub eta                 : F,

    pub hit_dist            : F,

    pub fhp                 : F3,
    pub normal              : F3,
    pub ffnormal            : F3,

    pub is_emitter          : bool,

    pub material            : Material,
    pub medium              : Medium,
}

impl State {
    pub fn new() -> Self {
        Self {
            depth           : 4,
            eta             : 0.0,

            hit_dist        : -1.0,

            fhp             : F3::zeros(),
            normal          : F3::zeros(),
            ffnormal        : F3::zeros(),

            is_emitter      : false,

            material        : Material::new(),
            medium          : Medium::new(),
        }
    }

    /// Calculate tangent and bitangent
    pub fn onb(&mut self, n: F3, t: &mut F3, b: &mut F3) {
        let up = if n.z.abs() < 0.999 { F3::new(0.0, 0.0, 1.0) } else { F3::new(1.0, 0.0, 0.0) };

        *t = up.cross(&n).normalize();
        *b = n.cross(&t);
    }

    /// State post-processing, called by the tracer after calling Scene::closest_hit()
    pub fn finalize(&mut self, ray: &Ray) {

        self.fhp = ray[0] + ray[1].mult_f(&self.hit_dist);

        if dot(&self.normal, &ray[1]) <= 0.0 {
            self.ffnormal = self.normal;
        } else {
            self.ffnormal = -self.normal;
        }

        self.material.finalize();
        self.eta = if dot(&ray[1], &self.normal) < 0.0 { 1.0 / self.material.ior } else { self.material.ior };
    }

}

// Light

#[derive(PartialEq, Clone, Debug)]
pub enum LightType {
    Rectangular,
    Spherical,
    Distant,
}

#[derive(PartialEq, Clone, Debug)]
pub struct Light {
    pub light_type          : LightType,
    pub position            : F3,
    pub emission            : F3,
    pub u                   : F3,
    pub v                   : F3,
    pub radius              : F,
    pub area                : F,
}

// ScatterSampleRec

#[derive(PartialEq, Clone, Debug)]
pub struct ScatterSampleRec {

    pub l                   : F3,
    pub f                   : F3,
    pub pdf                 : F,
}

impl ScatterSampleRec {
    pub fn new() -> Self {
        Self {
            l               : F3::zeros(),
            f               : F3::zeros(),
            pdf             : 0.0,
        }
    }
}

// LightSampleRec

#[derive(PartialEq, Clone, Debug)]
pub struct LightSampleRec {

    pub normal              : F3,
    pub emission            : F3,
    pub direction           : F3,

    pub dist                : F,
    pub pdf                 : F,
}

impl LightSampleRec {
    pub fn new() -> Self {
        Self {
            normal          : F3::zeros(),
            emission        : F3::zeros(),
            direction       : F3::zeros(),

            dist            : 0.0,
            pdf             : 0.0,
        }
    }
}

