use crate::prelude::*;

// State

#[derive(PartialEq, Clone, Debug)]
pub struct State {
    pub depth               : u16,
    pub eta                 : PTF,

    pub hit_dist            : PTF,

    pub fhp                 : PTF3,
    pub normal              : PTF3,
    pub ffnormal            : PTF3,
    pub tangent             : PTF3,
    pub bitangent           : PTF3,

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

            fhp             : PTF3::new(0.0, 0.0, 0.0),
            normal          : PTF3::new(0.0, 0.0, 0.0),
            ffnormal        : PTF3::new(0.0, 0.0, 0.0),
            tangent         : PTF3::new(0.0, 0.0, 0.0),
            bitangent       : PTF3::new(0.0, 0.0, 0.0),

            is_emitter      : false,

            material        : Material::new(PTF3::new(1.0, 1.0, 1.0)),
            medium          : Medium::new(),
        }
    }

    /// Sets the hit distance and the hit normal
    pub fn set_distance_and_normal(&mut self, ray: &Ray, hit_dist: PTF, normal: PTF3) {
        self.hit_dist = hit_dist;
        self.fhp = ray[0] + ray[1] * hit_dist;

        self.normal = normal.clone();

        if glm::dot(&normal, &ray[1]) <= 0.0 {
            self.ffnormal = normal;
        } else {
            self.ffnormal = -normal;
        }
    }

    /// Calculate tangent and bitangent
    pub fn onb(&mut self, n: PTF3, t: &mut PTF3, b: &mut PTF3) {
        let up = if n.z.abs() < 0.999 { PTF3::new(0.0, 0.0, 1.0) } else { PTF3::new(1.0, 0.0, 0.0) };

        *t = glm::normalize(&glm::cross(&up, &n));
        *b = glm::cross(&n, &t);
    }

    /// State post-processing, called by the tracer after calling Scene::closest_hit()
    pub fn finalize(&mut self, ray: &Ray) {
        let mut t = PTF3::new(0.0, 0.0, 0.0);
        let mut b = PTF3::new(0.0, 0.0, 0.0);
        self.onb(self.normal, &mut t, &mut b);
        self.tangent = t;
        self.bitangent = b;
        self.material.finalize();
        self.eta = if glm::dot(&ray[1], &self.normal) < 0.0 { 1.0 / self.material.ior } else { self.material.ior };
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
    pub position            : PTF3,
    pub emission            : PTF3,
    pub u                   : PTF3,
    pub v                   : PTF3,
    pub radius              : PTF,
    pub area                : PTF,
}

// ScatterSampleRec

#[derive(PartialEq, Clone, Debug)]
pub struct ScatterSampleRec {

    pub l                   : PTF3,
    pub f                   : PTF3,
    pub pdf                 : PTF,
}

impl ScatterSampleRec {
    pub fn new() -> Self {
        Self {
            l               : PTF3::new(0.0, 0.0, 0.0),
            f               : PTF3::new(0.0, 0.0, 0.0),
            pdf             : 0.0,
        }
    }
}

// LightSampleRec

#[derive(PartialEq, Clone, Debug)]
pub struct LightSampleRec {

    pub normal              : PTF3,
    pub emission            : PTF3,
    pub direction           : PTF3,

    pub dist                : PTF,
    pub pdf                 : PTF,
}

impl LightSampleRec {
    pub fn new() -> Self {
        Self {
            normal          : PTF3::new(0.0, 0.0, 0.0),
            emission        : PTF3::new(0.0, 0.0, 0.0),
            direction       : PTF3::new(0.0, 0.0, 0.0),

            dist            : 0.0,
            pdf             : 0.0,
        }
    }
}

