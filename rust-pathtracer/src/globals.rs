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

    pub mat                 : Material,
    pub medium              : Medium,
}

impl State {
    pub fn new() -> Self {
        Self {
            depth           : 0,
            eta             : 0.0,

            hit_dist        : -1.0,

            fhp             : PTF3::new(0.0, 0.0, 0.0),
            normal          : PTF3::new(0.0, 0.0, 0.0),
            ffnormal        : PTF3::new(0.0, 0.0, 0.0),
            tangent         : PTF3::new(0.0, 0.0, 0.0),
            bitangent       : PTF3::new(0.0, 0.0, 0.0),

            is_emitter      : false,

            mat             : Material::new(PTF3::new(1.0, 1.0, 1.0)),
            medium          : Medium::new(),
        }
    }
}