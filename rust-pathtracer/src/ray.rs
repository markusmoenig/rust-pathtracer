use crate::prelude::*;
use rhai::Engine;

/// Ray
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Ray {
    pub origin          : F3,
    pub direction       : F3,
}

impl Ray {
    pub fn new(origin: F3, direction: F3) -> Self {

        Self {
            origin,
            direction
        }
    }

    pub fn at(&self, dist: &F) -> F3 {
        self.origin + *dist * self.direction
    }

    pub fn get_origin(&mut self) -> F3 {
        self.origin
    }

    pub fn get_direction(&mut self) -> F3 {
        self.direction
    }

    /// Register to the engine
    pub fn register(engine: &mut Engine) {
        engine.register_type_with_name::<Ray>("Ray")
            .register_get("origin", Ray::get_origin)
            .register_get("direction", Ray::get_direction);
    }
}