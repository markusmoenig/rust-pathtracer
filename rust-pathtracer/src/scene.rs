use crate::prelude::*;

/// A trait based scene abstraction.
#[allow(unused)]
pub trait Scene : Sync + Send {

    fn new() -> Self where Self: Sized;

    /// Background color for the given ray
    fn background(&self, ray: &Ray) -> F3;

    /// Closest hit should return the state.hit_dist, state.normal and fill out the state.material as needed
    fn closest_hit(&self, ray: &Ray, state: &mut State, light: &mut LightSampleRec) -> bool;

    /// Used for shadow rays.
    fn any_hit(&self, ray: &Ray, max_dist: F) -> bool;

    /// Return the camera for the scene
    fn camera(&self) -> &Box<dyn Camera3D>;

    /// Return the number of lights in the scene
    fn number_of_lights(&self) -> usize;

    /// Return a reference for the light at the given index
    fn light_at(&self, index: usize) -> &AnalyticalLight;

    /// The recursion depth for the path tracer
    fn recursion_depth(&self) -> u16 {
        4
    }

    fn to_linear(&self, c: F3) -> F3 {
        F3::new(c.x.powf(2.2), c.y.powf(2.2), c.z.powf(2.2))
    }
}
