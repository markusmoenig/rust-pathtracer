use crate::prelude::*;

#[allow(unused)]
pub trait Scene : Sync + Send {

    fn new() -> Self where Self: Sized;

    fn closest_hit(&self, ray: &Ray, state: &mut State, light: &mut LightSampleRec) -> bool;

    fn any_hit(&self, ray: &Ray, max_dist: PTF) -> bool;

    fn light_at(&self, index: usize) -> &AnalyticalLight;
    fn number_of_lights(&self) -> usize;
}
