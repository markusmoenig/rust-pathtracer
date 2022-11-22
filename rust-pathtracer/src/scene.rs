use crate::prelude::*;

#[allow(unused)]
pub trait Scene : Sync + Send {

    fn new() -> Self where Self: Sized;

    fn background(&self, ray: &Ray) -> PTF3;

    fn closest_hit(&self, ray: &Ray, state: &mut State, light: &mut LightSampleRec) -> bool;
    fn any_hit(&self, ray: &Ray, max_dist: PTF) -> bool;

    fn camera(&self) -> &Box<dyn Camera3D>;

    fn light_at(&self, index: usize) -> &AnalyticalLight;
    fn number_of_lights(&self) -> usize;

    fn to_linear(&self, c: PTF3) -> PTF3 {
        PTF3::new(c.x.powf(2.2), c.y.powf(2.2), c.z.powf(2.2))
    }
}
