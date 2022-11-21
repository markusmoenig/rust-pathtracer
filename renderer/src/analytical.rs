use rust_pathtracer::prelude::*;
pub use nalgebra::*;

pub struct AnalyticalScene {
    lights              : Vec<AnalyticalLight>,
}

// The Scene

impl Scene for AnalyticalScene {

    fn new() -> Self {

        let lights = vec![AnalyticalLight::spherical(PTF3::new(3.0, 2.0, 2.0), 1.5, PTF3::new(1.0, 1.0, 1.0))];

        Self {
            lights,
        }
    }

    /// The closest hit, includes light sources.
    fn closest_hit(&self, ray: &Ray, state: &mut State, light: &mut LightSampleRec) -> bool {

        let center = PTF3::new(0.0, 0.0, 0.0);

        if let Some(dist) = self.sphere(ray, PTF3::new(0.0, 0.0, 0.0), 1.0) {

            let hp = ray[0] + ray[1] * dist;
            let normal = glm::normalize(&(center - hp));

            state.set_distance_and_normal(ray, dist, normal);

            return true;
        }

        false
    }

    /// Any hit
    fn any_hit(&self, ray: &Ray, max_dist: PTF) -> bool {

        if let Some(dist) = self.sphere(ray, PTF3::new(0.0, 0.0, 0.0), 1.0) {
            return true;
        }

        false
    }

    /// Returns the light at the given index
    fn light_at(&self, index: usize) -> &AnalyticalLight {
        &self.lights[index]
    }

    fn number_of_lights(&self) -> usize {
        self.lights.len()
    }

}

// Analytical Intersections

impl AnalyticalIntersections for AnalyticalScene {

    // Based on https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    fn sphere(&self, ray: &Ray, center: PTF3, radius: PTF) -> Option<PTF> {
        let l = center - ray[0];
        let tca = l.dot(&ray[1]);
        let d2 = l.dot(&l) - tca * tca;
        let radius2 = radius * radius;
        if d2 > radius2 {
            return None;
        }
        let thc = (radius2 - d2).sqrt();
        let mut t0 = tca - thc;
        let mut t1 = tca + thc;

        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }

        if t0 < 0.0 {
            t0 = t1;
            if t0 < 0.0 {
                return None;
            }
        }

        Some(t0)
   }
}

#[allow(unused)]
pub trait AnalyticalIntersections : Sync + Send {

    fn sphere(&self, ray: &Ray, center: PTF3, radius: PTF) -> Option<PTF>;
}
