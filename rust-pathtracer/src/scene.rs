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

    fn sample_lights(&self, ray: &Ray, state: &mut State, light_sample: &mut LightSampleRec, lights: &Vec<AnalyticalLight>) -> bool {

        // Based on https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
        fn sphere(ray: &Ray, center: F3, radius: F) -> Option<F> {
            let l = center - ray.origin;
            let tca = l.dot(&ray.direction);
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

        let mut hit = false;
        let mut dist = state.hit_dist;

        for light in lights {
            if light.light.light_type == LightType::Spherical {
                if let Some(d) = sphere(ray, light.light.position, light.light.radius) {
                    if d < dist {
                        dist = d;
                        let hit_point = ray.at(&d);
                        let cos_theta = dot(&-ray.direction, &normalize(&(hit_point - light.light.position)));
                        light_sample.pdf = (dist * dist) / (light.light.area * cos_theta * 0.5);
                        light_sample.emission = light.light.emission;
                        state.is_emitter = true;
                        state.hit_dist = d;
                        hit = true;
                    }
                }
            }
        }

        hit
    }
}
