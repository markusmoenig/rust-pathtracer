use rust_pathtracer::prelude::*;
pub use nalgebra::*;

pub struct AnalyticalScene {
    lights              : Vec<AnalyticalLight>,
    pinhole             : Box<dyn Camera3D>,
}

// The Scene

impl Scene for AnalyticalScene {

    fn new() -> Self {

        let lights = vec![AnalyticalLight::spherical(PTF3::new(3.0, 2.0, 2.0), 1.0, PTF3::new(3.0, 3.0, 3.0))];

        Self {
            lights,
            pinhole     : Box::new(Pinhole::new()),
        }
    }

    fn camera(&self) -> &Box<dyn Camera3D> {
        &self.pinhole
    }

    fn background(&self, ray: &Ray) -> PTF3 {
        // Taken from https://raytracing.github.io/books/RayTracingInOneWeekend.html, a source of great knowledge
        let t = 0.5 * (ray[1].y + 1.0);
        self.to_linear((1.0 - t) * PTF3::new(1.0, 1.0, 1.0) + t * PTF3::new(0.5, 0.7, 1.0)) * 0.1
    }


    /// The closest hit, includes light sources.
    fn closest_hit(&self, ray: &Ray, state: &mut State, light: &mut LightSampleRec) -> bool {

        let mut dist = PTF::MAX;
        let mut hit = false;

        let mut center = PTF3::new(-1.1, 0.0, 0.0);

        if let Some(d) = self.sphere(ray, center, 1.0) {

            let hp = ray[0] + ray[1] * d;
            let normal = glm::normalize(&(center - hp));

            state.hit_dist = d;
            state.normal = normal;

            // state.material.base_color = PTF3::new(1.0,0.4, 0.0);
            // state.material.clearcoat = 1.0;
            // state.material.clearcoat_gloss = 1.0;
            //state.material.roughness = 1.0;

            state.material.base_color = PTF3::new(5.0, 5.0, 5.0);//PTF3::new(0.815, 0.418501512, 0.00180012);
            state.material.roughness = 0.05;
            state.material.metallic = 1.0;

            // state.material.base_color = PTF3::new(1.0,1.0, 1.0);
            // state.material.spec_trans = 1.0;
            // state.material.roughness = 0.01;
            // state.material.ior = 1.45;

            hit = true;
            dist = d;
        }

        center = PTF3::new(1.1, 0.0, 0.0);

        if let Some(d) = self.sphere(ray, center, 1.0) {

            if d < dist {

                let hp = ray[0] + ray[1] * d;
                let normal = glm::normalize(&(center - hp));

                state.hit_dist = d;
                state.normal = normal;

                state.material.base_color = PTF3::new(1.0,0.186, 0.0);
                state.material.clearcoat = 1.0;
                state.material.clearcoat_gloss = 1.0;
                state.material.roughness = 0.6;

                // state.material.base_color = PTF3::new(0.9,0.9, 0.9);
                // state.material.roughness = 0.2;
                // state.material.metallic = 1.0;

                // state.material.base_color = PTF3::new(1.0,1.0, 1.0);
                // state.material.spec_trans = 1.0;
                // state.material.roughness = 0.01;
                // state.material.ior = 1.45;

                hit = true;
                dist = d;
            }
        }

        if let Some(d) = self.plane(ray) {

            if d < dist {
                state.hit_dist = d;
                state.normal = PTF3::new(0.0, 1.0, 0.0);

                fn checker(x: PTF, y: PTF) -> PTF {
                    let x1 = x.floor() % 2.0;
                    let y1 = y.floor() % 2.0;
                    if (x1 + y1) % 2.0 < 1.0 { 0.25 } else { 0.1 }
                }

                let c = checker(ray[1].x / ray[1].y * 0.5 + 100.0, ray[1].z / ray[1].y * 0.5 + 100.0);

                state.material.base_color = PTF3::new(c, c, c);
                state.material.roughness = 1.0;

                hit = true;
            }
        }

        hit
    }

    /// Any hit
    fn any_hit(&self, ray: &Ray, _max_dist: PTF) -> bool {

        if let Some(_d) = self.sphere(ray, PTF3::new(-1.1, 0.0, 0.0), 1.0) {
            return true;
        }

        if let Some(_d) = self.sphere(ray, PTF3::new(1.1, 0.0, 0.0), 1.0) {
            return true;
        }

        if let Some(_d) = self.plane(ray) {
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

    // Ray plane intersection
    fn plane(&self, ray: &Ray) -> Option<PTF> {
        let normal = PTF3::new(0.0, 1.0, 0.0);
        let denom = glm::dot(&normal, &ray[1]);

        if denom.abs() > 0.0001 {
            let t = glm::dot(&(PTF3::new(0.0, -1.0, 0.0) - ray[0]), &normal) / denom;
            if t >= 0.0 {
                return Some(t);
            }
        }
        None
    }
}

#[allow(unused)]
pub trait AnalyticalIntersections : Sync + Send {

    fn sphere(&self, ray: &Ray, center: PTF3, radius: PTF) -> Option<PTF>;
    fn plane(&self, ray: &Ray) -> Option<PTF>;

}
