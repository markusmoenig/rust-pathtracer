use rust_pathtracer::prelude::*;

pub struct AnalyticalScene {

}

// The Scene

impl Scene for AnalyticalScene {

    fn new() -> Self {

        Self {

        }
    }

    fn hit(&self, ray: &Ray) -> bool {

        if let Some(dist) = self.sphere(ray, PTF3::new(0.0, 0.0, 0.0), 1.0) {
            return true;
        }

        false
    }

}

// Analytical Intersections
// Based on https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

impl AnalyticalIntersections for AnalyticalScene {

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
