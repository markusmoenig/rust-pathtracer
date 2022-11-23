
use crate::prelude::*;

/// A traditional pinhole camera with an origin, center and fov.
pub struct Pinhole {
    origin          : PTF3,
    center          : PTF3,

    fov             : PTF,
}

impl Camera3D for Pinhole {

    fn new() -> Self {

        let origin = PTF3::new(0.0, 0.0, 3.0);
        let center = PTF3::new(0.0, 0.0, 0.0);

        Self {
            origin,
            center,

            fov     : 80.0,
        }
    }

    fn set(&mut self, origin: PTF3, center: PTF3) {
        self.origin = origin;
        self.center = center;
    }

    fn set_fov(&mut self, fov: PTF) {
        self.fov = fov;
    }


    #[inline(always)]
    fn gen_ray(&self, p: Vector2<PTF>, offset: PTF2, width: PTF, height: PTF) -> [Vector3<PTF>; 2] {
        let ratio = width / height;

        let pixel_size = PTF2::new( 1.0 / width, 1.0 / height);

        let half_width = (self.fov.to_radians() * 0.5).tan();
        let half_height = half_width / ratio;

        let up_vector = PTF3::new(0.0, 1.0, 0.0);

        let w = glm::normalize(&(self.origin - self.center));
        let u = glm::cross(&up_vector, &w);
        let v = glm::cross(&w, &u);

        let lower_left = self.origin - half_width * u - half_height * v - w;
        let horizontal = u * half_width * 2.0;
        let vertical = v * half_height * 2.0;

        let mut rd = lower_left - self.origin;
        rd += horizontal * (pixel_size.x * offset.x + p.x);
        rd += vertical * (pixel_size.y * offset.y + p.y);

        [self.origin, glm::normalize(&rd)]
    }
}
