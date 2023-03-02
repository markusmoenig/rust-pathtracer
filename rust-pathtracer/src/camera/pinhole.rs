
use crate::prelude::*;

/// A traditional pinhole camera with an origin, center and fov.
pub struct Pinhole {
    origin          : F3,
    center          : F3,

    fov             : F,
}

impl Camera3D for Pinhole {

    fn new() -> Self {

        let origin = F3::new(0.0, 0.0, 3.0);
        let center = F3::new(0.0, 0.0, 0.0);

        Self {
            origin,
            center,

            fov     : 80.0,
        }
    }

    fn set(&mut self, origin: F3, center: F3) {
        self.origin = origin;
        self.center = center;
    }

    fn set_fov(&mut self, fov: F) {
        self.fov = fov;
    }


    #[inline(always)]
    fn gen_ray(&self, p: F2, offset: F2, width: F, height: F) -> Ray {
        let ratio = width / height;

        let pixel_size = F2::new( 1.0 / width, 1.0 / height);

        let half_width = (self.fov.to_radians() * 0.5).tan();
        let half_height = half_width / ratio;

        let up_vector = F3::new(0.0, 1.0, 0.0);

        let w = (self.origin - self.center).normalize();
        let u = (&up_vector).cross(&w);
        let v = w.cross(&u);

        let lower_left = self.origin - u.mult_f(&half_width)  - v.mult_f(&half_height) - w;
        let horizontal = u.mult_f(&(half_width * 2.0));
        let vertical = v.mult_f(&(half_height * 2.0));

        let mut rd = lower_left - self.origin;
        rd += horizontal.mult_f(&(pixel_size.x * offset.x + p.x));
        rd += vertical.mult_f(&(pixel_size.y * offset.y + p.y));

        Ray::new(self.origin, rd.normalize())
    }
}
