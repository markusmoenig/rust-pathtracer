use rayon::{slice::ParallelSliceMut, iter::{IndexedParallelIterator, ParallelIterator}};
use rand::{thread_rng, Rng, rngs::ThreadRng};
use crate::prelude::*;

pub struct Tracer {
    eps                 : PTF,

    camera              : Box<dyn Camera3D>,
    scene               : Box<dyn Scene>,
}

impl Tracer {

    pub fn new(scene: Box<dyn Scene>) -> Self {
        Self {
            eps             : 0.0001,

            camera      : Box::new(Pinhole::new()),
            scene,
        }
    }

    /// Render one frame and accumulate into the pixels buffer
    pub fn render(&mut self, buffer: &mut ColorBuffer) {

        const LINES: usize = 1;

        let width = buffer.width;
        let height = buffer.height as PTF;

        buffer.pixels
            .par_rchunks_exact_mut(width * LINES * 4)
            .enumerate()
            .for_each(|(j, line)| {
                for (i, pixel) in line.chunks_exact_mut(4).enumerate() {
                    let i = j * width * LINES + i;

                    let x = (i % width) as PTF;
                    let y = height - (i / width) as PTF;

                    let xx = (x as PTF) / width as PTF;
                    let yy = ((y) as PTF) / height as PTF;

                    // Camera

                    let mut rng = thread_rng();
                    let cam_offset = PTF2::new(rng.gen(), rng.gen());
                    let coord = Vector2::new(xx, 1.0 - yy);
                    let ray = self.camera.gen_ray(coord, 80.0, cam_offset, width as PTF, height);

                    // -

                    let mut radiance = PTF3::new(0.0, 0.0, 0.0);
                    let mut throughput = PTF3::new(1.0, 1.0, 1.0);
                    let mut state = State::new();
                    let mut light_sample = LightSampleRec::new();
                    let mut scatter_sample = ScatterSampleRec::new();

                    let alpha = 1.0;

                    for _i in 0..state.depth {

                        let hit = self.scene.closest_hit(&ray, &mut state, &mut light_sample);

                        if !hit {
                            break;
                        }

                        radiance += state.material.emission.component_mul(&throughput);

                        radiance += self.direct_light(&ray, &state, true, &mut rng).component_mul(&throughput);

                    }

                    let color = [radiance.x, radiance.y, radiance.z, alpha];

                    #[inline(always)]
                    pub fn mix_color(a: &[PTF], b: &[PTF], v: PTF) -> [PTF; 4] {
                        [   (1.0 - v) * a[0] + b[0] * v,
                            (1.0 - v) * a[1] + b[1] * v,
                            (1.0 - v) * a[2] + b[2] * v,
                            (1.0 - v) * a[3] + b[3] * v ]
                    }

                    let mix = mix_color(pixel, &color, 1.0 / (buffer.frames + 1) as PTF);

                    pixel.copy_from_slice(&mix);
                }
            });

        buffer.frames += 1;

    }

    /// Sample direct lights
    fn direct_light(&self, ray: &Ray, state: &State, is_surface: bool, rng: &mut ThreadRng) -> PTF3 {

        let mut ld = PTF3::new(0.0, 0.0, 0.0);
        let mut li;

        let scatter_pos = state.fhp + state.ffnormal * self.eps;
        let mut scatter_sample = ScatterSampleRec::new();

        let number_lights = self.scene.number_of_lights();

        if number_lights > 0 {
            let mut random : PTF = rng.gen();
            random *= self.scene.number_of_lights() as PTF;
            let index = random as usize;

            let analytical = self.scene.light_at(index);

            let mut light_sample = LightSampleRec::new();

            self.sample_light(&analytical.light, &scatter_pos, &mut light_sample, rng);
            li = light_sample.emission;

            if glm::dot(&light_sample.direction, &light_sample.normal) < 0.0 {// Required for quad lights with single sided emission

                let shadow_ray: Ray = [scatter_pos, light_sample.direction];

                let in_shawdow = self.scene.any_hit(&shadow_ray, light_sample.dist - self.eps);

                if in_shawdow == false {

                    ld += PTF3::new(0.5, 0.5, 0.5);

                }
            }
        }

        ld
    }

    /// Sample one light
    fn sample_light(&self, light: &Light, scatter_pos: &PTF3, light_sample: &mut LightSampleRec, rng: &mut ThreadRng) {

        match light.light_type {
            LightType::Spherical => {

                fn uniform_sample_hemisphere(r1: &PTF, r2: &PTF) -> PTF3 {
                    let r = 0.0.max(1.0 - r1 * r1).sqrt();
                    let phi = crate::TWO_PI * r2;
                    PTF3::new(r * phi.cos(), r * phi.sin(), *r1)
                }

                pub fn onb(n: PTF3, t: &mut PTF3, b: &mut PTF3) {
                    let up = if n.z.abs() < 0.999 { PTF3::new(0.0, 0.0, 1.0) } else { PTF3::new(1.0, 0.0, 0.0) };

                    *t = glm::normalize(&glm::cross(&up, &n));
                    *b = glm::cross(&n, &t);
                }

                let r1 : PTF = rng.gen();
                let r2 : PTF = rng.gen();

                let mut sphere_center_to_surface = scatter_pos - light.position;
                let dist_to_sphere_center = glm::length(&sphere_center_to_surface);
                let mut sampled_dir = uniform_sample_hemisphere(&r1, &r2);

                sphere_center_to_surface /= dist_to_sphere_center;

                let mut t = PTF3::new(0.0, 0.0, 0.0);
                let mut b = PTF3::new(0.0, 0.0, 0.0);

                onb(sphere_center_to_surface, &mut t, &mut b);
                sampled_dir = t * sampled_dir.x + b * sampled_dir.y + sphere_center_to_surface * sampled_dir.z;

                let light_surface_pos = light.position + sampled_dir * light.radius;

                light_sample.direction = light_surface_pos - scatter_pos;
                light_sample.dist = glm::length(&light_sample.direction);
                let dist_sq = light_sample.dist * light_sample.dist;

                light_sample.direction /= light_sample.dist;
                light_sample.normal = glm::normalize(&(light_surface_pos - light.position));
                light_sample.emission = light.emission * self.scene.number_of_lights() as PTF;
                light_sample.pdf = dist_sq / (light.area * 0.5 * glm::dot(&light_sample.normal, &light_sample.direction).abs());
            },
            _ => {},
        }

    }

}