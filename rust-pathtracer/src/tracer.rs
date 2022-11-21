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

                    //ld += PTF3::new(0.1, 0.1, 0.1);

                    scatter_sample.f = self.disney_eval(state, -ray[1], &state.ffnormal, &light_sample.direction, &mut scatter_sample.pdf);

                    let mut mis_weight = 1.0;
                    if analytical.light.area > 0.0 {// No MIS for distant light
                        mis_weight = self.power_heuristic(&light_sample.pdf, &scatter_sample.pdf);
                    }

                    if scatter_sample.pdf > 0.0 {
                        ld += mis_weight * li.component_mul(&scatter_sample.f) / light_sample.pdf;
                    }
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

    #[inline(always)]
    fn power_heuristic(&self, a: &PTF, b: &PTF) -> PTF {
        let t = a * a;
        return t / (b * b + t);
    }

    #[inline(always)]
    fn mix_ptf(&self, a: &PTF, b: &PTF, v: &PTF) -> PTF {
        (1.0 - v) * a + b * v
    }

    // Disney

    fn luminance(&self, c: &PTF3) -> PTF {
        0.212671 * c.x + 0.715160 * c.y + 0.072169 * c.z
    }

    fn schlick_fresnel(&self, u: PTF) -> PTF {
        let m = (1.0 - u).clamp(0.0, 1.0);
        let m2 = m * m;
        m2 * m2 * m
    }

    fn dielectric_fresnel(&self, cos_theta_i: PTF, eta: PTF) -> PTF {
        let sin_theta_tsq = eta * eta * (1.0 - cos_theta_i * cos_theta_i);

        // Total internal reflection
        if sin_theta_tsq > 1.0 {
            return 1.0;
        }

        let cos_theta_t = (1.0 - sin_theta_tsq).max(0.0).sqrt();

        let rs = (eta * cos_theta_t - cos_theta_i) / (eta * cos_theta_t + cos_theta_i);
        let rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

        0.5 * (rs * rs + rp * rp)
    }

    fn get_spec_color(&self, material: &Material, eta: PTF, spec_col: &mut PTF3, sheen_col: &mut PTF3) {
        let lum = self.luminance(&material.base_color);
        let ctint = if lum > 0.0 { material.base_color / lum } else { PTF3::new(1.0, 1.0, 1.0) };
        let f0 = (1.0 - eta) / (1.0 + eta);
        *spec_col = glm::mix(&(f0 * f0 * glm::mix(&PTF3::new(1.0, 1.0, 1.0), &ctint, material.specular_tint)), &material.base_color, material.metallic);
        *sheen_col = glm::mix(&PTF3::new(1.0, 1.0, 1.0), &ctint, material.sheen_tint);
    }

    fn get_lobe_probabilities(&self, material: &Material, _eta: &PTF, spec_col: &PTF3, approx_fresnel: PTF, diffuse_wt: &mut PTF, spec_reflect_wt: &mut PTF,  spec_refract_wt: &mut PTF, clearcoat_wt: &mut PTF)
    {
        *diffuse_wt = self.luminance(&material.base_color) * (1.0 - material.metallic) * (1.0 - material.spec_trans);
        *spec_reflect_wt = self.luminance(&glm::mix(&spec_col, &PTF3::new(1.0, 1.0, 1.0), approx_fresnel));
        *spec_refract_wt = (1.0 - approx_fresnel) * (1.0 - material.metallic) * material.spec_trans * self.luminance(&material.base_color);
        *clearcoat_wt = 0.25 * material.clearcoat * (1.0 - material.metallic);
        let total_wt = *diffuse_wt + *spec_reflect_wt + *spec_refract_wt + *clearcoat_wt;

        *diffuse_wt /= total_wt;
        *spec_reflect_wt /= total_wt;
        *spec_refract_wt /= total_wt;
        *clearcoat_wt /= total_wt;
    }

    fn disney_fresnel(&self, material: &Material, eta: PTF, ldot_h: PTF, vdot_h: PTF) -> PTF {
        let metallic_fresnel = self.schlick_fresnel(ldot_h);
        let  dielectric_fresnel = self.dielectric_fresnel(vdot_h.abs(), eta);
        self.mix_ptf(&dielectric_fresnel, &metallic_fresnel, &material.metallic)
    }

    fn disney_eval(&self, state: &State, v: PTF3, n: &PTF3, l: &PTF3, bsdf_pdf: &mut PTF) -> PTF3 {
        *bsdf_pdf = 0.0;
        let f = PTF3::zeros();

        fn onb(n: &PTF3, t: &mut PTF3, b: &mut PTF3) {
            let up = if n.z.abs() < 0.999 { PTF3::new(0.0, 0.0, 1.0) } else { PTF3::new(1.0, 0.0, 0.0) };

            *t = glm::normalize(&glm::cross(&up, &n));
            *b = glm::cross(&n, &t);
        }

        fn to_local(x: &PTF3, y: &PTF3, z: &PTF3, v: &PTF3) -> PTF3 {
                PTF3::new(glm::dot(v, x), glm::dot(v, y), glm::dot(v, z))
        }

        let mut t = PTF3::zeros();
        let mut b = PTF3::zeros();

        onb(n, &mut t, &mut b);
        let v = to_local(&t, &b, n, &v); // NDotL = L.z; NDotV = V.z; NDotH = H.z
        let l = to_local(&t, &b, n, l);

        let mut h;
        if l.z > 0.0 {
            h = glm::normalize(&(l + v));
        } else {
            h = glm::normalize(&(l + v * state.eta));
        }

        if h.z < 0.0 {
            h = -h;
        }

        // Specular and sheen color
        let mut spec_col = PTF3::zeros();
        let mut sheen_col = PTF3::zeros();
        self.get_spec_color(&state.material, state.eta, &mut spec_col, &mut sheen_col);

        let mut diffuse_wt = 0.0; let mut spec_reflect_wt = 0.0; let mut spec_refract_wt = 0.0; let mut clearcoat_wt = 0.0;

        let fresnel = self.disney_fresnel(&state.material, state.eta, glm::dot(&l, &h), glm::dot(&v, &h));
        self.get_lobe_probabilities(&state.material, &state.eta, &spec_col, fresnel, &mut diffuse_wt, &mut spec_reflect_wt, &mut spec_refract_wt, &mut clearcoat_wt);


        f
    }


}