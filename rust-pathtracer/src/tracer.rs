use rayon::{slice::ParallelSliceMut, iter::{IndexedParallelIterator, ParallelIterator}};
use rand::{thread_rng, Rng, rngs::ThreadRng};
use crate::prelude::*;

pub struct Tracer {
    eps                 : F,

    scene               : Box<dyn Scene>,
}

impl Tracer {

    pub fn new(scene: Box<dyn Scene>) -> Self {
        Self {

            eps         : 0.0001,
            scene,
        }
    }

    /// Render one frame and accumulate into the pixels buffer
    pub fn render(&mut self, buffer: &mut ColorBuffer) {

        const LINES: usize = 1;

        let width = buffer.width;
        let height = buffer.height as F;

        buffer.pixels
            .par_rchunks_exact_mut(width * LINES * 4)
            .enumerate()
            .for_each(|(j, line)| {
                for (i, pixel) in line.chunks_exact_mut(4).enumerate() {
                    let i = j * width * LINES + i;

                    let x = (i % width) as F;
                    let y = height - (i / width) as F;

                    let xx = (x as F) / width as F;
                    let yy = ((y) as F) / height as F;

                    // Camera

                    let mut rng = thread_rng();
                    let cam_offset = F2::new(rng.gen(), rng.gen());
                    let coord = F2::new(xx, 1.0 - yy);
                    let mut ray = self.scene.camera().gen_ray(coord, cam_offset, width as F, height);

                    // -

                    let mut radiance = F3::new(0.0, 0.0, 0.0);
                    let mut throughput = F3::new(1.0, 1.0, 1.0);
                    let mut state = State::new();
                    let mut light_sample = LightSampleRec::new();
                    let mut scatter_sample = ScatterSampleRec::new();

                    state.depth = self.scene.recursion_depth();

                    let alpha = 1.0;

                    for _i in 0..state.depth {

                        state.material = Material::new();
                        let hit = self.scene.closest_hit(&ray, &mut state, &mut light_sample);

                        if !hit {
                            radiance += self.scene.background(&ray) * throughput;
                            break;
                        }

                        // State and material post-processing
                        state.finalize(&ray);

                        radiance += state.material.emission * throughput;

                        radiance += self.direct_light(&ray, &state, true, &mut rng) * throughput;

                        // Sample BSDF for color and outgoing direction
                        scatter_sample.f = self.disney_sample(&state, -ray.direction, &state.ffnormal, &mut scatter_sample.l, &mut scatter_sample.pdf, &mut rng);
                        if scatter_sample.pdf > 0.0 {
                           throughput = throughput * (scatter_sample.f / F3::new_x(scatter_sample.pdf));
                        } else {
                           break;
                        }

                        // Move ray origin to hit point and set direction for next bounce
                        ray.direction = scatter_sample.l;
                        ray.origin = state.fhp + self.eps * ray.direction;

                    }

                    let color = [radiance.x, radiance.y, radiance.z, alpha];

                    #[inline(always)]
                    pub fn mix_color(a: &[F], b: &[F], v: F) -> [F; 4] {
                        [   (1.0 - v) * a[0] + b[0] * v,
                            (1.0 - v) * a[1] + b[1] * v,
                            (1.0 - v) * a[2] + b[2] * v,
                            (1.0 - v) * a[3] + b[3] * v ]
                    }

                    let mix = mix_color(pixel, &color, 1.0 / (buffer.frames + 1) as F);

                    pixel.copy_from_slice(&mix);
                }
            });

        buffer.frames += 1;

    }

    /// Sample direct lights
    fn direct_light(&self, ray: &Ray, state: &State, _is_surface: bool, rng: &mut ThreadRng) -> F3 {

        let mut ld = F3::zeros();
        let li;

        let scatter_pos = state.fhp + self.eps * state.ffnormal;
        let mut scatter_sample = ScatterSampleRec::new();

        let number_lights = self.scene.number_of_lights();

        if number_lights > 0 {
            let mut random : F = rng.gen();
            random *= self.scene.number_of_lights() as F;
            let index = random as usize;

            let analytical = self.scene.light_at(index);

            let mut light_sample = LightSampleRec::new();

            self.sample_light(&analytical.light, &scatter_pos, &mut light_sample, rng);
            li = light_sample.emission;

            if dot(&light_sample.direction, &light_sample.normal) < 0.0 {// Required for quad lights with single sided emission

                let shadow_ray: Ray = Ray::new(scatter_pos, light_sample.direction);

                let in_shawdow = self.scene.any_hit(&shadow_ray, light_sample.dist - self.eps);

                if in_shawdow == false {
                    scatter_sample.f = self.disney_eval(state, -ray.direction, &state.ffnormal, &light_sample.direction, &mut scatter_sample.pdf);

                    let mut mis_weight = 1.0;
                    if analytical.light.area > 0.0 {// No MIS for distant light
                        mis_weight = self.power_heuristic(&light_sample.pdf, &scatter_sample.pdf);
                    }

                    if scatter_sample.pdf > 0.0 {
                        ld += mis_weight * li * (scatter_sample.f / F3::new_x(light_sample.pdf));
                    }
                }
            }
        }

        ld
    }

    /// Sample one light
    fn sample_light(&self, light: &Light, scatter_pos: &F3, light_sample: &mut LightSampleRec, rng: &mut ThreadRng) {

        match light.light_type {
            LightType::Spherical => {

                fn uniform_sample_hemisphere(r1: &F, r2: &F) -> F3 {
                    let r = 0.0_f64.max(1.0 - r1 * r1).sqrt();
                    let phi = crate::TWO_PI * r2;
                    F3::new(r * phi.cos(), r * phi.sin(), *r1)
                }

                pub fn onb(n: F3, t: &mut F3, b: &mut F3) {
                    let up = if n.z.abs() < 0.999 { F3::new(0.0, 0.0, 1.0) } else { F3::new(1.0, 0.0, 0.0) };

                    *t = normalize(&cross(&up, &n));
                    *b = cross(&n, &t);
                }

                let r1 : F = rng.gen();
                let r2 : F = rng.gen();

                let mut sphere_center_to_surface = *scatter_pos - light.position;
                let dist_to_sphere_center = length(&sphere_center_to_surface);
                let mut sampled_dir = uniform_sample_hemisphere(&r1, &r2);

                sphere_center_to_surface /= F3::new_x(dist_to_sphere_center);

                let mut t = F3::new(0.0, 0.0, 0.0);
                let mut b = F3::new(0.0, 0.0, 0.0);

                onb(sphere_center_to_surface, &mut t, &mut b);
                sampled_dir = sampled_dir.x * t + sampled_dir.y * b + sampled_dir.z * sphere_center_to_surface;

                let light_surface_pos = light.position + light.radius * sampled_dir;

                light_sample.direction = light_surface_pos - *scatter_pos;
                light_sample.dist = length(&light_sample.direction);
                let dist_sq = light_sample.dist * light_sample.dist;

                light_sample.direction /= F3::new_x(light_sample.dist);
                light_sample.normal = normalize(&(light_surface_pos - light.position));
                light_sample.emission = self.scene.number_of_lights() as F * light.emission;
                light_sample.pdf = dist_sq / (light.area * 0.5 * dot(&light_sample.normal, &light_sample.direction).abs());
            },
            _ => {},
        }

    }

    #[inline(always)]
    fn power_heuristic(&self, a: &F, b: &F) -> F {
        let t = a * a;
        return t / (b * b + t);
    }

    #[inline(always)]
    fn mix_ptf(&self, a: &F, b: &F, v: F) -> F {
        (1.0 - v) * a + b * v
    }

    fn gtr1(&self, ndoth: &F, a: F) -> F {
        if a >= 1.0 {
            return crate::INV_PI;
        }
        let a2 = a * a;
        let t = 1.0 + (a2 - 1.0) * ndoth * ndoth;
        return (a2 - 1.0) / (crate::PI * (a2).log2() * t);
    }

    fn sample_gtr1(&self, rgh: F, r1: F, _r2: F) -> F3 {
        let a = 0.001_f64.max(rgh);
        let a2 = a * a;

        let phi = r1 * crate::TWO_PI;

        let cos_theta = ((1.0 - a2.powf(1.0 - r1)) / (1.0 - a2)).sqrt();
        let sin_theta = (1.0 - (cos_theta * cos_theta)).sqrt().clamp(0.0, 1.0);
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        F3::new(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta)
    }

    fn sample_ggxvndf(&self, v: &F3, ax: F, ay: F, r1: F, r2: F) -> F3
    {
        let vh = normalize(&F3::new(ax * v.x, ay * v.y, v.z));

        let lensq = vh.x * vh.x + vh.y * vh.y;
        let t_1 = if lensq > 0.0 { F3::new(-vh.y, vh.x, 0.0).mult_f(&(1.0 / lensq.sqrt())) } else { F3::new(1.0, 0.0, 0.0) };
        let t_2 = cross(&vh, &t_1);

        let r = r1.sqrt();
        let phi = 2.0 * crate::PI * r2;
        let t1 = r * phi.cos();
        let mut t2 = r * phi.sin();
        let s = 0.5 * (1.0 + vh.z);
        t2 = (1.0 - s) * (1.0 - t1 * t1).sqrt() + s * t2;

        let nh = t1 * t_1 + t2 * t_2 + (0.0_f64.max(1.0 - t1 * t1 - t2 * t2)).sqrt() * vh;

        normalize(&F3::new(ax * nh.x, ay * nh.y, 0.0_f64.max(nh.z)))
    }

    fn smithg(&self, ndotv: &F, alphag: F) -> F {
        let a = alphag * alphag;
        let b = ndotv * ndotv;
        (2.0 * ndotv) / (ndotv + (a + b - a * b).sqrt())
    }

    // Disney

    fn luminance(&self, c: &F3) -> F {
        0.212671 * c.x + 0.715160 * c.y + 0.072169 * c.z
    }

    fn schlick_fresnel(&self, u: F) -> F {
        let m = (1.0 - u).clamp(0.0, 1.0);
        let m2 = m * m;
        m2 * m2 * m
    }

    fn gtr2aniso(&self, ndoth: &F, hdotx: &F, hdoty: &F, ax: &F, ay: &F) -> F {
        let a = hdotx / ax;
        let b = hdoty / ay;
        let c = a * a + b * b + ndoth * ndoth;
        1.0 / (crate::PI * ax * ay * c * c)
    }

    fn smithganiso(&self, ndotv: &F, vdotx: &F, vdoty: &F, ax: &F, ay: &F) -> F {
        let a = vdotx * ax;
        let b = vdoty * ay;
        let c = ndotv;
        (2.0 * ndotv) / (ndotv + (a * a + b * b + c * c).sqrt())
    }

    fn dielectric_fresnel(&self, cos_theta_i: F, eta: F) -> F {
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

    fn cosine_sample_hemisphere(&self, r1: F, r2: F) -> F3
    {
        let mut dir = F3::zeros();
        let r = r1.sqrt();
        let phi = crate::TWO_PI * r2;
        dir.x = r * phi.cos();
        dir.y = r * phi.sin();
        dir.z = 0.0_f64.max(1.0 - dir.x * dir.x - dir.y * dir.y).sqrt();
        dir
    }

    fn get_spec_color(&self, material: &Material, eta: F, spec_col: &mut F3, sheen_col: &mut F3) {
        let lum = self.luminance(&material.rgb);
        let ctint = if lum > 0.0 { material.rgb / F3::new_x(lum) } else { F3::new(1.0, 1.0, 1.0) };
        let f0 = (1.0 - eta) / (1.0 + eta);
        *spec_col = mix(&(f0 * f0 * mix(&F3::new(1.0, 1.0, 1.0), &ctint, &material.specular_tint)), &material.rgb, &material.metallic);
        *sheen_col = mix(&F3::new(1.0, 1.0, 1.0), &ctint, &material.sheen_tint);
    }

    fn eval_diffuse(&self, material: &Material, c_sheen: &F3, v: &F3, l: &F3, h: &F3, pdf: &mut F) -> F3 {
        *pdf = 0.0;
        if l.z <= 0.0 {
            return F3::zeros();
        }

        // Diffuse
        let fl = self.schlick_fresnel(l.z);
        let fv = self.schlick_fresnel(v.z);
        let fh = self.schlick_fresnel(dot(&l, &h));
        let fd90 = 0.5 + 2.0 * dot(&l, &h) * dot(&l, &h) * material.roughness;
        let fd = self.mix_ptf(&1.0, &fd90, fl) * self.mix_ptf(&1.0, &fd90, fv);

        // Fake Subsurface TODO: Replace with volumetric scattering
        let fss90 = dot(&l, &h) * dot(&l, &h) * material.roughness;
        let fss = self.mix_ptf(&1.0, &fss90, fl) * self.mix_ptf(&1.0, &fss90, fv);
        let ss = 1.25 * (fss * (1.0 / (l.z + v.z) - 0.5) + 0.5);

        // Sheen
        let fsheen = fh * material.sheen * *c_sheen;

        *pdf = l.z * crate::INV_PI;
        (1.0 - material.metallic) * (1.0 - material.spec_trans) * (crate::INV_PI * self.mix_ptf(&fd, &ss, material.subsurface) * material.rgb + fsheen)
    }

    fn eval_spec_reflection(&self, material: &Material, eta: F, spec_col: &F3, v: &F3, l: &F3, h: &F3, pdf: &mut F) -> F3 {
        *pdf = 0.0;
        if l.z <= 0.0 {
            return F3::zeros();
        }

        let fm = self.disney_fresnel(material, eta, dot(&l, &h), dot(&v, &h));
        let f = mix(&spec_col, &F3::new(1.0, 1.0, 1.0), &fm);
        let d = self.gtr2aniso(&h.z, &h.x, &h.y, &material.ax, &material.ay);
        let g1 = self.smithganiso(&v.z.abs(), &v.x, &v.y, &material.ax, &material.ay);
        let g2 = g1 * self.smithganiso(&l.z.abs(), &l.x, &l.y, &material.ax, &material.ay);

        *pdf = g1 * d / (4.0 * v.z);
        d * g2 * f / F3::new_x(4.0 * l.z * v.z)
    }

    fn eval_spec_refraction(&self, material: &Material, eta: F, v: &F3, l: &F3, h: &F3, pdf: &mut F) -> F3 {
        *pdf = 0.0;
        if l.z >= 0.0 {
            return F3::zeros();
        }

        let f = self.dielectric_fresnel(dot(&v, &h).abs(), eta);
        let d = self.gtr2aniso(&h.z, &h.x, &h.y, &material.ax, &material.ay);
        let g1 = self.smithganiso(&v.z.abs(), &v.x, &v.y, &material.ax, &material.ay);
        let g2 = g1 * self.smithganiso(&l.z.abs(), &l.x, &l.y, &material.ax, &material.ay);
        let mut denom = dot(&l, &h) + dot(&v, &h) * eta;
        denom *= denom;
        let eta2 = eta * eta;
        let jacobian = dot(&l, &h).abs() / denom;

        *pdf = g1 * 0.0_f64.max(dot(&v, &h)) * d * jacobian / v.z;

        (1.0 - material.metallic) * material.spec_trans * (1.0 - f) * d * g2 * dot(&v, &h).abs() * jacobian * eta2 / (l.z * v.z).abs() * pow(&material.rgb, &F3::new(0.5, 0.5, 0.5))
    }

    fn eval_clearcoat(&self, material: &Material, v: &F3, l: &F3, h: &F3, pdf: &mut F) -> F3
    {
        *pdf = 0.0;
        if l.z <= 0.0 {
            return F3::zeros();
        }

        let fh = self.dielectric_fresnel(dot(&v, &h), 1.0 / 1.5);
        let f = self.mix_ptf(&0.04, &1.0, fh);
        let d = self.gtr1(&h.z, material.clearcoat_roughness);
        let g = self.smithg(&l.z, 0.25) * self.smithg(&v.z, 0.25);
        let jacobian = 1.0 / (4.0 * dot(&v, &h));

        *pdf = d * h.z * jacobian;
        material.clearcoat * f * d * g / (4.0 * l.z * v.z) * F3::new(0.25, 0.25, 0.25)
    }

    fn get_lobe_probabilities(&self, material: &Material, _eta: &F3, spec_col: &F3, approx_fresnel: F, diffuse_wt: &mut F, spec_reflect_wt: &mut F,  spec_refract_wt: &mut F, clearcoat_wt: &mut F)
    {
        *diffuse_wt = self.luminance(&material.rgb) * (1.0 - material.metallic) * (1.0 - material.spec_trans);
        *spec_reflect_wt = self.luminance(&mix(&spec_col, &F3::new(1.0, 1.0, 1.0), &approx_fresnel));
        *spec_refract_wt = (1.0 - approx_fresnel) * (1.0 - material.metallic) * material.spec_trans * self.luminance(&material.rgb);
        *clearcoat_wt = 0.25 * material.clearcoat * (1.0 - material.metallic);
        let total_wt = *diffuse_wt + *spec_reflect_wt + *spec_refract_wt + *clearcoat_wt;

        *diffuse_wt /= total_wt;
        *spec_reflect_wt /= total_wt;
        *spec_refract_wt /= total_wt;
        *clearcoat_wt /= total_wt;
    }

    fn disney_fresnel(&self, material: &Material, eta: F, ldot_h: F, vdot_h: F) -> F {
        let metallic_fresnel = self.schlick_fresnel(ldot_h);
        let  dielectric_fresnel = self.dielectric_fresnel(vdot_h.abs(), eta);
        self.mix_ptf(&dielectric_fresnel, &metallic_fresnel, material.metallic)
    }

    fn disney_sample(&self, state: &State, mut v: F3, n: &F3, l: &mut F3, pdf: &mut F, rng: &mut ThreadRng) -> F3 {

        *pdf = 0.0;
        let f;

        let mut r1 : F = rng.gen();
        let r2 : F = rng.gen();

        fn onb(n: &F3, t: &mut F3, b: &mut F3) {
            let up = if n.z.abs() < 0.999 { F3::new(0.0, 0.0, 1.0) } else { F3::new(1.0, 0.0, 0.0) };

            *t = normalize(&cross(&up, n));
            *b = cross(n, t);
        }

        fn to_local(x: &F3, y: &F3, z: &F3, v: &F3) -> F3 {
                F3::new(dot(v, x), dot(v, y), dot(v, z))
        }

        fn to_world(x: &F3, y: &F3, z: &F3, v: &F3) -> F3 {
            v.x * *x + v.y * *y + v.z * *z
        }

        fn reflect(i: F3, n: F3) -> F3 {
            i - F3::new(2.0, 2.0, 2.0) * n * F3::new_x(dot(&n, &i))
        }

        fn refract(i: F3, n: F3, eta: F) -> F3 {
            let k = 1.0 - eta * eta * (1.0 - dot(&n, &i) * dot(&n, &i));
            if k < 0.0 {
                F3::zeros()
            } else {
                eta * i - (eta * dot(&n, &i) + k.sqrt()) * n
            }
        }

        let mut t = F3::zeros();
        let mut b = F3::zeros();

        onb(&n, &mut t, &mut b);
        v = to_local(&t, &b, n, &v); // NDotL = L.z; NDotV = V.z; NDotH = H.z

        // Specular and sheen color
        let mut spec_col = F3::zeros();
        let mut sheen_col = F3::zeros();
        self.get_spec_color(&state.material, state.eta, &mut spec_col, &mut sheen_col);

        // Lobe weights
        let mut diffuse_wt = 0.0; let mut spec_reflect_wt = 0.0; let mut spec_refract_wt = 0.0; let mut clearcoat_wt = 0.0;

        let approx_fresnel = self.disney_fresnel(&state.material, state.eta, v.z, v.z);
        self.get_lobe_probabilities(&state.material, &F3::new_x(state.eta), &spec_col, approx_fresnel, &mut diffuse_wt, &mut spec_reflect_wt, &mut spec_refract_wt, &mut clearcoat_wt);

        // CDF for picking a lobe
        let mut cdf = [0.0, 0.0, 0.0, 0.0];
        cdf[0] = diffuse_wt;
        cdf[1] = cdf[0] + clearcoat_wt;
        cdf[2] = cdf[1] + spec_reflect_wt;
        cdf[3] = cdf[2] + spec_refract_wt;

        if r1 < cdf[0] { // Diffuse Reflection Lobe
            r1 /= cdf[0];
            *l = self.cosine_sample_hemisphere(r1, r2);

            let h = normalize(&(*l + v));
            f = self.eval_diffuse(&state.material, &sheen_col, &v, &l, &h, pdf);
            *pdf *= diffuse_wt;
        } else
        if r1 < cdf[1] {// Clearcoat Lobe
            r1 = (r1 - cdf[0]) / (cdf[1] - cdf[0]);

            let mut h = self.sample_gtr1(state.material.clearcoat_roughness, r1, r2);

            if h.z < 0.0 {
                h = -h;
            }

            *l = normalize(&reflect(-v, h));
            f = self.eval_clearcoat(&state.material, &v, l, &h, pdf);
            *pdf *= clearcoat_wt;
        } else  // Specular Reflection/Refraction Lobes
        {
            r1 = (r1 - cdf[1]) / (1.0 - cdf[1]);
            let mut h = self.sample_ggxvndf(&v, state.material.ax, state.material.ay, r1, r2);

            if h.z < 0.0 {
                h = -h;
            }

            // TODO: Refactor into metallic BRDF and specular BSDF
            let fresnel = self.disney_fresnel(&state.material, state.eta, dot(l, &h), dot(&v, &h));
            let ff = 1.0 - ((1.0 - fresnel) * state.material.spec_trans * (1.0 - state.material.metallic));

            let rand : F = rng.gen();

            if rand < ff {
                *l = normalize(&reflect(-v, h));

                f = self.eval_spec_reflection(&state.material, state.eta, &spec_col, &v, l, &h, pdf);
                *pdf *= ff;
            } else {
                *l = normalize(&refract(-v, h, state.eta));

                f = self.eval_spec_refraction(&state.material, state.eta, &v, l, &h, pdf);
                *pdf *= 1.0 - ff;
            }

            *pdf *= spec_reflect_wt + spec_refract_wt;
        }

        *l = to_world(&t, &b, &n, &l);
        dot(n, l).abs() * f
    }

    fn disney_eval(&self, state: &State, v: F3, n: &F3, l: &F3, bsdf_pdf: &mut F) -> F3 {
        *bsdf_pdf = 0.0;
        let mut f = F3::zeros();

        fn onb(n: &F3, t: &mut F3, b: &mut F3) {
            let up = if n.z.abs() < 0.999 { F3::new(0.0, 0.0, 1.0) } else { F3::new(1.0, 0.0, 0.0) };

            *t = normalize(&cross(&up, &n));
            *b = cross(&n, &t);
        }

        fn to_local(x: &F3, y: &F3, z: &F3, v: &F3) -> F3 {
                F3::new(dot(v, x), dot(v, y), dot(v, z))
        }

        let mut t = F3::zeros();
        let mut b = F3::zeros();

        onb(n, &mut t, &mut b);
        let v = to_local(&t, &b, n, &v); // NDotL = L.z; NDotV = V.z; NDotH = H.z
        let l = to_local(&t, &b, n, l);

        let mut h;
        if l.z > 0.0 {
            h = normalize(&(l + v));
        } else {
            h = normalize(&(l + state.eta * v));
        }

        if h.z < 0.0 {
            h = -h;
        }

        // Specular and sheen color
        let mut spec_col = F3::zeros();
        let mut sheen_col = F3::zeros();
        self.get_spec_color(&state.material, state.eta, &mut spec_col, &mut sheen_col);

        // Lobe weights
        let mut diffuse_wt = 0.0; let mut spec_reflect_wt = 0.0; let mut spec_refract_wt = 0.0; let mut clearcoat_wt = 0.0;

        let fresnel = self.disney_fresnel(&state.material, state.eta, dot(&l, &h), dot(&v, &h));
        self.get_lobe_probabilities(&state.material, &F3::new_x(state.eta), &spec_col, fresnel, &mut diffuse_wt, &mut spec_reflect_wt, &mut spec_refract_wt, &mut clearcoat_wt);

        let mut pdf = 0.0;

        // Diffuse
        if diffuse_wt > 0.0 && l.z > 0.0 {
            f += self.eval_diffuse(&state.material, &sheen_col, &v, &l, &h, &mut pdf);
            *bsdf_pdf += pdf * diffuse_wt;
        }

        // Specular Reflection
        if spec_reflect_wt > 0.0 && l.z > 0.0 && v.z > 0.0 {
            f += self.eval_spec_reflection(&state.material, state.eta, &spec_col, &v, &l, &h, &mut pdf);
            *bsdf_pdf += pdf * spec_reflect_wt;
        }

        // Specular Refraction
        if spec_refract_wt > 0.0 && l.z < 0.0 {
            f += self.eval_spec_refraction(&state.material, state.eta, &v, &l, &h, &mut pdf);
            *bsdf_pdf += pdf * spec_refract_wt;
        }

        // Clearcoat
        if clearcoat_wt > 0.0 && l.z > 0.0 && v.z > 0.0 {
            f +=self.eval_clearcoat(&state.material, &v, &l, &h, &mut pdf);
            *bsdf_pdf += pdf * clearcoat_wt;
        }

        l.z.abs() * f
    }


}