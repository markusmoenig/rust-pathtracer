use crate::prelude::*;

use rhai::Engine;

// Medium

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum MediumType {
    None,
    Absorb,
    Scatter,
    Emissive
}

/// The Medium struct for volumetric objects.
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Medium {
    pub medium_type             : MediumType,
    pub density                 : F,
    pub color                   : F3,
    pub anisotropy              : F,
}

impl Medium {

    pub fn new() -> Self {
        Self {
            medium_type         : MediumType::None,
            density             : 0.0,
            color               : F3::new(0.0, 0.0, 0.0),
            anisotropy          : 0.0,
        }
    }
}

// Material

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum AlphaMode
{
    Opaque,
    Blend,
    Mask
}

/// The material struct holds all BSDF properties as well as the Medium.
#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Material {
    pub rgb                     : F3,
    pub anisotropic             : F,
    pub emission                : F3,

    pub metallic                : F,
    pub roughness               : F,
    pub subsurface              : F,
    pub specular_tint           : F,

    pub sheen                   : F,
    pub sheen_tint              : F,
    pub clearcoat               : F,
    pub clearcoat_gloss         : F,
    /// Do not use clearcoat_roughness directly, it is for internal use only. Use clearcoat_gloss.
    pub clearcoat_roughness     : F,

    pub spec_trans              : F,
    pub ior                     : F,

    pub opacity                 : F,
    pub alpha_mode              : AlphaMode,
    pub alpha_cutoff            : F,

    pub ax                      : F,
    pub ay                      : F,

    pub medium                  : Medium,
}

impl Material {

    pub fn new() -> Self {

        Self {
            rgb                 : F3::new(0.5, 0.5, 0.5),
            emission            : F3::new(0.0, 0.0, 0.0),

            anisotropic         : 0.0,
            metallic            : 0.0,
            roughness           : 0.5,
            subsurface          : 0.0,
            specular_tint       : 0.0,

            sheen               : 0.0,
            sheen_tint          : 0.0,

            clearcoat           : 0.0,
            clearcoat_gloss     : 0.0,
            clearcoat_roughness : 0.0,
            spec_trans          : 0.0,
            ior                 : 1.5,

            opacity             : 1.0,
            alpha_mode          : AlphaMode::Opaque,
            alpha_cutoff        : 0.0,

            medium              : Medium::new(),

            ax                  : 0.0,
            ay                  : 0.0
        }
    }

    /// Material post-processing, called by the tracer after calling Scene::closest_hit()
    pub fn finalize(&mut self) {

        self.roughness = self.roughness.max(0.01);

        fn mix_ptf(a: &F, b: &F, v: F) -> F {
            (1.0 - v) * a + b * v
        }

        self.clearcoat_roughness = mix_ptf(&0.1, &0.001, self.clearcoat_gloss); // Remapping from gloss to roughness
        self.medium.anisotropy = self.medium.anisotropy.clamp(-0.9, 0.9);

        let aspect = (1.0 - self.anisotropic * 0.9).sqrt();
        self.ax = (self.roughness / aspect).max(0.001);
        self.ay = (self.roughness * aspect).max(0.001);
    }

    /// Mixes this material with another material
    pub fn mix(&self, other: &Material, v: F) -> Material {
        let mut m = Material::new();

        m.rgb = mix(&self.rgb, &other.rgb, &v);
        m.emission = mix(&self.emission, &other.emission, &v);

        m.anisotropic = mix_f(&self.anisotropic, &other.anisotropic, &v);
        m.metallic = mix_f(&self.metallic, &other.metallic, &v);
        m.roughness = mix_f(&self.roughness, &other.roughness, &v);
        m.subsurface = mix_f(&self.subsurface, &other.subsurface, &v);
        m.specular_tint = mix_f(&self.specular_tint, &other.specular_tint, &v);

        m.sheen = mix_f(&self.sheen, &other.sheen, &v);
        m.sheen_tint = mix_f(&self.sheen_tint, &other.sheen_tint, &v);

        m.clearcoat = mix_f(&self.clearcoat, &other.clearcoat, &v);
        m.clearcoat_gloss = mix_f(&self.clearcoat_gloss, &other.clearcoat_gloss, &v);
        m.spec_trans = mix_f(&self.spec_trans, &other.spec_trans, &v);
        m.ior = mix_f(&self.ior, &other.ior, &v);

        m
    }

    // --------- Getter / Setter

    pub fn get_rgb(&mut self) -> F3 {
        self.rgb
    }

    pub fn set_rgb(&mut self, new_val: F3) {
        self.rgb = new_val;
    }

    pub fn get_roughness(&mut self) -> F {
        self.roughness
    }

    pub fn set_roughness(&mut self, new_val: F) {
        self.roughness = new_val;
    }

    pub fn get_metallic(&mut self) -> F {
        self.metallic
    }

    pub fn set_metallic(&mut self, new_val: F) {
        self.metallic = new_val;
    }

    /// Register to the engine
    pub fn register(engine: &mut Engine) {
        engine.register_type_with_name::<Material>("Material")
            .register_fn("Material", Material::new)

            .register_get_set("roughness", Material::get_roughness, Material::set_roughness)
            .register_get_set("metallic", Material::get_metallic, Material::set_metallic)

            .register_get_set("rgb", Material::get_rgb, Material::set_rgb);
    }
}