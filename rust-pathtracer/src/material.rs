use crate::prelude::*;

// Medium

#[derive(PartialEq, Clone, Debug)]
pub enum MediumType {
    None,
    Absorb,
    Scatter,
    Emissive
}

#[derive(PartialEq, Clone, Debug)]
pub struct Medium {
    pub medium_type             : MediumType,
    pub density                 : PTF,
    pub color                   : PTF3,
    pub anisotropy              : PTF,
}

impl Medium {

    pub fn new() -> Self {
        Self {
            medium_type         : MediumType::None,
            density             : 0.0,
            color               : PTF3::new(0.0, 0.0, 0.0),
            anisotropy          : 0.0,
        }
    }
}

// Material

#[derive(PartialEq, Clone, Debug)]
pub enum AlphaMode
{
    Opaque,
    Blend,
    Mask
}

#[derive(PartialEq, Clone, Debug)]
pub struct Material {
    pub base_color              : PTF3,
    pub anisotropic             : PTF,
    pub emission                : PTF3,

    pub metallic                : PTF,
    pub roughness               : PTF,
    pub subsurface              : PTF,
    pub specular_tint           : PTF,

    pub sheen                   : PTF,
    pub sheen_tint              : PTF,
    pub clearcoat               : PTF,
    pub clearcoat_gloss         : PTF,

    pub spec_trans              : PTF,
    pub ior                     : PTF,

    pub opacity                 : PTF,
    pub alpha_mode              : AlphaMode,
    pub alpha_cutoff            : PTF,

    pub ax                      : PTF,
    pub ay                      : PTF,

    pub medium                  : Medium,
}

impl Material {

    pub fn new(base_color: PTF3) -> Self {

        let anisotropic = 0.0;
        let roughness = 0.5;

        let aspect = (1.0 - anisotropic * 0.9).sqrt();
        let ax = 0.001.max(roughness / aspect);
        let ay = 0.001.max(roughness * aspect);

        Self {
            base_color,
            emission            : PTF3::new(0.0, 0.0, 0.0),

            anisotropic         : 0.0,
            metallic            : 0.0,
            roughness           : 0.5,
            subsurface          : 0.0,
            specular_tint       : 0.0,

            sheen               : 0.0,
            sheen_tint          : 0.0,

            clearcoat           : 0.0,
            clearcoat_gloss     : 0.0,
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
        let aspect = (1.0 - self.anisotropic * 0.9).sqrt();
        self.ax = 0.001.max(self.roughness / aspect);
        self.ay = 0.001.max(self.roughness * aspect);
    }

}