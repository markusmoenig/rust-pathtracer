
use crate::prelude::*;

/// A color buffer holding an array of either f32 or f64 pixels (as defined by PTF).
pub struct ColorBuffer {

    pub width               : usize,
    pub height              : usize,

    pub pixels              : Vec<PTF>,

    pub frames              : usize,
}

impl ColorBuffer {

    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,

            pixels      : vec![0.0; width * height * 4],
            frames      : 0,
        }
    }

    /// Convert the pixel buffer to an Vec<u8> and converts gamma the colors from linear into gamma space.
    pub fn convert_to_u8(&self, frame: &mut [u8]) {
        for y in 0..self.height {
            for x in 0..self.width {
                let o = x * 4 + y * self.width * 4;
                let c = [(self.pixels[o].powf(0.4545) * 255.0) as u8, (self.pixels[o+1].powf(0.4545) * 255.0) as u8, (self.pixels[o+2].powf(0.4545) * 255.0) as u8, (self.pixels[o+3] * 255.0) as u8];
                // let c = [(self.pixels[o]* 255.0) as u8, (self.pixels[o+1] * 255.0) as u8, (self.pixels[o+2]* 255.0) as u8, (self.pixels[o+3] * 255.0) as u8];
                frame[o..o + 4].copy_from_slice(&c);
            }
        }
    }
}