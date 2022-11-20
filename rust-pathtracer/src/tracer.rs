use crate::prelude::*;

pub struct Tracer {

    width               : usize,
    height              : usize,

    pixels              : Vec<PTF>,
}

impl Tracer {

    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,

            pixels      : vec![0.0; width * height * 4],
        }
    }

    /// Render one frame and accumulate into the pixels buffer
    pub fn render(&mut self) {

    }

    /// Convert the pixel buffer to u8
    pub fn convert_to_u8(&self, frame: &mut [u8]) {
        for y in 0..self.height {
            for x in 0..self.width {
                let d = x * 4 + y  * self.width * 4;
                let s = x * 4 + y * self.width * 4;
                let c = [(self.pixels[s] * 255.0) as u8, (self.pixels[s+1] * 255.0) as u8,  (self.pixels[s+2] * 255.0) as u8,  (self.pixels[s+3] * 255.0) as u8];
                frame[d..d + 4].copy_from_slice(&c);
            }
        }
    }
}