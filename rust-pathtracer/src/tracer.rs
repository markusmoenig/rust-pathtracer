use rayon::{slice::ParallelSliceMut, iter::{IndexedParallelIterator, ParallelIterator}};
use rand::{thread_rng, Rng};
use crate::prelude::*;

pub struct Tracer {

    width               : usize,
    height              : usize,

    camera              : Box<dyn Camera3D>,
    scene               : Box<dyn Scene>,

    pixels              : Vec<PTF>,

    frames              : usize,
}

impl Tracer {

    pub fn new(width: usize, height: usize, scene: Box<dyn Scene>) -> Self {
        Self {
            width,
            height,

            camera      : Box::new(Pinhole::new()),
            scene,

            pixels      : vec![0.0; width * height * 4],
            frames      : 0,
        }
    }

    /// Render one frame and accumulate into the pixels buffer
    pub fn render(&mut self) {

        const LINES: usize = 1;

        let width = self.width;
        let height = self.height as PTF;

        self.pixels
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

                    let mut color = [0.0, 0.0, 0.0, 0.0];

                    if self.scene.hit(&ray) {
                        color[0] = 1.0;
                        color[3] = 1.0;
                    }

                    #[inline(always)]
                    pub fn mix_color(a: &[PTF], b: &[PTF], v: PTF) -> [PTF; 4] {
                        [   (1.0 - v) * a[0] + b[0] * v,
                            (1.0 - v) * a[1] + b[1] * v,
                            (1.0 - v) * a[2] + b[2] * v,
                            (1.0 - v) * a[3] + b[3] * v ]
                    }

                    let mix = mix_color(pixel, &color, 1.0 / (self.frames + 1) as PTF);

                    pixel.copy_from_slice(&mix);
                }
            });

        self.frames += 1;

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