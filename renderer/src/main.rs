#![deny(clippy::all)]
#![forbid(unsafe_code)]

mod analytical;

use analytical::AnalyticalScene;
use rust_pathtracer::prelude::*;

use log::error;
use pixels::{Error, Pixels, SurfaceTexture};
use tao::{
    dpi::PhysicalPosition,
    dpi::LogicalSize,
    event::{Event, DeviceEvent, ElementState, KeyEvent, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    keyboard::KeyCode,
    menu::{MenuBar, MenuItem},
    window::WindowBuilder,
    //keyboard::{Key},
};

extern crate nalgebra_glm as glm;


fn main() -> Result<(), Error> {

    let width     : usize = 800;
    let height    : usize = 600;

    let mut buffer = ColorBuffer::new(width, height);

    let scene = Box::new(AnalyticalScene::new());
    let mut pt = Tracer::new(scene);

    env_logger::init();
    let event_loop = EventLoop::new();
    let window = {
        let mut file_menu = MenuBar::new();
        file_menu.add_native_item(MenuItem::Quit);

        let mut menu = MenuBar::new();
        menu.add_submenu("File", true, file_menu);

        let size = LogicalSize::new(width as f64, height as f64);
        WindowBuilder::new()
            .with_title("Rust-Pathtracer")
            .with_menu(menu)
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .unwrap()
    };

    let mut pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width, window_size.height, &window);
        Pixels::new(width as u32, height as u32, surface_texture)?
    };

    // Init the code editor

    let mut coords = PhysicalPosition::new(0.0, 0.0);
    let mut is_pressed = false;

    event_loop.run(move |event, _, control_flow| {
        match event {
            Event::WindowEvent { event, .. } => match event {
                // Close events
                WindowEvent::CloseRequested
                | WindowEvent::KeyboardInput {
                    event:
                        KeyEvent {
                            physical_key: KeyCode::Escape,
                            ..
                        },
                    ..
                } => {
                    *control_flow = ControlFlow::Exit;
                }

                // Resize the window
                WindowEvent::Resized(size) => {
                    pixels.resize_surface(size.width, size.height);
                    // let scale = window.scale_factor() as u32;
                    // pixels.resize_buffer(size.width / scale, size.height / scale);
                    // width = size.width as usize / scale as usize;
                    // height = size.height as usize / scale as usize;
                    window.request_redraw();
                }

                WindowEvent::CursorMoved { position, .. } => {
                    coords = position;
                }

                _ => (),
            },

            // Update internal state and request a redraw
            Event::MainEventsCleared => {
                window.request_redraw();
            }

            // Draw the current frame
            Event::RedrawRequested(_) => {

                let frame = pixels.get_frame();
                pt.render(&mut buffer);
                buffer.convert_to_u8(frame);
                //editor_lib::rust_draw(frame.as_mut_ptr(), width as u32, height as u32);

                if pixels
                    .render()
                    .map_err(|e| error!("pixels.render() failed: {}", e))
                    .is_err()
                {
                    *control_flow = ControlFlow::Exit;
                }
            },

            Event::DeviceEvent { event, .. } => match event {
                DeviceEvent::MouseMotion { /*delta,*/ .. } => {
                    //println!("mouse moved: {:?}", delta),
                    if let Some(_pixel_pos) = pixels.window_pos_to_pixel((coords.x as f32, coords.y as f32)).ok() {
                        if is_pressed {
                        //     if editor_lib::rust_touch_dragged(pixel_pos.0 as f32, pixel_pos.1 as f32) {
                        //         window.request_redraw();
                        //     }
                        } else {

                        }
                        // if editor_lib::rust_hover(pixel_pos.0 as f32, pixel_pos.1 as f32) {
                        //     window.request_redraw();
                        // }
                    }
                }
                DeviceEvent::Button {state, .. } => match state {
                    ElementState::Pressed => {
                        //println!("mouse button {} pressed", button);
                        if let Some(_pixel_pos) = pixels.window_pos_to_pixel((coords.x as f32, coords.y as f32)).ok() {
                            is_pressed = true;
                            // if editor_lib::rust_touch_down(pixel_pos.0 as f32, pixel_pos.1 as f32) {
                            //     window.request_redraw();
                            // }
                        }
                    }
                    ElementState::Released => {
                        //println!("mouse button {} released", button),
                        if let Some(_pixel_pos) = pixels.window_pos_to_pixel((coords.x as f32, coords.y as f32)).ok() {
                            is_pressed = false;
                            // if editor_lib::rust_touch_up(pixel_pos.0 as f32, pixel_pos.1 as f32) {
                            //     window.request_redraw();
                            // }
                        }
                    }
                    _ => (),
                },

                DeviceEvent::MouseWheel { delta, .. } => match delta {
                    // tao::event::MouseScrollDelta::LineDelta(x, y) => {
                    //     println!("mouse wheel Line Delta: ({},{})", x, y);
                    //     let pixels_per_line = 120.0;
                    //     let mut pos = window.outer_position().unwrap();
                    //     pos.x -= (x * pixels_per_line) as i32;
                    //     pos.y -= (y * pixels_per_line) as i32;
                    //     window.set_outer_position(pos)
                    // }
                    tao::event::MouseScrollDelta::PixelDelta(_p) => {
                        //println!("mouse wheel Pixel Delta: ({},{})", p.x, p.y);
                        //if editor.mouse_wheel((p.x as isize, p.y as isize)) {
                        //    window.request_redraw();
                            //mouse_wheel_ongoing = true;
                        //}
                    }
                    _ => (),
                },
                _ => (),
            }
            _ => {}
        }
    });
}