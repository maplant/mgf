// Copyright 2017 Matthew Plant. This file is part of MGF.
//
// MGF is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MGF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with MGF. If not, see <http://www.gnu.org/licenses/>.

extern crate mgf;

extern crate rand;

#[macro_use]
extern crate gfx;
extern crate gfx_window_glutin;
extern crate glutin;
extern crate genmesh;

use std::f32;

use mgf::cgmath;

mod input;
#[macro_use]
mod world;

use std::time;

//use cgmath::prelude::*;
use cgmath::*;
use mgf::*;

use input::*;
use world::*;

use gfx::Device;
use gfx_window_glutin as gfx_glutin;
use glutin::{GlContext, GlRequest, VirtualKeyCode};
use glutin::Api::OpenGl;

fn main() {
    let mut events_loop = glutin::EventsLoop::new();
    let windowbuilder = glutin::WindowBuilder::new()
        .with_title("MGF Demo: \"capsules\" (tm)".to_string())
        .with_dimensions(SCREEN_WIDTH, SCREEN_HEIGHT);
    let contextbuilder = glutin::ContextBuilder::new()
        .with_gl(GlRequest::Specific(OpenGl,(3,2)))
        .with_vsync(true);
    let (window, mut device, mut factory, color_view, depth_view) =
        gfx_glutin::init::<gfx::format::Srgba8, gfx::format::DepthStencil>(
            windowbuilder, contextbuilder, &events_loop
        );

    let mut world = World::new(&mut factory);


    // The following code and parameters are taken from nphysics3D in order to
    // make a reasonable comparison.

     let comp = Component::from(
        Component::from(
            Capsule{
                a: Point3::new(-0.5, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0
            }
        ),
     );

    let num     = 1500.0f32.powf(1.0f32 / 3.0) as usize;
    let rad     = 2.0;
    let shift   = 2.5 * rad;
    let centerx = shift * (num as f32) / 2.0;
    let centery = shift * (num as f32) / 2.0;

    for i in 0usize .. num {
        for j in 0usize .. num {
            for k in 0usize .. num {
                let x = i as f32 * 2.5 * rad - centerx;
                let y = 10.0 + j as f32 * 2.5 * rad + centery * 2.0;
                let z = k as f32 * 2.5 * rad - centerx;

                let mut rb = comp.clone();
                rb.set_pos(Point3::new(x, y, z));
                world.add_body(rb, 1.0, 0.3, 0.6, Vector3::new(0.0, -9.8, 0.0));
            }
        }
    }

    let mut input = Input::new();
    input.bind_key(VirtualKeyCode::W, INPUT_UP);
    input.bind_key(VirtualKeyCode::S, INPUT_DOWN);
    input.bind_key(VirtualKeyCode::A, INPUT_LEFT);
    input.bind_key(VirtualKeyCode::D, INPUT_RIGHT);

    let mut encoder: gfx::Encoder<_, _> = factory.create_command_buffer().into();

    while input.gather(&mut events_loop) {
        let start = time::Instant::now();
        world.enter_frame(&input, 1.0 / 60.0);
        let elapsed = start.elapsed();
        print!("Physics step elapsed, took {} ms                \r",
               elapsed.as_secs() * 1000 +
               elapsed.subsec_nanos() as u64 / 1_000_000);

        world.render(&mut encoder, color_view.clone(), depth_view.clone());
        window.swap_buffers().unwrap();
        device.cleanup();
        encoder.flush(&mut device);
    }
    println!("\nDemo finished");
}
