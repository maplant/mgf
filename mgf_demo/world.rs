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

use cgmath;
use cgmath::*;
use mgf;
use mgf::*;

use rand::{SeedableRng, StdRng};
use rand::distributions::{Range, IndependentSample};

use genmesh;
use genmesh::{Vertices, Triangulate};
use genmesh::generators::{SphereUV, SharedVertex, IndexedPolygon};

use gfx;
use gfx::{CommandBuffer, Encoder, PipelineState, Primitive, Resources, Slice, ShaderSet};
use gfx::handle::{Buffer, RenderTargetView, DepthStencilView};
use gfx::traits::{Factory, FactoryExt};

use input::*;

type ColorFormat = gfx::format::Srgba8;
type DepthFormat =  gfx::format::DepthStencil;


macro_rules! shader_file {
    ($resource:expr) => (concat!(env!("CARGO_MANIFEST_DIR"), "/shaders/", $resource));
}

gfx_defines! {
    vertex Vertex {
        pos: [f32; 3] = "v_pos",
    }

    constant Locals {
        color: [f32; 4] = "u_color",
        model: [[f32; 4]; 4] = "u_model",
        view: [[f32; 4]; 4] = "u_view",
        proj: [[f32; 4]; 4] = "u_proj",
    }

    pipeline pipe {
        vbuf: gfx::VertexBuffer<Vertex> = (),
        locals: gfx::ConstantBuffer<Locals> = "Locals",
        out_color: gfx::RenderTarget<ColorFormat> = "Target0",
        out_depth: gfx::DepthTarget<DepthFormat> =
            gfx::preset::depth::LESS_EQUAL_WRITE,
    }
}

pub struct World<R: Resources> {
    rot_x: f32,
    rot_y: f32,
    cam_pos: Point3<f32>,
    cam_dir: Vector3<f32>,
    cam_up: Vector3<f32>,
    bodies: Vec<(SimpleDynamicBody<Sphere>, usize)>,
    bvh: BVH<AABB, usize>,
    floor_terrain: Mesh,
    wall_terrain: Mesh,
    locals: Buffer<R, Locals>,
    sphere_model: (Buffer<R, Vertex>, Slice<R>),
    terrain_model: (Buffer<R, Vertex>, Slice<R>),    
    pipe_state: PipelineState<R, pipe::Meta>,
}

pub const SCREEN_WIDTH: u32 = 1920;
pub const SCREEN_HEIGHT: u32 = 1080;
pub const ZFAR: f32 = 1024.0; 
pub const ZNEAR: f32 = 0.1;

impl<R: Resources> World<R> {
    pub fn new<F: Factory<R>>(factory: &mut F) -> Self {
        // Generate sphere
        let sphere = SphereUV::new(25, 25);
        let vertex_data: Vec<Vertex> = sphere.shared_vertex_iter()
            .map(|genmesh::Vertex{ pos, .. }| {
                Vertex{ pos }
            })
            .collect();
        let index_data: Vec<u32> = sphere.indexed_polygon_iter()
            .triangulate()
            .vertices()
            .map(|i| i as u32)
            .collect();
        let shaders = ShaderSet::Simple(
            factory.create_shader_vertex(
                include_bytes!(shader_file!("balls_vs.glsl")))
                .expect("failed to compile vertex shader"),
            factory.create_shader_pixel(
                include_bytes!(shader_file!("balls_fs.glsl")))
                .expect("failed to compile fragment shader")
        );
        // Generate terrain
        let mut floor_terrain_mesh = Mesh::new();
        let mut wall_terrain_mesh = Mesh::new();
        let terrain_verts = [
            Vertex{ pos: [ -10.0, 0.0, -10.0 ] },
            Vertex{ pos: [ -10.0, 0.0, 10.0 ] },
            Vertex{ pos: [ 10.0, 0.0, 10.0 ] },
            Vertex{ pos: [ 10.0, 0.0, -10.0 ] },
            Vertex{ pos: [ -10.0, 10.0, -10.0 ] },
            Vertex{ pos: [ -10.0, 10.0, 10.0 ] },
            Vertex{ pos: [ 10.0, 10.0, 10.0 ] },
            Vertex{ pos: [ 10.0, 10.0, -10.0 ] },

        ];
        for vert_i in 0..4 { 
            floor_terrain_mesh.push_vert(mgf::Vertex{
                p: Point3::from(terrain_verts[vert_i].pos),
                n: Vector3::zero() // Not used
            });
        }
        for vert in terrain_verts.iter() {
            wall_terrain_mesh.push_vert(mgf::Vertex{
                p: Point3::from(vert.pos),
                n: Vector3::zero()
            });
        }
        // Don't show the walls.
        let terrain_inds = vec![
            0u32, 1, 3,
            1, 2, 3,
        ];
        // It is extremely important to ensure that the triangle formed has the correct,
        // intended normal vector. This does not come from the mesh's stored value per vertex,
        floor_terrain_mesh.push_face((0, 1, 3));
        floor_terrain_mesh.push_face((1, 2, 3));
        floor_terrain_mesh.set_pos(Point3::new(0.0, -10.0, 0.0));
        // Extremely discontinuous surfaces, i.e. changes in slope of 90 degrees, is not
        // handled particularly well with one StaticBody.
        wall_terrain_mesh.push_face((0, 5, 1));
        wall_terrain_mesh.push_face((0, 4, 5));
        wall_terrain_mesh.push_face((0, 3, 7));
        wall_terrain_mesh.push_face((0, 7, 4));
        wall_terrain_mesh.push_face((2, 6, 3));
        wall_terrain_mesh.push_face((3, 6, 7));
        wall_terrain_mesh.push_face((1, 5, 2));
        wall_terrain_mesh.push_face((2, 5, 6));
        wall_terrain_mesh.set_pos(Point3::new(0.0, -10.0, 0.0));
        World {
            rot_x: 0.0,
            rot_y: 0.0,
            cam_pos: Point3::new(-20.0, 5.0, 0.0),
            cam_dir: Vector3::unit_x(),
            cam_up: Vector3::unit_y(),
            bodies: Vec::new(),
            bvh: BVH::new(),
            floor_terrain:  floor_terrain_mesh,
            wall_terrain: wall_terrain_mesh,
            locals: factory.create_constant_buffer(1),
            sphere_model:  factory.create_vertex_buffer_with_slice(
                &vertex_data, &index_data[..]
            ),
            terrain_model: factory.create_vertex_buffer_with_slice(
                &terrain_verts, &terrain_inds[..]
            ),
            pipe_state: factory.create_pipeline_state(
                &shaders, Primitive::TriangleList, gfx::state::Rasterizer::new_fill(),
                pipe::new()
            ).unwrap(),
        }
    }

    pub fn insert_sphere(&mut self, r: SimpleDynamicBody<Sphere>) -> usize {
        let id = self.bodies.len();
        let bounds: AABB = r.bounds();
        let bvh_id = self.bvh.insert(&(bounds + 5.0), id);
        self.bodies.push((r, bvh_id));
        id
    }

    pub fn enter_frame(&mut self, input: &Input, dt: f32) {
        self.rot_x -= input.delta_x as f32 * 0.05;
        self.rot_y += input.delta_y as f32 * 0.05;
        self.rot_y = if self.rot_y > 90.0 { 90.0 }
        else if self.rot_y < -90.0 { -90.0 } else { self.rot_y };
        let q_x = Quaternion::from_axis_angle(Vector3::unit_y(), Deg(self.rot_x));
        let cam_dir = q_x.rotate_vector(Vector3::unit_x());
        let q_y = Quaternion::from_axis_angle(Vector3::unit_y().cross(cam_dir),
                                              Deg(self.rot_y));
        self.cam_dir = q_y.rotate_vector(cam_dir);
        self.cam_up = q_y.rotate_vector(Vector3::unit_y());
        self.cam_pos +=
             self.cam_dir * 0.5 *
        // Determine if we are moving forwards or back
            if input.move_forward {
                if input.move_backward {
                    0.0
                } else {
                    1.0
                }
            } else if input.move_backward {
                -1.0
            } else {
                0.0
            } +
            self.cam_up.cross(self.cam_dir) * 0.5 *
        // Determine if we are strafing
            if input.strafe_left {
                if input.strafe_right {
                    0.0
                } else {
                    1.0
                }
            } else if input.strafe_right {
                -1.0
            } else {
                0.0
            };
        self.step(dt);
    }

    fn step(&mut self, dt: f32) {
        let mut floor_body = StaticBody::new(0.5, &self.floor_terrain);
        let mut wall_body = StaticBody::new(0.5, &self.wall_terrain);
        let mut contact_solver: ContactSolver = ContactSolver::new();
        for body_i in 0..self.bodies.len() {
            let (left, right) = (&mut self.bodies).split_at_mut(body_i);
            right[0].0.integrate(dt);
            let bounds: AABB = right[0].0.bounds();
            if !self.bvh[right[0].1].contains(&bounds) {
                self.bvh.remove(right[0].1);
                right[0].1 = self.bvh.insert(&(bounds + 5.0), body_i);
            }
            let bvh = &self.bvh;
            bvh.query(&bounds, |&id| {
                if id >= body_i {
                    return;
                }
                let mut pruner: ContactPruner = ContactPruner::new();
                right[0].0.local_contacts(&left[id].0, |lc|{ pruner.push(lc); });
                unsafe {
                    let body_a = &mut right[0].0 as *mut SimpleDynamicBody<Sphere>;
                    let body_b = &mut left[id].0 as *mut SimpleDynamicBody<Sphere>;
                    contact_solver.add_constraint(
                        &mut *body_a,
                        &mut *body_b,
                        Manifold::from(pruner),
                        dt
                    );
                }
            });
            let mut pruner: ContactPruner = ContactPruner::new();
            right[0].0.local_contacts(&floor_body, |lc|{ pruner.push(lc); });
            unsafe {
                let body = &mut right[0].0 as *mut SimpleDynamicBody<Sphere>;
                let floor_body = &mut floor_body as *mut StaticBody<Mesh>;
                contact_solver.add_constraint(
                    &mut *body,
                    &mut *floor_body,
                    Manifold::from(pruner),
                    dt
                );
            }
            let mut pruner: ContactPruner = ContactPruner::new();
            right[0].0.local_contacts(&wall_body, |lc|{ pruner.push(lc); });
            unsafe {
                let body = &mut right[0].0 as *mut SimpleDynamicBody<Sphere>;
                let wall_body = &mut wall_body as *mut StaticBody<Mesh>;
                contact_solver.add_constraint(
                    &mut *body,
                    &mut *wall_body,
                    Manifold::from(pruner),
                    dt
                );
            }
        }
        contact_solver.solve(20);
    }
    
    pub fn render<C>(
        &mut self,
        encoder: &mut Encoder<R, C>,
        color: RenderTargetView<R, ColorFormat>,
        depth: DepthStencilView<R, DepthFormat>,
    ) where
        C: CommandBuffer<R>
    {
        let seed: &[_] = &[1, 2, 3, 4];
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let between = Range::new(0f32, 1.0);

        encoder.clear(&color, [1.0, 1.0, 1.0, 1.0]);
        encoder.clear_depth(&depth, 1.0);

        let aspect_ratio = (SCREEN_WIDTH as f32) / (SCREEN_HEIGHT as f32);
        let proj = cgmath::perspective(Deg(90.0), aspect_ratio, ZNEAR, ZFAR);

        let view = Matrix4::look_at(
            self.cam_pos,
            self.cam_pos + self.cam_dir,
            self.cam_up 
        );

        let mut data = pipe::Data {
            vbuf: self.sphere_model.0.clone(),
            locals: self.locals.clone(),
            out_color: color,
            out_depth: depth,
        };

        for body in self.bodies.iter() {
            let locals = Locals {
                color: [ between.ind_sample(&mut rng),
                         between.ind_sample(&mut rng),
                         between.ind_sample(&mut rng),
                         1.0 ],
                model: (Matrix4::from_translation(body.0.center().to_vec())
                        * Matrix4::from_scale(body.0.collider.0.r)).into(),
                view: view.into(),
                proj: proj.into(),
            };
            encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
            encoder.draw(&self.sphere_model.1, &self.pipe_state, &data);
        }
        data.vbuf = self.terrain_model.0.clone();
        let locals = Locals {
            color: [ 0.3, 0.25, 0.55, 1.0 ],
            model: Matrix4::from_translation(self.floor_terrain.center().to_vec()).into(),
            view: view.into(),
            proj: proj.into(),
        };
        encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
        encoder.draw(&self.terrain_model.1, &self.pipe_state, &data);
    }
}
            
            
