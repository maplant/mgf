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
use mgf::*;

use rand::{SeedableRng, StdRng};
use rand::distributions::{Range, IndependentSample};

use genmesh;
use genmesh::{Vertices, Triangulate};
use genmesh::generators::{SphereUV, Cylinder, SharedVertex, IndexedPolygon};

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
    bodies: RigidBodyVec,
    bvh_ids: Vec<usize>,
    bvh: BVH<AABB, usize>,
    terrain: Mesh,
    locals: Buffer<R, Locals>,
    sphere_model: (Buffer<R, Vertex>, Slice<R>),
    cylinder_model: (Buffer<R, Vertex>, Slice<R>),
    terrain_model: (Buffer<R, Vertex>, Slice<R>),    
    pipe_state: PipelineState<R, pipe::Meta>,
}

pub const SCREEN_WIDTH: u32 = 1920;
pub const SCREEN_HEIGHT: u32 = 1080;
pub const ZFAR: f32 = 1024.0; 
pub const ZNEAR: f32 = 0.1;

impl<R: Resources> World<R> {
    pub fn new<F: Factory<R>>(factory: &mut F) -> Self {
        let shaders = ShaderSet::Simple(
            factory.create_shader_vertex(
                include_bytes!(shader_file!("balls_vs.glsl")))
                .expect("failed to compile vertex shader"),
            factory.create_shader_pixel(
                include_bytes!(shader_file!("balls_fs.glsl")))
                .expect("failed to compile fragment shader")
        );

        let sphere = SphereUV::new(25, 25);
        let sphere_verts: Vec<Vertex> = sphere.shared_vertex_iter()
            .map(|genmesh::Vertex{ pos, .. }| {
                Vertex{ pos }
            })
            .collect();
        let sphere_inds: Vec<u32> = sphere.indexed_polygon_iter()
            .triangulate()
            .vertices()
            .map(|i| i as u32)
            .collect();
        let cylinder = Cylinder::new(25);
        let cylinder_verts: Vec<Vertex> = cylinder.shared_vertex_iter()
            .map(|genmesh::Vertex{ pos, .. }| {
                Vertex{ pos: [pos[0], pos[2] / 2.0,  pos[1]] }
            })
            .collect();
        let cylinder_inds: Vec<u32> = cylinder.indexed_polygon_iter()
            .triangulate()
            .vertices()
            .map(|i| i as u32)
            .collect();
        let mut terrain_mesh = Mesh::new();
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
        for vert in terrain_verts.iter() {
            terrain_mesh.push_vert(Point3::from(vert.pos));
        }
        let terrain_inds = vec![
            0u32, 1, 3,
            1, 2, 3,
        ];
        // It is extremely important to ensure that the triangle formed has the correct,
        // intended normal vector. This does not come from the mesh's stored value per vertex,
        // but from the ordering of the points
        terrain_mesh.push_face((0, 1, 3));
        terrain_mesh.push_face((1, 2, 3));
        terrain_mesh.push_face((0, 5, 1));
        terrain_mesh.push_face((0, 4, 5));
        terrain_mesh.push_face((0, 3, 7));
        terrain_mesh.push_face((0, 7, 4));
        terrain_mesh.push_face((2, 6, 3));
        terrain_mesh.push_face((3, 6, 7));
        terrain_mesh.push_face((1, 5, 2));
        terrain_mesh.push_face((2, 5, 6));
        terrain_mesh.set_pos(Point3::new(0.0, -10.0, 0.0));
        World {
            rot_x: 0.0,
            rot_y: 0.0,
            cam_pos: Point3::new(-20.0, 5.0, 0.0),
            cam_dir: Vector3::unit_x(),
            cam_up: Vector3::unit_y(),
            bodies: RigidBodyVec::new(),
            bvh_ids: Vec::new(),
            bvh: BVH::new(),
            terrain:  terrain_mesh,
            locals: factory.create_constant_buffer(1),
            sphere_model:  factory.create_vertex_buffer_with_slice(
                &sphere_verts, &sphere_inds[..]
            ),
            cylinder_model: factory.create_vertex_buffer_with_slice(
                &cylinder_verts, &cylinder_inds[..]
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

    pub fn add_body(&mut self, collider: Component, mass: f32, restitution: f32, friction: f32, world_force: Vector3<f32>) -> usize {
        let id: usize = self.bodies.add_body(collider, mass, restitution, friction, world_force).into();
        let bounds: AABB = self.bodies.collider[id].bounds();
        let bvh_id = self.bvh.insert(&(bounds + 0.25), id);
        self.bvh_ids.push(bvh_id);
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
        let mut solver = Solver::<ContactConstraint<_>>::new();

        self.bodies.complete_motion();
        self.bodies.integrate(dt);

        for (i, collider) in self.bodies.colliders().enumerate() {
            let bounds: AABB = collider.bounds();
            if !self.bvh[self.bvh_ids[i]].contains(&bounds) {
                self.bvh.remove(self.bvh_ids[i]);
                self.bvh_ids[i] = self.bvh.insert(&(bounds + 0.25), i);
            }

            collider.local_contacts(
                &self.terrain,
                | lc | {
                    solver.add_constraint(
                        ContactConstraint::new(
                            &self.bodies,
                            RigidBodyRef::Dynamic(i),
                            RigidBodyRef::Static{ center: self.terrain.center(), friction: 0.0 },
                            Manifold::from(lc),
                            dt,
                        )
                    )
                }
            );
            
            // If we're the first body integrated we don't need to do anything.
            if i == 0 {
                continue;
            }

            let bvh = &self.bvh;
            bvh.query(
                &bounds,
                |&collider_i| {
                    // For rigid body collisions, collect the contacts into a pruner
                    // and then put that into a manifold.
                    if collider_i >= i {
                        return;
                    }
                    let mut pruner: ContactPruner = ContactPruner::new();
                    collider.local_contacts(
                        &self.bodies.collider[collider_i],
                        |lc| {
                            pruner.push(lc);
                        }
                    );
                    let manifold = Manifold::from(pruner);
                    if manifold.len() == 0 {
                        return;
                    }
                    solver.add_constraint(
                        ContactConstraint::new(
                            &self.bodies,
                            RigidBodyRef::Dynamic(i),
                            RigidBodyRef::Dynamic(collider_i),
                            manifold,
                            dt,
                        )
                    );
                }
            );
        }

        solver.solve(&mut self.bodies, 20);
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

        for (i, &collider) in self.bodies.collider.iter().enumerate() {
            match collider {
                Moving(Component::Sphere(s),_) => {
                    let locals = Locals {
                        color: [ between.ind_sample(&mut rng),
                                 between.ind_sample(&mut rng),
                                 between.ind_sample(&mut rng),
                                 1.0 ],
                        model: (Matrix4::from_translation(self.bodies.x[i].to_vec())
                                * Matrix4::from_scale(s.r)).into(),
                        view: view.into(),
                        proj: proj.into(),
                    };
                    encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
                    encoder.draw(&self.sphere_model.1, &self.pipe_state, &data);
                },

                Moving(Component::Capsule(c),_) => {
                    let color = [ between.ind_sample(&mut rng),
                                  between.ind_sample(&mut rng),
                                  between.ind_sample(&mut rng),
                                  1.0 ];
                    let locals = Locals {
                        color,
                        model: (Matrix4::from_translation(c.a.to_vec())
                                * Matrix4::from_scale(c.r)).into(),
                        view: view.into(),
                        proj: proj.into(),
                    };
                    encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
                    encoder.draw(&self.sphere_model.1, &self.pipe_state, &data);
                    let locals = Locals {
                        color,
                        model: (Matrix4::from_translation(c.a.to_vec() + c.d)
                                * Matrix4::from_scale(c.r)).into(),
                        view: view.into(),
                        proj: proj.into(),
                    };
                    encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
                    encoder.draw(&self.sphere_model.1, &self.pipe_state, &data);
                    let locals = Locals {
                        color,
                        model: (Matrix4::from_translation(c.center().to_vec())
                                * Matrix4::from(self.bodies.q[i])
                                * Matrix4::from_nonuniform_scale(c.r, c.d.magnitude(), c.r)).into(),
                        view: view.into(),
                        proj: proj.into(),
                    };
                    data.vbuf = self.cylinder_model.0.clone();
                    encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
                    encoder.draw(&self.cylinder_model.1, &self.pipe_state, &data);
                    data.vbuf = self.sphere_model.0.clone();
                },
            }
        }
        data.vbuf = self.terrain_model.0.clone();
        let locals = Locals {
            color: [ 0.3, 0.25, 0.55, 1.0 ],
            model: Matrix4::from_translation(self.terrain.center().to_vec()).into(),
            view: view.into(),
            proj: proj.into(),
        };
        encoder.update_buffer(&data.locals, &[locals], 0).unwrap();
        encoder.draw(&self.terrain_model.1, &self.pipe_state, &data);
    }
}
            
        
