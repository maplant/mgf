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

use std::ops::{AddAssign, SubAssign};
use std::vec::Vec;

use bvh::*;
use geom::*;
use collision::*;
use bounds::{BoundedBy};
use cgmath::prelude::*;
use cgmath::{Vector3, Point3, Zero, Rotation3};

/// A triangle mesh is a set of triangles that forms some sort of mesh. There
/// are no requirements on the convexivity of the mesh.
#[derive(Clone)]
pub struct Mesh {
    pub x: Vector3<f32>,
    pub verts: Vec<Point3<f32>>,
    pub faces: Vec<(usize, usize, usize)>,
    pub bvh: BVH<AABB, usize>, 
}

impl Mesh {
    pub fn new() -> Self {
        Mesh {
            x: Vector3::zero(),
            verts: Vec::new(),
            faces: Vec::new(),
            bvh: BVH::new(),
        }
    }

    pub fn with_capacity(cap_verts: usize, cap_faces: usize) -> Self {
        Mesh {
            x: Vector3::zero(),
            verts: Vec::with_capacity(cap_verts),
            faces: Vec::with_capacity(cap_faces),
            bvh: BVH::with_capacity(cap_faces),
        }
    }

    pub fn push_vert(&mut self, p: Point3<f32>) -> usize {
        let id = self.verts.len();
        self.verts.push(p);
        id
    }

    pub fn push_face(&mut self, f: (usize, usize, usize)) -> usize {
        let a = self.verts[f.0];
        let b = self.verts[f.1];
        let c = self.verts[f.2];
        let tri = Triangle::from((a, b, c));
        let index = self.faces.len();
        self.faces.push(f);
        self.bvh.insert(&tri, index);
        index
    }
}

impl AddAssign<Vector3<f32>> for Mesh {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.x += v;
    }
}

impl SubAssign<Vector3<f32>> for Mesh {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.x -= v;
    }
}

impl Shape for Mesh {
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.x)
    }
}

/// Rotating meshes is not a fast operation.
impl Volumetric for Mesh {
    fn rotate<R: Rotation3<f32>>(mut self, rot: R) -> Mesh {
        for vert in self.verts.iter_mut() {
            *vert = rot.rotate_point(*vert);
        }
        self.bvh.clear();
        for (i, &(a, b, c)) in self.faces.iter().enumerate() {
            let tri = Triangle::from(
                (self.verts[a], self.verts[b], self.verts[c])
            );
            self.bvh.insert(&tri, i);
        }
        self
    }
}

impl<RHS> Contacts<RHS> for Mesh
where
    RHS: Contacts<Triangle> + Contacts<Rectangle> + BoundedBy<AABB>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        let mut collided = false;
        self.bvh.query(&(rhs.bounds() - self.x), |&face_index| {
            let (a, b, c) = self.faces[face_index];
            let a = self.verts[a] + self.x;
            let b = self.verts[b] + self.x;
            let c = self.verts[c] + self.x;
            let tri = Triangle::from((a, b, c));
            rhs.contacts(&tri, |c| {
                collided = true;
                callback(Contact {
                    a: c.b,
                    b: c.a,
                    t: c.t,
                    n: -c.n,
                });
            });
        });
        collided
    }
}
