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
use cgmath::{EuclideanSpace, Vector3, Point3, Zero};

/// A point and a normal 
pub struct Vertex {
    pub p: Point3<f32>,
    pub n: Vector3<f32>,
}

/// A triangle mesh is a set of triangles that forms some sort of mesh. There
/// are no requirements on the convexivity of the mesh.
pub struct Mesh {
    pub disp: Vector3<f32>,
    pub verts: Vec<Vertex>,
    pub faces: Vec<(usize, usize, usize)>,
    pub bvh: BVH<AABB, usize>, 
}

impl Mesh {
    pub fn new() -> Self {
        Mesh {
            disp: Vector3::zero(),
            verts: Vec::new(),
            faces: Vec::new(),
            bvh: BVH::new(),
        }
    }

    pub fn with_capacity(cap_verts: usize, cap_faces: usize) -> Self {
        Mesh {
            disp: Vector3::zero(),
            verts: Vec::with_capacity(cap_verts),
            faces: Vec::with_capacity(cap_faces),
            bvh: BVH::with_capacity(cap_faces),
        }
    }

    pub fn push_vert(&mut self, v: Vertex) -> usize {
        let id = self.verts.len();
        self.verts.push(v);
        id
    }

    pub fn push_face(&mut self, f: (usize, usize, usize)) -> usize {
        let a = self.verts[f.0].p;
        let b = self.verts[f.1].p;
        let c = self.verts[f.2].p;
        let tri = Triangle::from((a, b, c));
        let index = self.faces.len();
        self.faces.push(f);
        self.bvh.insert(&tri, index);
        index
    }
}

impl AddAssign<Vector3<f32>> for Mesh {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.disp += v;
    }
}

impl SubAssign<Vector3<f32>> for Mesh {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.disp -= v;
    }
}

impl Shape for Mesh {
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.disp)
    }
}

impl<RHS> Contacts<RHS> for Mesh
where
    RHS: Contacts<Triangle> + Contacts<Rectangle> + BoundedBy<AABB>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        let mut collided = false;
        self.bvh.query(&(rhs.bounds() - self.disp), |&face_index| {
            let (a, b, c) = self.faces[face_index];
            let a = self.verts[a].p + self.disp;
            let b = self.verts[b].p + self.disp;
            let c = self.verts[c].p + self.disp;
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
