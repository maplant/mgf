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

use crate::bvh::*;
use crate::geom::*;
use crate::collision::*;
use crate::bounds::{BoundedBy};
use cgmath::prelude::*;
use cgmath::{Vector3, Point3, Zero, Rotation3};

use serde::{Serialize, Deserialize};

/// A triangle mesh is a set of triangles that forms some sort of mesh. There
/// are no requirements on the convexivity of the mesh.
#[derive(Clone)]
#[derive(Serialize, Deserialize)]
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

    fn closest_point(&self, _to: Point3<f32>) -> Point3<f32> {
        unimplemented!();
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

/// A closed convex mesh. Represented by a point soup.
#[derive(Clone)]
#[derive(Serialize, Deserialize)]
pub struct ConvexMesh {
    pub x: Vector3<f32>,
    pub sum: Vector3<f32>,
    pub verts: Vec<Point3<f32>>,
}

impl ConvexMesh {
    pub fn new() -> Self {
        ConvexMesh {
            x: Vector3::zero(),
            sum: Vector3::zero(),
            verts: vec![],
        }
    }

    pub fn with_capacity(cap: usize) -> Self {
        ConvexMesh {
            x: Vector3::zero(),
            sum: Vector3::zero(),
            verts: Vec::with_capacity(cap),
        }
    }

    pub fn push(&mut self, p: Point3<f32>) {
        let prev_center = self.sum / self.verts.len() as f32;
        self.sum += p.to_vec();
        self.verts.push(p);
        let new_center = self.sum / self.verts.len() as f32;
        let disp = new_center - prev_center;
        self.x += disp;
    }
}

impl From<Vec<Point3<f32>>> for ConvexMesh {
    fn from(verts: Vec<Point3<f32>>) -> Self {
        let mut sum = Vector3::zero();
        for p in verts.iter() {
            sum += p.to_vec();
        }
        ConvexMesh {
            x: Vector3::zero(),
            sum,
            verts,
        }
    }
}

impl AddAssign<Vector3<f32>> for ConvexMesh {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.x += v;
    }
}

impl SubAssign<Vector3<f32>> for ConvexMesh {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.x -= v;
    }
}

impl Shape for ConvexMesh {
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.x + self.sum / self.verts.len() as f32)
    }

    fn closest_point(&self, _to: Point3<f32>) -> Point3<f32> {
        unimplemented!();
    }
}

impl Volumetric for ConvexMesh {
    fn rotate<R: Rotation3<f32>>(mut self, rot: R) -> ConvexMesh {
        let center = self.sum / self.verts.len() as f32;
        for vert in self.verts.iter_mut() {
            *vert = rot.rotate_point(*vert + -center) + center;
        }
        self
    }
}

impl Convex for ConvexMesh {
    fn support(&self, d: Vector3<f32>) -> Point3<f32> {
        let mut best_vert = self.verts[0];
        let mut best_norm = d.dot(self.verts[0].to_vec());
        for vert in self.verts[1..].iter() {
            let norm = d.dot(vert.to_vec());
            if norm > best_norm {
                best_vert = *vert;
                best_norm = norm;
            }
        }
        best_vert
    }
}
