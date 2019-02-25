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

use std::f32;
use std::fmt;
use std::mem;
use std::collections::HashMap;

use cgmath::prelude::*;
use cgmath::{EuclideanSpace, InnerSpace, Point3, Vector3};
use crate::geom::*;
use crate::collision::*;
use crate::pool::*;

pub struct Simplex<Point = Point3<f32>>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    points: [Point; 4],
    state: &'static SimplexState<Point>,
}

impl<Point> fmt::Debug for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static + fmt::Debug
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        #[derive(Debug)]
        enum StateName {
            Vertex,
            Edge,
            Face,
            Volume
        };
        let (state, points) = if self.state as *const _ == &VERTEX_DATAPTRLOC as *const _ {
            (StateName::Vertex, &self.points[..1])
        } else if self.state as *const _ == &EDGE_DATAPTRLOC as * const _  {
            (StateName::Edge, &self.points[..2])
        } else if self.state as * const _  == &FACE_DATAPTRLOC as * const _ {
            (StateName::Face, &self.points[..3])
        } else {
            (StateName::Volume, &self.points[..4])
        };
        write!(f, "Simplex {{ state: {:?}, points: {:?} }}", state, points)
    }
}

impl<Point> From<[Point; 1]> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static 
{
    fn from(p: [Point; 1]) -> Self {
        Self::from(p[0])
    }
}
          
impl<Point> From<Point> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: Point) -> Self {
        Simplex {
            points: [
                p,
                Point::from(Point3::new(0.0, 0.0, 0.0)),
                Point::from(Point3::new(0.0, 0.0, 0.0)),
                Point::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &VERTEX_DATAPTRLOC,
        }
    }
}

impl<Point> From<[Point; 2]> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Point; 2]) -> Self {
        Self::from((p[0], p[1]))
    }
}

impl<Point> From<(Point, Point)> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Point, Point)) -> Self {
        Simplex {
            points: [
                p.0,
                p.1,
                Point::from(Point3::new(0.0, 0.0, 0.0)),
                Point::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &EDGE_DATAPTRLOC,
        }
    }
}

impl<Point> From<[Point; 3]> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Point; 3]) -> Self {
        Self::from((p[0], p[1], p[2]))
    }
}
          
impl<Point> From<(Point, Point, Point)> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Point, Point, Point)) -> Self {
        Simplex {
            points: [
                p.0,
                p.1,
                p.2,
                Point::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &FACE_DATAPTRLOC,
        }
    }
}

impl<Point> From<[Point; 4]> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Point; 4]) -> Self {
        Self::from((p[0], p[1], p[2], p[3]))
    }
}

impl<Point> From<(Point, Point, Point, Point)> for Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Point, Point, Point, Point)) -> Self {
        Simplex {
            points: [
                p.0,
                p.1,
                p.2,
                p.3,
            ],
            state: &VOLUME_DATAPTRLOC,
        }
    }
}

impl<Point> Simplex<Point>
where
    Point: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static + fmt::Debug
{
    /// Finds the closest point on the shape to the origin
    pub fn closest_point_to_origin<S>(&mut self, shape: &S) -> Point3<f32>
    where
        S: Convex<Point>
    {
        let mut prev_norm = Vector3::zero();
        loop {
            let (min_norm, next_state) = self.state.min_norm(&mut self.points);
            if min_norm.magnitude2() < COLLISION_EPSILON  {
                // If the simplex is not a tetrahedron, we want to sample more 
                // axis until it is one.
                for i in self.state.len()..4 {
                    let min_norm = -Vector3::new(prev_norm.z, prev_norm.x, prev_norm.y);
                    let support = shape.support(-min_norm.normalize());
                    prev_norm = -min_norm.normalize();
                    self.points[i] = support;
                }
                self.state = &VOLUME_DATAPTRLOC;
                return Point3::new(0.0, 0.0, 0.0);
            }
            let support = shape.support(-min_norm.normalize());
            let support_v = support.into().to_vec();
            prev_norm = min_norm;
            if min_norm.magnitude2() >= support_v.magnitude2() {
                return Point3::from_vec(min_norm);
            }
            self.state = next_state;
            self.state.add_point(&mut self.points, support);
        }
    }
}

trait SimplexState<Point>
where
    Point: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, _: &mut [Point; 4]) -> (Vector3<f32>, &'static SimplexState<Point>);

    fn add_point(&self, _: &mut [Point; 4], _: Point);

    fn len(&self) -> usize;
}

struct VertexSimplex{}
struct EdgeSimplex{}
struct FaceSimplex{}
struct VolumeSimplex{}

static VERTEX_DATAPTRLOC: VertexSimplex = VertexSimplex{};
static EDGE_DATAPTRLOC: EdgeSimplex = EdgeSimplex{};
static FACE_DATAPTRLOC: FaceSimplex = FaceSimplex{};
static VOLUME_DATAPTRLOC: VolumeSimplex = VolumeSimplex{};

impl<Point> SimplexState<Point> for VertexSimplex
where
    Point: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Point; 4]) -> (Vector3<f32>, &'static SimplexState<Point>) {
        (simp[0].into().to_vec(), &EDGE_DATAPTRLOC)
    }

    fn add_point(&self, simp:  &mut [Point; 4], p: Point) {
        simp[0] = p
    }

    fn len(&self) -> usize { 1 }
}

impl<Point> SimplexState<Point> for EdgeSimplex
where
    Point: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Point; 4]) -> (Vector3<f32>, &'static SimplexState<Point>) {
        let ab = simp[1].into() - simp[0].into();
        let t = ab.dot(-simp[0].into().to_vec());
        if t <= 0.0 {
            (simp[0].into().to_vec(), &EDGE_DATAPTRLOC)
        } else {
            let denom = ab.dot(ab);
            if t >= denom {
                simp[0] = simp[1];
                (simp[1].into().to_vec(), &EDGE_DATAPTRLOC)
            } else {
                (simp[0].into().to_vec() + ab * (t / denom), &FACE_DATAPTRLOC)
            }
        }
    }

    fn add_point(&self, simp: &mut [Point; 4], p: Point) {
        simp[1] = p
    }

    fn len(&self) -> usize { 2 }
}    


impl<Point> SimplexState<Point> for FaceSimplex
where
    Point: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Point; 4]) -> (Vector3<f32>, &'static SimplexState<Point>) {
        let (a, b, c): (Point3<f32>, Point3<f32>, Point3<f32>) = (
            simp[0].into(), simp[1].into(), simp[2].into()
        );
        let ab = b - a;
        let ac = c - a;
        let ap = -a.to_vec();
        let d1 = ab.dot(ap);
        let d2 = ac.dot(ap);

        // Vertex region A
        if d1 <= 0.0 && d2 <= 0.0 {
            return (simp[0].into().to_vec(), &EDGE_DATAPTRLOC);
        }

        // Vertex region B
        let bp = -b.to_vec();
        let d3 = ab.dot(bp);
        let d4 = ac.dot(bp);
        if d3 >= 0.0 && d4 <= d3 {
            simp[0] = simp[1];
            return (simp[1].into().to_vec(), &EDGE_DATAPTRLOC);
        }

        // Edge region AB
        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return ((simp[0].into().to_vec() + ab * v), &FACE_DATAPTRLOC);
        }

        // Vertex region C
        let cp = -c.to_vec();
        let d5 = ab.dot(cp);
        let d6 = ac.dot(cp);
        if d6 >= 0.0 && d5 <= d6 {
            simp[0] = simp[2];
            return (simp[2].into().to_vec(), &EDGE_DATAPTRLOC);
        }

        // Edge region AC
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            simp[1] = simp[2];
            return ((simp[0].into().to_vec() + ac * w), &FACE_DATAPTRLOC);
        }

        // Edge region BC
        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            simp[0] = simp[2];
            return ((simp[1].into().to_vec() + (simp[2].into() - simp[1].into()) * w), &FACE_DATAPTRLOC);
        }

        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        ((simp[0].into().to_vec() + ab * v + ac * w), &VOLUME_DATAPTRLOC)
    }

    fn add_point(&self, simp:  &mut [Point; 4], p: Point) {
        simp[2] = p;
    }

    fn len(&self) -> usize { 3 }
}

fn origin_outside_plane(
    a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>, d: Vector3<f32>
) -> bool {
    let ab_x_ac = (b - a).cross(c - a);
    let sign_p = (-a).dot(ab_x_ac);
    let sign_d = (d - a).dot(ab_x_ac);
    sign_p * sign_d < 0.0
}

impl<Point> SimplexState<Point> for VolumeSimplex
where
    Point: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Point; 4]) -> (Vector3<f32>, &'static SimplexState<Point>) {
        let mut closest_pt: Vector3<f32> = Vector3::zero();
        let mut best_dist: f32 = f32::INFINITY;
        let mut next_state: &'static SimplexState<Point> = &VERTEX_DATAPTRLOC;
        let (a, b, c, d) = (simp[0], simp[1], simp[2], simp[3]);
        let (av, bv, cv, dv) = (a.into().to_vec(), b.into().to_vec(), c.into().to_vec(), d.into().to_vec());
        // Test face abc
        if origin_outside_plane(av, bv, cv, dv) {
            let mut new_simp = [ a, b, c, d ];
            let (p, new_state) = FACE_DATAPTRLOC.min_norm(&mut new_simp);
            let new_dist = p.magnitude2();
            if new_dist < best_dist {
                closest_pt = p;
                best_dist = new_dist;
                next_state = new_state;
                *simp = new_simp;
            }
        }
        // Test face acd
        if origin_outside_plane(av, cv, dv, bv) {
            let mut new_simp = [ a, c, d, b ];
            let (p, new_state) = FACE_DATAPTRLOC.min_norm(&mut new_simp);
            let new_dist = p.magnitude2();
            if new_dist < best_dist {
                closest_pt = p;
                best_dist = new_dist;
                next_state = new_state;
                *simp = new_simp;
            }
        }
        // Test face adb
        if origin_outside_plane(av, dv, bv, cv) {
            let mut new_simp = [ a, d, b, c ];
            let (p, new_state) = FACE_DATAPTRLOC.min_norm(&mut new_simp);
            let new_dist = p.magnitude2();
            if new_dist < best_dist {
                closest_pt = p;
                best_dist = new_dist;
                next_state = new_state;
                *simp = new_simp;
            }
        }
        // Test face bdc
        if origin_outside_plane(bv, dv, cv, av) {
            let mut new_simp = [ b, d, c, a ];
            let (p, new_state) = FACE_DATAPTRLOC.min_norm(&mut new_simp);
            let new_dist = p.magnitude2();
            if new_dist < best_dist {
                closest_pt = p;
                next_state = new_state;
                *simp = new_simp;
            }
        }

        (closest_pt, next_state)
    }

    fn add_point(&self, simp: &mut [Point; 4], p: Point) {
        simp[3] = p
    }

    fn len(&self) -> usize { 4 }
}

#[derive(Default)]
struct EdgeMap {
    map: HashMap<[ (u32, u32, u32); 2 ], [ (Point3<f32>, Point3<f32>); 2 ]>,
}

impl EdgeMap {
    fn add_edge(&mut self, a: SupportPoint, b: SupportPoint) {
        let la = (a.a, a.b);
        let lb = (b.a, b.b);
        let a = unsafe {
            (
                mem::transmute::<f32, u32>(a.p.x),
                mem::transmute::<f32, u32>(a.p.y),
                mem::transmute::<f32, u32>(a.p.z),
            )
        };
        let b = unsafe {
            (
                mem::transmute::<f32, u32>(b.p.x),
                mem::transmute::<f32, u32>(b.p.y),
                mem::transmute::<f32, u32>(b.p.z),
            )
        };
        let ba = [ b, a ];
        if self.map.contains_key(&ba) {
            self.map.remove(&ba);
            return;
        }
        self.map.insert(
            [ a, b ],
            [ la, lb ]
        );
    }
}
        

impl Simplex<SupportPoint> {
    /// Generates a contact from a simplex. This uses the EPA algorithm. Based on
    /// the description here: http://hacktank.net/blog/?p=119
    pub fn compute_contact<S1, S2>(&self, s1: &S1, s2: &S2) -> Contact
    where
        S1: Convex,
        S2: Convex
    {
        if self.state.len() != 4 {
            panic!("simplex is too small");
        }
        let diff = MinkowskiDiff{ s1, s2 };
        let [ a, b, c, d ] = self.points;
        let mut tris = Pool::from(
            vec![
                (a, b, c),
                (a, c, d),
                (a, d, b),
                (b, d, c)
            ]
        );
        let mut edges = EdgeMap::default();
        const MAX_ITERATIONS: usize = 100;
        for iter in 0..=MAX_ITERATIONS {
            let (closest_dist, closest_i, closest_n) = {
                let mut closest_dist = f32::INFINITY;
                let mut closest_i = 0usize;
                let mut closest_n = Vector3::zero();
                for (i, (a, b, c)) in tris.iter() {
                    let tri = Triangle::from((a.p, b.p, c.p));
                    let n = tri.normal();
                    let dist = n.dot(a.p.to_vec()).abs();
                    if closest_dist > dist {
                        closest_dist = dist;
                        closest_i = i;
                        closest_n = n;
                    }
                }
                (closest_dist, closest_i, closest_n)
            };
            let closest_tri = {
                let (a, b, c) = tris[closest_i];
                ( Triangle::from((a.p, b.p, c.p)), Triangle::from((a.a, b.a, c.a)) )
            };
            let support: SupportPoint = diff.support(closest_n);
            let v = closest_n.dot(support.p.to_vec()) - closest_dist;
            if v < COLLISION_EPSILON || iter == MAX_ITERATIONS {
                let (u, v, w ) = closest_tri.0.barycentric(Point3::from_vec(closest_dist * closest_n));
                let a = u * closest_tri.1.a + v * closest_tri.1.b + w * closest_tri.1.c;
                return Contact {
                    a: Point3::from_vec(a),
                    b: Point3::from_vec(a - closest_dist*closest_n),
                    n: closest_n,
                    t: 0.0
                };
            }
            // TODO: fix this being a thing.
            let mut to_remove = Vec::new();
            for (i, (a, b, c)) in tris.iter() {
                let n = Triangle::from((a.p, b.p, c.p)).normal();
                if n.dot(support.p - a.p) > 0.0 {
                    edges.add_edge(*a, *b);
                    edges.add_edge(*b, *c);
                    edges.add_edge(*c, *a);
                    to_remove.push(i);
                }
            }
            for i in to_remove {
                tris.remove(i);
            }
            for ([( apx, apy, apz ), ( bpx, bpy, bpz ) ], [ la, lb ]) in edges.map.iter() {
                let a = unsafe {
                    (
                        mem::transmute::<u32, f32>(*apx),
                        mem::transmute::<u32, f32>(*apy),
                        mem::transmute::<u32, f32>(*apz)
                    )
                };
                let b = unsafe {
                    (
                        mem::transmute::<u32, f32>(*bpx),
                        mem::transmute::<u32, f32>(*bpy),
                        mem::transmute::<u32, f32>(*bpz)
                    )
                };
                let a = SupportPoint {
                    p: Point3::new(a.0, a.1, a.2),
                    a: la.0,
                    b: la.1
                };
                let b = SupportPoint {
                    p: Point3::new(b.0, b.1, b.2),
                    a: lb.0,
                    b: lb.1
                };
                tris.push((support, a, b));
            }
            edges.map.clear();
        }
        unreachable!()
    }
}
    
