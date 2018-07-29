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

use cgmath::prelude::*;
use cgmath::{EuclideanSpace, InnerSpace, Point3, Vector3};
use geom::*;

pub struct Simplex<Support = Point3<f32>>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    points: [Support; 4],
    state: &'static SimplexState<Support>,
}

impl<Support> fmt::Debug for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static + fmt::Debug
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

impl<Support> From<[Support; 1]> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static 
{
    fn from(p: [Support; 1]) -> Self {
        Self::from(p[0])
    }
}
          
impl<Support> From<Support> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: Support) -> Self {
        Simplex {
            points: [
                p,
                Support::from(Point3::new(0.0, 0.0, 0.0)),
                Support::from(Point3::new(0.0, 0.0, 0.0)),
                Support::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &VERTEX_DATAPTRLOC,
        }
    }
}

impl<Support> From<[Support; 2]> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Support; 2]) -> Self {
        Self::from((p[0], p[1]))
    }
}

impl<Support> From<(Support, Support)> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Support, Support)) -> Self {
        Simplex {
            points: [
                p.0,
                p.1,
                Support::from(Point3::new(0.0, 0.0, 0.0)),
                Support::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &EDGE_DATAPTRLOC,
        }
    }
}

impl<Support> From<[Support; 3]> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Support; 3]) -> Self {
        Self::from((p[0], p[1], p[2]))
    }
}
          
impl<Support> From<(Support, Support, Support)> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Support, Support, Support)) -> Self {
        Simplex {
            points: [
                p.0,
                p.1,
                p.2,
                Support::from(Point3::new(0.0, 0.0, 0.0)),
            ],
            state: &FACE_DATAPTRLOC,
        }
    }
}

impl<Support> From<[Support; 4]> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: [Support; 4]) -> Self {
        Self::from((p[0], p[1], p[2], p[3]))
    }
}

impl<Support> From<(Support, Support, Support, Support)> for Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    fn from(p: (Support, Support, Support, Support)) -> Self {
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

impl<Support> Simplex<Support>
where
    Support: Into<Point3<f32>> + From<Point3<f32>> + Copy + Clone + 'static
{
    /// Finds the closest point on the shape to the origin
    pub fn closest_point_to_origin<S>(&mut self, shape: &S) -> Point3<f32>
    where
        S: Convex<Support>
    {
        loop {
            let (min_norm, next_state) = self.state.min_norm(&mut self.points);
            self.state = next_state;
            if min_norm.magnitude2() < COLLISION_EPSILON {
                return Point3::from_vec(min_norm);
            }
            let support = shape.support(-min_norm.normalize());
            let support_v = support.into().to_vec();
            if min_norm.magnitude2() == support_v.magnitude2() {
                return Point3::from_vec(min_norm);
            }
            self.state.add_point(&mut self.points, support);
        }
    }
}

trait SimplexState<Support>
where
    Support: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, &mut [Support; 4]) -> (Vector3<f32>, &'static SimplexState<Support>);

    fn add_point(&self, &mut [Support; 4], Support);
}

struct VertexSimplex{}
struct EdgeSimplex{}
struct FaceSimplex{}
struct VolumeSimplex{}

static VERTEX_DATAPTRLOC: VertexSimplex = VertexSimplex{};
static EDGE_DATAPTRLOC: EdgeSimplex = EdgeSimplex{};
static FACE_DATAPTRLOC: FaceSimplex = FaceSimplex{};
static VOLUME_DATAPTRLOC: VolumeSimplex = VolumeSimplex{};

impl<Support> SimplexState<Support> for VertexSimplex
where
    Support: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Support; 4]) -> (Vector3<f32>, &'static SimplexState<Support>) {
        (simp[0].into().to_vec(), &EDGE_DATAPTRLOC)
    }

    fn add_point(&self, simp:  &mut [Support; 4], p: Support) {
        simp[0] = p
    }
}

impl<Support> SimplexState<Support> for EdgeSimplex
where
    Support: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Support; 4]) -> (Vector3<f32>, &'static SimplexState<Support>) {
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

    fn add_point(&self, simp: &mut [Support; 4], p: Support) {
        simp[1] = p
    }
}    


impl<Support> SimplexState<Support> for FaceSimplex
where
    Support: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Support; 4]) -> (Vector3<f32>, &'static SimplexState<Support>) {
        let ab = simp[1].into() - simp[0].into();
        let ac = simp[2].into() - simp[0].into();
        let ap = -simp[0].into().to_vec();
        let d1 = ab.dot(ab);
        let d2 = ac.dot(ap);

        // Vertex region A
        if d1 <= 0.0 && d2 <= 0.0 {
            return (simp[0].into().to_vec(), &EDGE_DATAPTRLOC);
        }

        // Vertex region B
        let bp = -simp[1].into().to_vec();
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
        let cp = -simp[2].into().to_vec();
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

    fn add_point(&self, simp:  &mut [Support; 4], p: Support) {
        simp[2] = p;
    }
}

fn origin_outside_plane(
    a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>, d: Vector3<f32>
) -> bool {
    let ab_x_ac = (b - a).cross(c - a);
    let sign_p = (-a).dot(ab_x_ac);
    let sign_d = (d - a).dot(ab_x_ac);
    sign_p * sign_d < 0.0
}

impl<Support> SimplexState<Support> for VolumeSimplex
where
    Support: Into<Point3<f32>> + Copy + Clone + 'static
{
    fn min_norm(&self, simp: &mut [Support; 4]) -> (Vector3<f32>, &'static SimplexState<Support>) {
        let mut closest_pt: Vector3<f32> = Vector3::zero();
        let mut best_dist: f32 = f32::INFINITY;
        let mut next_state: &'static SimplexState<Support> = &VERTEX_DATAPTRLOC;
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

    fn add_point(&self, simp: &mut [Support; 4], p: Support) {
        simp[3] = p
    }
}
