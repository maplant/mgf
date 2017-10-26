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
use std::ops::{Add, AddAssign, Sub, SubAssign};

use cgmath::{EuclideanSpace, InnerSpace, Point3, Quaternion, Rotation, Vector3};

use collision;

/// Maximum tolerence for error, i.e. what we consider the x86 floating
/// point epsilon.
pub const COLLISION_EPSILON: f32 = 0.000001;

/// A normal vector and a distance.
#[derive(Copy, Clone, Debug)]
pub struct Plane {
    /// The normal vector of the plane.
    pub n: Vector3<f32>,
    /// The plane's distance from the origin.
    pub d: f32,
}

impl From<(Point3<f32>, Point3<f32>, Point3<f32>)> for Plane {
    fn from(p: (Point3<f32>, Point3<f32>, Point3<f32>)) -> Self {
        let (a, b, c) = p;
        let n = (b - a).cross(b - c).normalize();
        Plane {
            n: n,
            d: n.dot(a.to_vec()),
        }
    }
}

/// A point and a direction with infinite distance.
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    /// The origin of the ray.
    pub p: Point3<f32>,
    /// The direction of the ray. Does not need to be normalized.
    pub d: Vector3<f32>,
}

/// A point and a direction with a finite distance.
#[derive(Copy, Clone, Debug)]
pub struct Segment {
    /// The starting point of the segment.
    pub a: Point3<f32>,
    /// The end point of the segment.
    pub b: Point3<f32>,
}

impl From<(Point3<f32>, Point3<f32>)> for Segment {
    fn from(p: (Point3<f32>, Point3<f32>)) -> Self {
        Segment {
            a: p.0,
            b: p.1
        }
    }
}

impl From<Segment> for Ray {
    fn from(s: Segment) -> Self {
        Ray {
            p: s.a,
            d: s.b - s.a,
        }
    }
}

/// Three points in space.
#[derive(Copy, Clone, Debug)]
pub struct Triangle {
    /// The first point in the triangle.
    pub a: Vector3<f32>,
    /// The second point in the triangle.
    pub b: Vector3<f32>,
    /// The third point in the triangle.
    pub c: Vector3<f32>,
}

// TODO: in the future, check that the triangle is clockwise and orient it if it
// is not.
impl From<(Point3<f32>, Point3<f32>, Point3<f32>)> for Triangle {
    fn from(p: (Point3<f32>, Point3<f32>, Point3<f32>)) -> Self {
        Triangle {
            a: p.0.to_vec(),
            b: p.1.to_vec(),
            c: p.2.to_vec(),
        }
    }
}

impl From<Triangle> for Plane {
    fn from(t: Triangle) -> Self {
        Plane::from(
            (
                Point3::from_vec(t.a),
                Point3::from_vec(t.b),
                Point3::from_vec(t.c)
            )
        )
    }
}

/// A center point, two direction, and two half widths.
#[derive(Copy, Clone)]
pub struct Rectangle {
    /// The center of the rectangle.
    pub c: Point3<f32>,
    /// The directions of the rectangle.
    pub u: [Vector3<f32>; 2],
    /// Half the lengths of each side of the rectangle.
    pub e: [f32; 2],
}

pub type Rect = Rectangle;

impl Into<Plane> for Rectangle {
    fn into(self) -> Plane {
        let n = self.u[1].cross(self.u[0]);
        let d = n.dot(self.c.to_vec());
        Plane{ n, d }
    }
}

/// An Axis Aligned Bounding Box.
///
/// AABBs are closed boxes aligned to the axes of the coordinate system. AABBs
/// are described by a point and three half widths.
///
/// AABB's being Closed means that a point lying on the surface of the AABB is
/// considered contained by the AABB.
#[derive(Copy, Clone, Debug)]
pub struct AABB {
    pub c: Point3<f32>,
    pub r: Vector3<f32>,
}

/// A point and a distance.
///
/// Like AABBs, spheres as closed volumes.
#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    /// The center point of the sphere.
    pub c: Point3<f32>,
    /// The radius of the sphere.
    pub r: f32,
}

impl Sphere {
    /// Rotates the sphere about the origin
    pub fn rotate(self, r: Quaternion<f32>) -> Self {
        Sphere {
            c: r.rotate_point(self.c),
            ..self
        }
    }
}

/// A sphere swept along a line.
///
/// Capsules are described as all of the spheres that have a center point
/// along their line segment.
///
/// A capsule can be constructed from a `Moving<Sphere>`
#[derive(Copy, Clone, Debug)]
pub struct Capsule {
    /// The starting point of the segment sweeping the sphere
    pub a: Point3<f32>,
    /// When added to `a` produces the end point of the segment.
    pub d: Vector3<f32>,
    /// Radius of the sphere.
    pub r: f32,
}

impl Capsule {
    /// Rotates the capsule about the origin
    pub fn rotate(self, r: Quaternion<f32>) -> Self {
        Capsule {
            a: r.rotate_point(self.a),
            d: r.rotate_vector(self.d),
            ..self
        }
    }
}

impl From<Capsule> for Segment {
    fn from(c: Capsule) -> Segment {
        // Starting to regret this choice of segment
        Segment{ a: c.a, b: c.a + c.d }
    }
}

impl From<Moving<Sphere>> for Capsule {
    fn from(m: Moving<Sphere>) -> Self {
        Capsule {
            a: m.0.c,
            d: m.1,
            r: m.0.r,
        }
    }
}

/// A geometry swept accross a given path of motion.
#[derive(Copy, Clone, Debug)]
pub struct Moving<T: Copy + Clone + Shape>(pub T, pub Vector3<f32>);

impl<T: Copy + Clone + Shape> Moving<T> {
    /// Create a moving object with velocity of vel
    pub fn sweep(obj: T, vel: Vector3<f32>) -> Self {
        Moving(obj, vel)
    }

    /// Return the velocity of the object.
    pub fn vel(&self) -> Vector3<f32> {
        self.1
    }
}

impl<T: Copy + Clone + Shape> AsRef<T> for Moving<T> {
    fn as_ref(&self) -> &T {
        &self.0
    }
}

/// An object with a positional derivative over a timestep
pub trait Delta {
    fn delta(&self) -> Vector3<f32>;
}

impl<T: Copy + Clone + Shape> Delta for Moving<T> {
    fn delta(&self) -> Vector3<f32> {
        self.1
    }
}


#[inline(always)]
fn clamp(n: f32, min: f32, max: f32) -> f32 {
    if n < min {
        min
    } else if n > max {
        max
    } else {
        n
    }
}

/// Finds the result corresponding to the minimum distance between two objects.
pub trait MinDistance<To = Point3<f32>, Result = Point3<f32>> {
    fn min_dist(&self, &To) -> Result;
}

impl MinDistance<Point3<f32>, f32> for Plane {
    /// Returns the smallest distance squared
    fn min_dist(&self, q: &Point3<f32>) -> f32 {
        self.n.dot(q.to_vec()) - self.d
    }
}

impl MinDistance<Point3<f32>> for Plane {
    /// Returns closest point on plane to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        q + -self.n * (self.n.dot(q.to_vec()) - self.d)
    }
}

impl MinDistance<Point3<f32>> for Ray {
    /// Returns closest point on ray to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        let p = (q - self.p).dot(self.d);
        if p < 0.0 {
            self.p
        } else {
            self.p + self.d * (p / self.d.magnitude2())
        }
    }
}

impl MinDistance<Point3<f32>, f32> for Segment {
    /// Returns squared distance between the segment and the point.
    fn min_dist(&self, c: &Point3<f32>) -> f32 {
        let ab = self.b - self.a;
        let ac = c - self.a;
        let e = ac.dot(ab);
        if e <= 0.0 {
            ac.magnitude2()
        } else {
            let f = ab.magnitude2();
            if e >= f {
                let bc = c - self.a;
                bc.magnitude2()
            } else {
                ac.magnitude2() - e * e / f
            }
        }
    }
}

impl MinDistance<Point3<f32>> for Segment {
    /// Returns closest point on segment to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        let ab = self.b - self.a;
        let t = ab.dot(q - self.a);
        if t <= 0.0 {
            self.a
        } else {
            let denom = ab.dot(ab);
            if t >= denom {
                self.b
            } else {
                self.a + ab * (t / denom)
            }
        }
    }
}

impl MinDistance<Segment, Option<(Point3<f32>, Point3<f32>)>> for Segment {
    /// Returns the a pair of points that are the minimum distance from each
    /// other. If the segments are parallel, return None.
    fn min_dist(&self, to: &Segment) -> Option<(Point3<f32>, Point3<f32>)> {
        let d1 = self.b - self.a;
        let d2 = to.b - to.a;
        let a = d1.magnitude2();
        let e = d2.magnitude2();
        let r = self.a - to.a;
        let f = d2.dot(r);
        let (s, t) = if a <= COLLISION_EPSILON {
            if e <= COLLISION_EPSILON {
                (0.5, 0.5)
            } else {
                (0.5, clamp(f / e, 0.0, 1.0))
            }
        } else {
            let c = d1.dot(r);
            if e <= COLLISION_EPSILON {
                (clamp(-c / a, 0.0, 1.0), 0.0)
            } else {
                let b = d1.dot(d2);
                let denom = a * e - b * b;
                let s = if denom != 0.0 {
                    clamp((b * f - c * e) / denom, 0.0, 1.0)
                } else {
                    return None;
                };
                let t = b * s + f;
                if t < 0.0 {
                    (clamp(-c / a, 0.0, 1.0), 0.0)
                } else if t > e {
                    (clamp((b - c) / a, 0.0, 1.0), 1.0)
                } else {
                    (s, t / e)
                }
            }
        };
        Some((self.a + d1 * s, to.a + d2 * t))
    }
}

impl MinDistance<Point3<f32>> for Triangle {
    /// Returns closest point on the triangle to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        let ap = q.to_vec() - self.a;
        let d1 = ab.dot(ab);
        let d2 = ac.dot(ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return Point3::from_vec(self.a);
        }

        let bp = q.to_vec() - self.b;
        let d3 = ab.dot(bp);
        let d4 = ac.dot(bp);
        if d3 >= 0.0 && d4 <= d3 {
                return Point3::from_vec(self.b);
        }

        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return Point3::from_vec(self.a + ab * v);
        }

        let cp = q.to_vec() - self.c;
        let d5 = ab.dot(cp);
        let d6 = ac.dot(cp);
        if d6 >= 0.0 && d5 <= d6 {
            return Point3::from_vec(self.c);
        }
        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return Point3::from_vec(self.a + ac * w);
        }

        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return Point3::from_vec(self.b + (self.c - self.b) * w);
        }

        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        Point3::from_vec(self.a + ab * v + ac * w)
    }
}

impl MinDistance<Point3<f32>> for Rectangle {
    /// Returns closest point on the rectangle to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        let d = q - self.c;
        let mut q = self.c;
        for i in 0..2 {
            let dist = d.dot(self.u[i]);
            q = q + self.u[i] * clamp(dist, -self.e[i], self.e[i]);
        }
        q
    }
}

impl MinDistance<Point3<f32>> for AABB {
    /// Returns closest point on the AABB to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        Point3::new(
            clamp(q.x, self.c.x - self.r.x, self.c.x + self.r.x),
            clamp(q.y, self.c.y - self.r.y, self.c.y + self.r.y),
            clamp(q.z, self.c.z - self.r.z, self.c.z + self.r.z),
        )
    }
}

impl MinDistance<Point3<f32>> for Sphere {
    /// Returns closest point on the sphere to q
    fn min_dist(&self, q: &Point3<f32>) -> Point3<f32> {
        let d = q - self.c;
        let rat = d.magnitude2() / (self.r * self.r);
        self.c + d * rat
    }
}

impl MinDistance<Point3<f32>> for Capsule {
    /// Returns closest point on the capsule to q
    fn min_dist(&self, _q: &Point3<f32>) -> Point3<f32> {
        unimplemented!();
    }
}

/// A type that describes a set of points in space.
///
/// Shape is usually used to describe static static objects such as a regular
/// `Sphere` or `Capsule`. It allows for these shapes to be decomposed into
/// some arbitrary center point and displaced by Vectors.
pub trait Shape
    : AddAssign<Vector3<f32>>
    + SubAssign<Vector3<f32>>
{
    /// Returns the center of mass of the geometry, assuming a regular density.
    fn center(&self) -> Point3<f32>;

    /// Sets the center of the shape to p.
    fn set_pos(&mut self, p: Point3<f32>) {
        let disp = p - self.center();
        *self += disp;
    }
}

macro_rules! impl_shape {
    (
        $name:ident, $center:ident
    ) => {
        impl Add<Vector3<f32>> for $name {
            type Output = Self;

            fn add(self, v: Vector3<f32>) -> Self {
                $name{ $center: self.$center + v, ..self }
            }
        }
        impl Sub<Vector3<f32>> for $name {
            type Output = Self;

            fn sub(self, v: Vector3<f32>) -> Self {
                $name{ $center: self.$center + -v, ..self }
            }
        }
        impl AddAssign<Vector3<f32>> for $name {
            fn add_assign(&mut self, v: Vector3<f32>) {
                self.$center += v
            }
        }
        impl SubAssign<Vector3<f32>> for $name {
            fn sub_assign(&mut self, v: Vector3<f32>) {
                self.$center += -v
            }
        }
        impl Shape for $name {
            fn center(&self) -> Point3<f32> {
                self.$center
            }
        }
    };
}

impl Add<Vector3<f32>> for Plane {
    type Output = Self;

    fn add(self, v: Vector3<f32>) -> Plane {
        Plane{ d: (self.n*self.d + v).dot(self.n), ..self }
    }
}

impl Sub<Vector3<f32>> for Plane {
    type Output = Self;

    fn sub(self, v: Vector3<f32>) -> Plane {
        Plane{ d: (self.n*self.d - v).dot(self.n), ..self }
    }
}

impl AddAssign<Vector3<f32>> for Plane {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.d = (self.n*self.d + v).dot(self.n)
    }
}

impl SubAssign<Vector3<f32>> for Plane {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.d = (self.n*self.d - v).dot(self.n)
    }
}

impl Shape for Plane {
    /// Conceptually, it's hard to consider any point of a plane the "center".
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.n * self.d)
    }
}

impl_shape!(Ray, p);

impl Add<Vector3<f32>> for Segment {
    type Output = Self;

    fn add(self, v: Vector3<f32>) -> Segment {
        Segment{ a: self.a + v, b: self.b + v }
    }
}

impl Sub<Vector3<f32>> for Segment {
    type Output = Self;

    fn sub(self, v: Vector3<f32>) -> Segment {
        Segment{ a: self.a + -v, b: self.b + -v }
    }
}

impl AddAssign<Vector3<f32>> for Segment {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.a += v;
        self.b += v;
    }
}

impl SubAssign<Vector3<f32>> for Segment {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.a += -v;
        self.b += -v;
    }
}

impl Shape for Segment {
    fn center(&self) -> Point3<f32> {
        self.a + (self.b-self.a) * 0.5
    }
}

impl Add<Vector3<f32>> for Triangle {
    type Output = Self;

    fn add(self, v: Vector3<f32>) -> Triangle {
        Triangle{ a: self.a + v, b: self.b + v, c: self.c + v }
    }
}

impl Sub<Vector3<f32>> for Triangle {
    type Output = Self;

    fn sub(self, v: Vector3<f32>) -> Triangle {
        Triangle{ a: self.a - v, b: self.b - v, c: self.c - v }
    }
}

impl AddAssign<Vector3<f32>> for Triangle {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.a += v;
        self.b += v;
        self.c += v;
    }
}

impl SubAssign<Vector3<f32>> for Triangle {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.a += -v;
        self.b += -v;
        self.c += -v;
    }
}

impl Shape for Triangle {
    fn center(&self) -> Point3<f32> {
        Point3::from_vec((self.a + self.b + self.c) / 3.0)
    }
}

impl_shape!(Rectangle, c);
impl_shape!(AABB, c);
impl_shape!(Sphere, c);

impl Add<Vector3<f32>> for Capsule {
    type Output = Self;

    fn add(self, v: Vector3<f32>) -> Capsule {
        Capsule{ a: self.a + v, ..self }
    }
}

impl Sub<Vector3<f32>> for Capsule {
    type Output = Self;

    fn sub(self, v: Vector3<f32>) -> Capsule {
        Capsule{ a: self.a + -v, ..self }
    }
}

impl AddAssign<Vector3<f32>> for Capsule {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.a += v;
    }
}

impl SubAssign<Vector3<f32>> for Capsule {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.a += -v;
    }
}

impl Shape for Capsule {
    fn center(&self) -> Point3<f32> {
        self.a + self.d * 0.5
    }
}

/// A type that is composed of vertices, edges and has a face.
///
/// Besides being able to present their vertices and edges, a Polygon can
/// produce a reference to some a face, which is a type that can be
/// decomposed into a plane and determine if it contains a point.
pub trait Polygon : Shape {
    /// Type of the face object the polygon returns. 
    type FaceType: collision::Contains<Point3<f32>> + Into<Plane> + Clone;

    /// The number of edges available to query.
    fn num_vertices(&self) -> usize;

    /// Returns the ith vertex as a Point.
    fn vertex(&self, usize) -> Point3<f32>;

    /// The number of edges available to query.
    fn num_edges(&self) -> usize;

    /// Returns the ith edge of the polygon as a pair of indices.
    fn edge(&self, usize) -> (usize, usize);

    fn face(&self) -> &Self::FaceType;
}

impl Polygon for Triangle {
    type FaceType = Self;

    fn num_vertices(&self) -> usize { 3 }

    fn vertex(&self, i: usize) -> Point3<f32> {
        // I really hope this becomes a single copy instead of three, Rust
        // should know how to optimize this.
        Point3::from_vec([self.a, self.b, self.c][i])
    }

    fn num_edges(&self) -> usize { 3 }

    fn edge(&self, i: usize) -> (usize, usize) {
        [(0, 1), (1, 2), (2, 0)][i]
    }

    fn face(&self) -> &Triangle { self }
}

impl Polygon for Rectangle {
    type FaceType = Self;

    fn num_vertices(&self) -> usize { 4 }

    fn vertex(&self, i: usize) -> Point3<f32> {
        match i {
            // top right (1, 1)
            0 => self.c + self.u[0]*self.e[0] + self.u[1]*self.e[1],
            // bottom right (1, -1)
            1 => self.c + self.u[0]*self.e[0] + -self.u[1]*self.e[1],
            // bottom left (-1, -1)
            2 => self.c + -self.u[0]*self.e[0] + -self.u[1]*self.e[1],
            // top left (1, 1)
            3 => self.c + -self.u[0]*self.e[0] + self.u[1]*self.e[1],
            _ => unreachable!(),
        }
    }

    fn num_edges(&self) -> usize { 4 }

    fn edge(&self, i: usize) -> (usize, usize) {
        [(0, 1), (1, 2), (2, 3), (3, 0)][i]
    }

    fn face(&self) -> &Rectangle { self }
}
