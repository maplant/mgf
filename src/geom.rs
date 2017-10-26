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
use std::ops::{Add, AddAssign, Sub, SubAssign, Neg};

use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Quaternion,  Rotation,
             Vector3, Zero};

use bitset::FixedSizeBitSet;

/// Maximum tolerence for error, i.e. what we consider the x86 floating
/// point epsilon.
pub const COLLISION_EPSILON: f32 = 0.000001;

/// Planes are a normal vector and a distance.
#[derive(Copy, Clone, Debug)]
pub struct Plane {
    pub n: Vector3<f32>,
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

/// Rays are points and direction with infinite distance.
/// Distance vector does not need to be normalized.
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub p: Point3<f32>,
    pub d: Vector3<f32>,
}

/// Segments are Rays with finite distances.
#[derive(Copy, Clone, Debug)]
pub struct Segment {
    pub a: Point3<f32>,
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

/// Triangles are three points in space.
#[derive(Copy, Clone, Debug)]
pub struct Triangle {
    pub a: Vector3<f32>,
    pub b: Vector3<f32>,
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

/// Rectangles are a center point, two directions and two half widths.
#[derive(Copy, Clone)]
pub struct Rectangle {
    pub c: Point3<f32>,
    pub u: [Vector3<f32>; 2],
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

/// Axis Aligned Bounding Boxes are closed boxes aligned to the axes of the
/// coordinate system. AABBs are described by a point and three half widths.
#[derive(Copy, Clone, Debug)]
pub struct AABB {
    pub c: Point3<f32>,
    pub r: Vector3<f32>,
}

/// Spheres are a point and a distance.
/// Like AABBs, spheres as bounds are closed.
#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    pub c: Point3<f32>,
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
#[derive(Copy, Clone, Debug)]
pub struct Capsule {
    pub a: Point3<f32>,  // Line start
    pub d: Vector3<f32>, // Line direction
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

/// Often times we want to determine how close to objects are, or what pair of
/// points on their surfaces are closest.
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

/// Set of points in space.
/// A common goal in collision detection is to determine the minimum distance
/// between two objects.
/// With this in mind we can determine one important thing from any shape: the
/// point contained within the set closest to any given point.
/// Additionally, shapes have some extra properties in that they can be moved
/// around and decompose into a single point (their center).
/// It's not always incredibly efficient to move or take the center point of a
/// shape but it's certainly possible.
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

/// Polygons are objects composed of edges and vertices that can accurately
/// describe whether or not they contain a point.
pub trait Polygon : Shape {
    /// Type of the face object the polygon returns. 
    type FaceType: Collider<Contains, Point3<f32>> + Into<Plane> + Clone;

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
            
/// Determine if objects collide and collect information on the collision if they
/// do. The type of the collision and the amount of information on the collision
/// returned depends on the CollisionType argument requested by the callback passed
/// to collide.
pub trait Collider<CollisionType: Clone, T = Self> {
    /// Collide with an object and call the callback for as many contacts there
    /// are. True is returned if any contact is found.
    fn collide<F: FnMut(CollisionType)>(&self, other: &T, mut callback: F) -> bool {
        if let Some(contact) = self.check_collision(other) {
            callback(contact);
            return true;
        }
        false
    }

    /// Returns the first collision found if any exists.
    fn check_collision(&self, other: &T) -> Option<CollisionType> {
        let mut collision = None;
        self.collide(other, |c| {
            collision = Some(c);
        });
        collision
    }
}

/// Specifies a collision in which two objects overlap. No further information
/// is provided.
#[derive(Copy, Clone, PartialEq)]
pub struct Overlaps;

impl PartialEq<Option<Overlaps>> for Overlaps {
    fn eq(&self, res: &Option<Overlaps>) -> bool {
        res.is_some()
    }
}

/*
macro_rules! commute_overlap {
    (
        $recv:ty, $arg:ty
    ) => {
        impl Collider<Overlaps, $arg> for $recv {
            fn check_collision(&self, : &$arg) -> Option<Overlaps> {
                rhs.check_collision(self)
            }
        }
    };
}
*/

/// Specifies a collision in which the receiver object subsumes the argument
/// object.
#[derive(Copy, Clone, PartialEq)]
pub struct Contains;

impl PartialEq<Option<Contains>> for Contains {
    fn eq(&self, res: &Option<Contains>) -> bool {
        res.is_some()
    }
}

////////////////////////////////////////////////////////////////////////////////
// Point containment:

impl Collider<Contains, Point3<f32>> for Plane {
    fn check_collision(&self, p: &Point3<f32>) -> Option<Contains> {
        if relative_eq!(self.n.dot(p.to_vec()), self.d, epsilon = COLLISION_EPSILON) {
            Option::from(Contains)
        } else {
            None
        }
    }
}

impl Collider<Contains, Point3<f32>> for Triangle {
    fn check_collision(&self, p: &Point3<f32>) -> Option<Contains> {
        let v = p.to_vec() - self.a;
        let ac = self.c - self.a;
        let ab = self.b - self.a;
        let dot1 = ac.dot(ac);
        let dot2 = ac.dot(ab);
        let dot3 = ac.dot(v);
        let dot4 = ab.dot(ab);
        let dot5 = ab.dot(v);
        let invd = 1.0 / (dot1 * dot4 - dot2 * dot2);
        let u = (dot4 * dot3 - dot2 * dot5) * invd;
        let v = (dot1 * dot5 - dot2 * dot3) * invd;
        if u >= 0.0 && v >= 0.0 && (u + v) < 1.0 {
            Contains.into()
        } else {
            None
        }
    }
}

impl Collider<Contains, Point3<f32>> for Rectangle {
    fn check_collision(&self, p: &Point3<f32>) -> Option<Contains> {
        // This might be sped up by simply finding the closest point and checking
        // if it is equal to p. I'm not sure though.
        let p = p.to_vec();
        let n = self.u[0].cross(self.u[1]);
        if relative_eq!(p.dot(n), n.dot(self.c.to_vec()), epsilon = COLLISION_EPSILON)
            && p.dot(self.u[0]).abs() <= self.e[0]
            && p.dot(self.u[1]).abs() <= self.e[1]
        {
            Contains.into()
        } else {
            None
        }
    }
}

impl Collider<Contains, Point3<f32>> for AABB {
    fn check_collision(&self, p: &Point3<f32>) -> Option<Contains> {
        if (self.c.x - p.x).abs() <= self.r.x
            && (self.c.y - p.y).abs() <= self.r.y
            && (self.c.z - p.z).abs() <= self.r.z
        {
            Contains.into()
        } else {
            None
        }
    }
}

impl Collider<Contains, Point3<f32>> for Sphere {
    fn check_collision(&self, p: &Point3<f32>) -> Option<Contains> {
        if (p - self.c).magnitude2() <= self.r*self.r {
            Contains.into()
        } else {
            None
        }
    }
}

/// Intersection models a collision that occurs at some time and has a single
/// point of contact. Useful for collisions that do not concern volume.
#[derive(Copy, Clone, Debug)]
pub struct Intersection {
    pub p: Point3<f32>,
    pub t: f32,
}

macro_rules! commute_intersection {
    (
        $recv:ty, $arg:ty
    ) => {
        impl Collider<Intersection, $arg> for $recv {
             fn collide<F: FnMut(Intersection)>(&self, rhs: &$arg, callback: F) -> bool {
                 rhs.collide(self, callback)
             }
        }
    };
}

impl Collider<Intersection, Ray> for Plane {
    fn collide<F: FnMut(Intersection)>(&self, r: &Ray, mut callback: F) -> bool {
        let denom = self.n.dot(r.d);
        if denom == 0.0 {
            return false;
        }
        let t = (self.d - self.n.dot(r.p.to_vec())) / denom;
        if t <= 0.0 {
            return false;
        }
        callback(Intersection {
            p: r.p + r.d * t,
            t,
        });
        true
    }
}

commute_intersection!(Ray, Plane);

impl Collider<Intersection, AABB> for Ray {
    fn collide<F: FnMut(Intersection)>(&self, a: &AABB, mut callback: F) -> bool {
        let (mut t_min, mut t_max): (f32, f32) = (0.0, f32::INFINITY);
        for dim in 0..3 {
            if self.d[dim].abs() < COLLISION_EPSILON {
                if (self.p[dim] - a.c[dim]).abs() > a.r[dim] {
                    return false;
                }
            } else {
                let ood = 1.0 / self.d[dim];
                let t1 = (a.c[dim] - a.r[dim] - self.p[dim]) * ood;
                let t2 = (a.c[dim] + a.r[dim] - self.p[dim]) * ood;
                if t1 > t2 {
                    t_min = t_min.max(t2);
                    t_max = t_max.min(t1);
                } else {
                    t_min = t_min.max(t1);
                    t_max = t_max.min(t2);
                }
                if t_min > t_max {
                    return false;
                }
            }
        }
        callback(Intersection {
            p: self.p + self.d * t_min,
            t: t_min,
        });
        true
    }
}

impl Collider<Intersection, Sphere> for Ray {
    fn collide<F: FnMut(Intersection)>(&self, s: &Sphere, mut callback: F) -> bool {
        let m = self.p - s.c;
        let a = self.d.dot(self.d);
        let b = m.dot(self.d);
        let c = m.dot(m) - s.r * s.r;
        if c > 0.0 && b > 0.0 {
            return false;
        }
        let discr = b * b - a * c;
        if discr < 0.0 {
            return false;
        }
        let t = ((-b - discr.sqrt()) / a).max(0.0);
        callback(Intersection {
            p: self.p + t * self.d,
            t,
        });
        true
    }
}

impl Collider<Intersection, Capsule> for Ray {
    fn collide<F: FnMut(Intersection)>(&self, cap: &Capsule, mut callback: F) -> bool {
        // This code is horrible and needs to be completely redone.
        let m = self.p - cap.a;
        let md = m.dot(cap.d);
        let nd = self.d.dot(cap.d);
        let dd = cap.d.dot(cap.d);
        let nn = self.d.dot(self.d);
        let mn = m.dot(self.d);
        let a = dd * nn - nd * nd;
        let k = m.dot(m) - cap.r * cap.r;
        if a.abs() < COLLISION_EPSILON {
            let (b, c) = if md < 0.0 {
                (mn, k)
            } else if md > dd {
                let m2 = self.p - (cap.a + cap.d);
                (m2.dot(self.d), m2.dot(m2) - cap.r * cap.r)
            } else {
                // Already colliding
                return false;
            };
            if c > 0.0 && b > 0.0 {
                 return false;
             }
            let discr = b * b - nn*c;
            if discr < 0.0 {
                return false;
            }
            let t = ((-b - discr.sqrt()) / nn).max(0.0);
            callback(Intersection {
                p: self.p + t * self.d,
                t,
            });
            return true;
        }
        let c = dd * k - md * md;
        let b = dd * mn - nd * md;
        let discr = b * b - a * c;
        if discr < 0.0 {
            return false;
        }
        let t = (-b - discr.sqrt()) / a;
        if t < 0.0 {
            // Intersection is behind ray.
            return false;
        }
        // I would like to avoid taking a square root twice, but for now that's
        // the price we pay.
        let t = if md + t * nd < 0.0 {
            if mn > 0.0 && k > 0.0 {
                return false;
            }
            let discr = mn * mn - nn*k;
            if discr < 0.0 {
                return false;
            }
            ((-mn - discr.sqrt()) / nn).max(0.0)
        } else if md + t * nd > dd {
            let m2 = self.p - (cap.a + cap.d);
            let b = m2.dot(self.d);
            let c = m2.magnitude2() - cap.r * cap.r;
            if c > 0.0 && b > 0.0 {
                return false;
            }
            let discr = b * b - nn*c;
            if discr < 0.0 {
                return false;
            }
            ((-b - discr.sqrt()) / nn).max(0.0)
        } else {
            t
        };
        callback(Intersection {
                p: self.p + t * self.d,
                t,
        });
        return true;
    }
}

impl Collider<Intersection, Moving<Sphere>> for Ray {
    fn collide<F: FnMut(Intersection)>(&self, s: &Moving<Sphere>, callback: F) -> bool {
        // ray intersection with moving spheres is exactly the same problem as
        // capsule intersection
        return self.collide(
            &Capsule {
                a: s.0.c,
                d: s.1,
                r: s.0.r,
            },
            callback,
        );
    }
}

/// Anything that can collide with a ray can collide with a line segment.
impl<T: Collider<Intersection, Ray>> Collider<Intersection, T> for Segment {
    fn collide<F: FnMut(Intersection)>(&self, other: &T, mut callback: F) -> bool {
        let mut collision = false;
        let r = Ray::from(*self);
        other.collide(&r, |intersect| if intersect.t <= 1.0 {
            collision = true;
            callback(Intersection {
                t: intersect.t,
                ..intersect
            })
        });
        collision
    }
}

impl Collider<Intersection, Ray> for Triangle {
    fn collide<F: FnMut(Intersection)>(&self, r: &Ray, mut callback: F) -> bool {
        let p = Plane::from(
            (
                Point3::from_vec(self.a),
                Point3::from_vec(self.b),
                Point3::from_vec(self.c),
            )
        );
        let mut collision = false;
        p.collide(r, |i| if Contains == self.check_collision(&i.p) {
            collision = true;
            callback(i);
        });
        collision
    }
}

impl Collider<Intersection, Ray> for Rectangle {
    fn collide<F: FnMut(Intersection)>(&self, r: &Ray, mut callback: F) -> bool {
        let n = self.u[1].cross(self.u[0]);
        let d = n.dot(self.c.to_vec());
        let p = Plane { n: n, d: d };
        let mut collision = false;
        p.collide(r, |i| {
            let d = i.p - self.c;
            let dist1 = d.dot(self.u[0]);
            if dist1 > self.e[0] || dist1 < -self.e[0] {
                return;
            }
            let dist2 = d.dot(self.u[1]);
            if dist2 > self.e[1] || dist2 < -self.e[1] {
                return;
            }
            collision = true;
            callback(i);
        });
        collision
    }
}

commute_intersection!(AABB, Ray);
commute_intersection!(Sphere, Ray);
commute_intersection!(Capsule, Ray);
commute_intersection!(Moving<Sphere>, Ray);

/// Contact models a collision that occurs at some time 0 <= t <= 1.0 for
/// geometries that have a three dimensional surface.
/// Contacts of these types share a convention where the collision normal
/// provided is from the surface of the object implementing the the trait (object A).
#[derive(Copy, Clone, Debug)]
pub struct Contact {
    pub a: Point3<f32>,  // Contact point for collider in global coordinates
    pub b: Point3<f32>,  // Contact point for collidee
    pub n: Vector3<f32>, // Collision normal
    pub t: f32,          // Time of impact
}

impl Contact {
    /// Computes an orthonormal basis for the contact. This is usually used to
    /// produce tangent vectors for friction contacts.
    /// Code taken from http://box2d.org/2014/02/computing-a-basis/
    pub fn compute_basis(&self) -> [Vector3<f32>; 2] {
        let b = if self.n.x.abs() >= 0.57735 {
            Vector3::new(self.n.y, -self.n.x, 0.0)
        } else {
            Vector3::new(0.0, self.n.z, -self.n.y)
        }.normalize();
        [b, self.n.cross(b)]
    }
}

impl Neg for Contact {
    type Output = Contact;

    /// Negate the normal and swap contact points
    fn neg(self) -> Self {
        Contact {
            a: self.b,
            b: self.a,
            n: -self.n,
            ..self
        }
    }
}

macro_rules! commute_contact {
    (
        $recv:ty, $arg:ty
    ) => {
        impl Collider<Contact, $arg> for $recv {
            fn collide<F: FnMut(Contact)>(&self, rhs: &$arg, mut callback: F) -> bool {
                rhs.collide(self, |c: Contact|{ callback(-c) })
            }
        }
    };
}

impl Collider<Contact, Moving<Sphere>> for Plane {
    fn collide<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
        let &Moving(s, v) = sphere;
        let dist = self.n.dot(s.c.to_vec()) - self.d;
        if dist.abs() <= s.r {
            callback(Contact {
                a: s.c + -self.n * dist,
                b: s.c + -self.n * s.r,
                t: 0.0,
                n: self.n,
            });
            return true;
        }
        let denom = self.n.dot(v);
        if denom * dist >= 0.0 {
            return false;
        }
        let r = if dist > 0.0 { s.r } else { -s.r };
        let t = (r - dist) / denom;
        if t <= 1.0 {
            let q = Point3::from_vec(s.c.to_vec() + t * v - r * self.n);
            callback(Contact {
                a: q,
                b: q,
                t: t,
                n: self.n,
            });
            true
        } else {
            false
        }
    }
}

impl Collider<Contact, Moving<Capsule>> for Plane {
    fn collide<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
        // The closest point on a line segment relative to a plane is always
        // either all of the points, if the line is parallel to the plane
        // one of the end points, if the line is not parallel, or
        // a point in the line segment, if the the line is collideing the
        // plane.
        // This does not change if the line segment is translated in any way.
        // Thus, we can use line - plane collision to determine our closest
        // point, put a sphere there, and then perform moving sphere plane
        // collision.
        let &Moving(c, v) = capsule;
        let denom = self.n.dot(c.d.normalize());
        // Line segment collision time
        let center = if denom.abs() < COLLISION_EPSILON {
            // The capsule is parallel to the plane, so we simply choose the
            // mid-point. For planes this contact is always fine since a
            // collision with the center implies a collision with the ends
            // if the capsule is parallel.
            c.a + c.d * 0.5
        } else {
            let t = (self.d - self.n.dot(c.a.to_vec())) / denom;
            // TODO: Fix this
            if t > 1.0 {
                c.a + c.d
            } else if t < 0.0 {
                c.a
            } else {
                // We already collide with the plane.
                let q = c.a + c.d * t;
                callback(Contact {
                    a: q,
                    b: {
                        let dist = self.n.dot(c.a.to_vec()) - self.d;
                        (if dist < 0.0 {
                            c.a
                        } else  {
                            c.a + c.d
                        }) + -self.n * c.r
                    },
                    t: 0.0,
                    n: self.n,
                });
                return true;
            }
        };
        let ms = Moving::sweep(Sphere { c: center, r: c.r }, v);
        self.collide(&ms, callback)
    }
}

commute_contact!(Moving<Sphere>, Plane);
commute_contact!(Moving<Capsule>, Plane);

impl<Poly: Polygon> Collider<Contact, Moving<Sphere>> for Poly {
    fn collide<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
        let &Moving(s, v) = sphere;
        let mut collision = false;
        let p: Plane = self.face().clone().into();
        p.collide(sphere, |contact: Contact| {
            // If the point lies on the face we know automatically that it is
            // a valid collision.
            if Contains == self.face().check_collision(&contact.a) {
                collision = true;
                callback(contact);
                return;
            }
            // Otherwise we need to raycast the sphere's center to the
            // capsules created at the edges of the polygon.
            let mut first_t = f32::INFINITY;
            let mut tri_p = Point3::new(0.0, 0.0, 0.0);
            if v.magnitude2() == 0.0 {
                // Not overlapping and not moving towards each other.
                return;
            }
            let ray = Ray {
                p: s.c,
                d: v,
            };
            for edge_i in 0..self.num_edges() {
                let (a, b) = self.edge(edge_i);
                let v1 = self.vertex(a);
                let v2 = self.vertex(b);
                let c = Capsule{ a: v1, d: v2 - v1, r: s.r };
                ray.collide(&c, |i| {
                    if i.t <= 1.0 && i.t < first_t {
                        first_t = i.t;
                        tri_p = Segment::from((v1, v2)).min_dist(&i.p);
                    }
                });
            }
            if first_t != f32::INFINITY {
                collision = true;
                callback(Contact {
                    a: tri_p,
                    b: tri_p,
                    n: p.n,
                    t: first_t,
                });
            }
        });
        collision
    }
}

commute_contact!(Moving<Sphere>, Triangle);
commute_contact!(Moving<Capsule>, Triangle);
commute_contact!(Moving<Sphere>, Rectangle);
commute_contact!(Moving<Capsule>, Rectangle);

/*
impl<Poly: Polygon + Copy> Collider<Contact, Sphere> for Moving<Poly> {
    fn collide<F: FnMut(Contact)>(&self, poly: &Poly, mut callback: F) -> bool {
        panic!("Fuck")
        // poly.collide(self, |c|callback(-c))
    }
}
*/

// Some 2D segment intersection code for Capsule/Terrain collisions.
fn signed_2d_tri_area(a: Point2<f32>, b: Point2<f32>, c: Point2<f32>) -> f32 {
    (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x)
}

fn seg_2d_intersect(
    a: Point2<f32>,
    b: Point2<f32>,
    c: Point2<f32>,
    d: Point2<f32>,
) -> Option<(Point2<f32>, f32)> {
    let a1 = signed_2d_tri_area(a, b, d);
    let a2 = signed_2d_tri_area(a, b, c);
    if a1 * a2 <= 0.0 {
        let a3 = signed_2d_tri_area(c, d, a);
        let a4 = a3 + a2 - a1;
        if a3 * a4 <= 0.0 {
            let t = a3 / (a3 - a4);
            return Some((a + t * (b - a), t));
        }
    }
    None
}


/// Colliding a moving capsule with a single-faced polygon can potentially
/// produce multiple contacts depending on the angle between the capsule's length
/// and the normal of the polygon's face.
impl<Poly: Polygon> Collider<Contact, Moving<Capsule>> for Poly {
    fn collide<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
        let &Moving(c, v) = capsule;
        let p: Plane = self.face().clone().into();
        // Check if the capsule is already colliding.
        let denom = p.n.dot(c.d.normalize());
        if denom.abs() > COLLISION_EPSILON {
            let t = (p.d - p.n.dot(c.a.to_vec())) / denom;
            if t <= 1.0 && t >= 0.0 {
                // Already colliding the plane
                let q = c.a + c.d * t;
                if Contains == self.face().check_collision(&q) {
                    // Already colliding with triangle
                    callback(Contact {
                        a: q,
                        b: if p.n.dot(c.a.to_vec()) - p.d < 0.0 {
                            c.a
                        } else {
                            c.a + c.d
                        } + -p.n * c.r,
                        t: 0.0,
                        n: p.n,
                    });
                    return true;
                }
            }
        }
        
        // Find a contact to begin our search.
        // The following code needs to be completely re-done. 
        let start_sphere = Moving::sweep(Sphere { c: c.a, r: c.r }, v);
        let end_sphere = Moving::sweep(
            Sphere {
                c: c.a + c.d,
                r: c.r,
            },
            v,
        );

        let found_contact: Option<(Contact, Vector3<f32>)> =
            p.check_collision(&start_sphere)
            .map_or_else(
                || { 
                    p.check_collision(&end_sphere)
                        .map(|c1: Contact| { (c1, -c.d) })
                },
                |c1: Contact| {
                    p.check_collision(&end_sphere)
                        .map_or(Some((c1, c.d)),
                                |c2: Contact| {
                                    if c2.t < c1.t {
                                        Some((c2, -c.d))
                                    } else {
                                        Some((c1, c.d))
                                    }
                                })
                });

        // The following needs to be completely refactored:
        if let Some((contact, dir)) = found_contact {
            // Next we determine which interval of spheres on he capsule
            // could possibly intersect with the triangle.

            // Project the capsule's height onto the plane and use this to
            // form a segment from the contact point. We can use these
            // lines, projected onto the triangle's plane, to determine the
            // minimum and maximum potential spheres to collide with the
            // triangle.
            let silhouette_v = dir - p.n * dir.dot(p.n) / p.n.magnitude2();

            // In order to speed up computation we transform the plane into
            // the x-y plane. This allows us to perform 2D line
            // intersection tests.
            let n_xy = Vector3::new(0.0, 0.0, 1.0);
            let plane_rot = Quaternion::from_arc(p.n, n_xy, None);

            let silhouette_a = Point2::from_vec(
                plane_rot
                    .rotate_vector(contact.a.to_vec() + -p.n * p.d)
                    .truncate(),
            );

            let silhouette_b = Point2::from_vec(
                plane_rot
                    .rotate_vector(contact.a.to_vec() + silhouette_v - p.n * p.d)
                    .truncate(),
            );
            
            // Find the first possible contact
            if Contains == self.face().check_collision(&contact.a) {
                // If the triangle contains the contact we know it is a 
                // valid contact.
                callback(contact);

                // If the capsule is not parallel with the triangle there's 
                // nothing left to do.
                if dir.dot(p.n).abs() >= COLLISION_EPSILON {
                    return true;
                }

                // Find the second contact point on the capsule
                let mut t_max = 0.0;
                for edge_i in 0..self.num_edges() {
                    let (a, b) = self.edge(edge_i);
                    let edge_a = plane_rot
                        .rotate_vector(self.vertex(a).to_vec() -p.n * p.d)
                        .truncate();
                    let edge_b = plane_rot
                        .rotate_vector(self.vertex(b).to_vec() - p.n * p.d)
                        .truncate();
                    if let Some((_, t)) = seg_2d_intersect(
                        silhouette_a, silhouette_b,
                        Point2::from_vec(edge_a), Point2::from_vec(edge_b)
                    ) {
                        if t_max < t {
                            t_max = t
                        }
                    }
                }
                
                let t_max = if t_max == 0.0 { 1.0 } else { t_max };

                // If the capsule is parallel we want to publish a second contact
                // for stability.
                let q = contact.a + silhouette_v * t_max;
                callback(Contact {
                    a: q,
                    b: q,
                    t: contact.t,
                    n: p.n,
                });
                return true;
            }
            if contact.t > 0.0 && dir.dot(p.n).abs() < COLLISION_EPSILON {
                // If the contact time is zero and the contact is not
                // contained in the triangle then we need to revert to 
                // testing the minkowski sum.
                // Otherwise, Determine if and where the capsule's silhouette
                // intersects the silhouette of the triangle.
                let (mut t_min, mut t_max): (f32, f32) = (f32::INFINITY, 0.0);
                let mut found = false;
                for edge_i in 0..self.num_edges() { //edge_base_i..edge_end {
                    let (a, b) = self.edge(edge_i);
                    let edge_a = plane_rot
                        .rotate_vector(self.vertex(a).to_vec() - p.n * p.d)
                        .truncate();
                    let edge_b = plane_rot
                        .rotate_vector(self.vertex(b).to_vec() - p.n * p.d)
                        .truncate();
                    if let Some((_, t)) = seg_2d_intersect(
                        silhouette_a, silhouette_b,
                        Point2::from_vec(edge_a), Point2::from_vec(edge_b)
                    ) {
                        found = true;
                        if t_min > t {
                            t_min = t
                        }
                        if t_max < t {
                            t_max = t
                        }
                    }
                }
                if found {
                    let t_max = if t_max == 0.0 { 1.0 } else { t_max };
                    let q = contact.a + silhouette_v * t_min;
                    let t = contact.t;
                    callback(Contact {
                        a: q,
                        b: q,
                        t,
                        n: p.n,
                    });
                    let q = contact.a + silhouette_v * t_max;
                    callback(Contact {
                        a: q,
                        b: q,
                        t,
                        n: p.n,
                    });
                    return true;
                }
            }
            // If we couldn't find an intersection between the silhouette and
            // the triangle, we need to intersect with the minkowski sum of 
            // the triangle.
        }
        // intersect the minkowski sum of the triangle and the capsule with
        // the ray originating at the capsule's origin.
        // Unfortunately, we have to use an vector here because rust does not
        // support stack allocation if we want to support any number of
        // edges. In order to subvert that as much as possilbe we're using
        // a pool Block64 as a bitset.
        // Thus, any face that wishes to compute here may only have 64
        // vertices (for now).
        // Honestly, that is a rather reasonable requirement. I'm certainly
        // never going to use more than four.
        if self.num_vertices() > 64 {
            return false;
        }
        let mut parallel_edge_vert: u64 = 0;

        // Find all edges that are parallel to the capsule.
        let mut best_par = (
            f32::INFINITY,
            Point3::new(0.0,0.0,0.0),
            Point3::new(0.0,0.0,0.0)
        );
        //        let mut second_contact: Option<Contact> = None;
        for edge_i in 0..self.num_edges() { //edge_base_i..edge_end {
            let (a, b) = self.edge(edge_i);
            let edge_a = self.vertex(a);
            let edge_b = self.vertex(b);
            let ab = edge_b - edge_a;
            let ab_cd = ab.dot(c.d);
            if ab_cd.abs() != c.d.magnitude()*ab.magnitude() {
                // Not parallel
                continue;
            }
            // Capsule is parallel to the triangle edge.
            parallel_edge_vert.insert(a);
            parallel_edge_vert.insert(b);
            // capsule is parallel to edge
            // Form a capsule at the edge
            let ray = Ray { p: c.a, d: v };
            let (edge_a, edge_b) = if ab_cd < 0.0 {
                (edge_b, edge_a)
            } else {
                (edge_a, edge_b)
            };
            let edge_sum = Capsule{ a: edge_a, d: edge_b - edge_a, r: c.r };
            // Here is a good example of where using both the bool
            // returned and the callback is nice. If we found a
            // collision, but it didn't beat the best time, we
            // still want to skip checking the next edge.
            let m_edge = ab.magnitude2();
            if !ray.collide(&edge_sum, |inter| {
                if inter.t > best_par.0.min(1.0) {
                    return;
                }
                let tri_p: Point3<f32> =
                    Segment::from((edge_a, edge_b)).min_dist(&inter.p);

                let m_proj = ((tri_p + c.d) - edge_a).magnitude2();
                let c_t = if  m_proj > m_edge {
                    (m_proj - m_edge)
                        / (m_proj - (tri_p - edge_a).magnitude2())
                } else {
                    1.0
                };
                let q = tri_p + c.d * c_t;
                best_par = (inter.t, tri_p, q);

            }) {
                let vert_sum = Capsule{ a: edge_a, d: -c.d, r: c.r };
                ray.collide(&vert_sum, |inter| {
                    if inter.t > best_par.0.min(1.0) {
                        return;
                    }
                    let d = inter.p - edge_a;
                    let capsule_t = -d.dot(c.d) / c.d.magnitude2();
                    let tri_p: Point3<f32> =
                        Segment::from((edge_a, edge_a + -c.d))
                        .min_dist(&inter.p);
                    let a = tri_p + c.d * capsule_t;
                    let m_proj = ((tri_p + c.d) - edge_a).magnitude2();
                    let b = if m_proj > m_edge {
                        edge_b
                    } else {
                        tri_p + c.d
                    };
                    best_par = (inter.t, a, b);
                });
            }
        }
        // Perform edge collisions.
        let mut best_sum = (
            f32::INFINITY,
            Point3::new(0.0,0.0,0.0),
        );

        for edge_i in 0..self.num_edges() { 
            let (a, b) = self.edge(edge_i);

            // If this edge is parallel, skip it.
            let a_on_parallel_edge = parallel_edge_vert.get(a);
            let b_on_parallel_edge = parallel_edge_vert.get(b);
            if a_on_parallel_edge && b_on_parallel_edge {
                continue;
            }
            let edge_a = self.vertex(a);
            let edge_b = self.vertex(b);
            let tris = [
                Triangle::from((edge_a + -c.d, edge_a, edge_b)),
                Triangle::from((edge_a + -c.d, edge_b, edge_b + -c.d)),
            ];
            let p: Plane = tris[1].into();
            let s = Sphere {
                c: c.a,
                r: c.r
            };
            // We need to perform a moving sphere/tri collision. We
            // optimize it slightly by only doing relevant capsule
            // ray casts.
            p.collide(&Moving::sweep(s, v), |contact: Contact| {
                if (Contains == tris[0].check_collision(&contact.a)
                    || Contains == tris[1].check_collision(&contact.a))
                    && best_sum.0 > contact.t
                {
                    let d = contact.a - edge_a;
                    let capsule_t = -d.dot(c.d) / c.d.magnitude2();
                    best_sum = (
                        contact.t,
                        contact.a + c.d * capsule_t
                    );
                } else {
                    let ray = Ray { p: c.a, d: v };
                    // We always need to check the top and bottom edges of
                    // the quad made with the tri edge and the capsule's
                    // direction
                    let bottom_edge = Capsule {
                        a: edge_a,
                        d: edge_b - edge_a,
                        r: c.r
                    };
                    ray.collide(&bottom_edge, |inter| {
                        if inter.t > 1.0 || inter.t > best_sum.0 {
                            return;
                        }
                        let q: Point3<f32> = Segment::from((edge_a, edge_b))
                            .min_dist(&inter.p);
                        best_sum = (inter.t, q);
                    });
                    let top_edge = Capsule {
                        a: edge_a + -c.d,
                        d: edge_b - edge_a,
                        r: c.r
                    };
                    ray.collide(&top_edge, |inter| {
                        if inter.t > 1.0 || inter.t > best_sum.0 {
                            return;
                        }
                        let plane_p = inter.p + c.d;
                        let q: Point3<f32> = Segment::from((edge_a, edge_b))
                            .min_dist(&plane_p);
                        best_sum = (inter.t, q);
                    });
                    // Hopefully this gets unrolled.
                    for &(vert, is_parallel) in [
                        (edge_a, a_on_parallel_edge),
                        (edge_b, b_on_parallel_edge),
                    ].iter() {
                        if is_parallel {
                            continue;
                        }
                        let cap = Capsule{ a: vert, d: -c.d, r: c.r };
                        ray.collide(&cap, |inter| {
                            if inter.t > 1.0 || inter.t >= best_sum.0 {
                                return;
                            }
                            best_sum = (inter.t, vert);
                        });
                    }

                }
            });
        }
        if best_sum.0 < best_par.0 {
            callback(Contact {
                t: best_sum.0,
                a: best_sum.1,
                b: best_sum.1,
                n: p.n,
            });
        } else if best_par.0 != f32::INFINITY {
            callback(Contact {
                t: best_par.0,
                a: best_par.1,
                b: best_par.1,
                n: p.n,
            });
            callback(Contact {
                t: best_par.0,
                a: best_par.2,
                b: best_par.2,
                n: p.n,
            });
        } else {
            return false;
        }
        return true;
    }
}

impl Collider<Contact, Moving<Sphere>> for Sphere {
    fn collide<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
        let &Moving(s, v) = sphere;
        // Check if spheres are already overlapping
        let r = self.r + s.r;
        let d = s.c.to_vec() - self.c.to_vec();
        let len = d.magnitude2();
        if len <= r * r {
            let n = if len == 0.0 {
                if v.is_zero() {
                    return false;
                } else {
                    -v.normalize()
                }
            } else {
                d / len.sqrt()
            };
            callback(Contact {
                a: self.c + n * self.r,
                b: s.c + -n * s.r,
                t: 0.0,
                n,
            });
            return true;
        }
        let l = v.magnitude2();
        if l == 0.0 {
            // Not overlapping and not moving towards each other.
            return false;
        }
        let ray = Ray {
            p: self.c,
            d: -v,
        };
        ray.collide(&Sphere{ c: s.c, r }, |intersect| if intersect.t <= 1.0 {
            let end_c = s.c + v * intersect.t;
            let ba = (end_c.to_vec() - self.c.to_vec()).normalize();
            let a = self.c + ba * self.r;
            callback(Contact {
                a,
                b: a,
                t: intersect.t,
                // The collision normal for two spheres is simply the
                // surface normal at the point of contact for either.
                n: ba,
            })
        })
    }
}

commute_contact!(Sphere, Moving<Capsule>);

impl Collider<Contact, Moving<Sphere>> for Capsule {
    fn collide<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
        // Collision for moving capsule and sphere is exactly the same as ray
        // capsule intersection except we must increase the radius of the capsule
        // by the radius of the sphere.
        let &Moving(s, v) = sphere;
        // Check if capsule and sphere are already overlapping
        let r = self.r + s.r;
        let closest_pt: Point3<f32> =
            Segment::from((self.a, self.a + self.d)).min_dist(&s.c);
        let d = s.c - closest_pt;
        let len = d.magnitude2();
        if len <= r * r {
            let n = if len == 0.0 {
                if v.is_zero() {
                    return false;
                } else {
                    -v.normalize()
                }
            } else {
                d / len.sqrt()
            };
            callback(Contact {
                a: closest_pt + n * self.r,
                b: s.c + -n * s.r,
                t: 0.0,
                n,
            });
            return true;
        }
        let l = v.magnitude2();
        if l == 0.0 {
            // Not overlapping and not moving towards each other.
            return false;
        }
        let ray = Ray {
            p: s.c,
            d: v,
        };
        ray.collide(&Capsule{ r: s.r + self.r, ..*self }, |intersect| if intersect.t <= 1.0 {
            let b = s.c + v * intersect.t;
            let a: Point3<f32> = Segment::from(*self).min_dist(&b);
            let ba = (b - a).normalize();
            let q = a + ba * self.r;
            callback(Contact {
                a: q,
                b: q,
                t: intersect.t,
                // The collision normal for two spheres is simply the
                // surface normal at the point of contact for either.
                n: ba,
            })
        })
    }
}

impl Collider<Contact, Moving<Capsule>> for Capsule {
    fn collide<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
        // TODO: reduce all of the normalizations and magnitudes
        let &Moving(c, v) = capsule;
        let self_seg: Segment = (self.a, self.a + self.d).into();
        let (p1, p2) = if let Some((p,_)) = 
            self_seg.min_dist(&Segment::from((c.a, c.a + v)))
        {
            (p, self_seg.min_dist(&Segment::from((c.a + c.d, c.a + c.d + v))).unwrap().0)
        } else {
            (self.a, self.a + self.d)
        };
        let self_seg: Segment = (p1, p2).into();
        if let Some((q, _)) = self_seg.min_dist(&Segment::from((c.a, c.a + c.d))) {
            let ss = Sphere {
                c: q,
                r: self.r,
            };
            return ss.collide(capsule, callback);
        }

        // The capsules are parallel, we need to use the power of
        // intervals to solve this problem
        // Project c.a - self.a and c.a + c.d - self.a onto self.d.
        let d_mag2 = self.d.magnitude2();
        let t1 = (c.a - self.a).dot(self.d) / d_mag2;
        let t2 = (c.a + c.d - self.a).dot(self.d) / d_mag2;
        let (t_min, t_max, c_a, c_d) =
            if t1 < t2 { (t1, t2, c.a, c.d) } else { (t2, t1, c.a + c.d, -c.d) };

        // Determine the height from the line formed by the two capsules.
        let h = self.a - (c_a + c_d * (-t_min / (t_max - t_min)));
        let h_len = h.magnitude();

        if h_len <= self.r + c.r {
            // If the height between the two capsules is less than the sum of
            // their radii, then either the capsules were already colliding or
            // they will only collide at the ends.
            if t_max <= 0.0 { 
                return self.collide(
                    &Moving::sweep(Sphere{ c: c_a + c_d, r: c.r }, v), callback
                );
            }
            if t_min >= 1.0 {
                return self.collide(
                    &Moving::sweep(Sphere { c: c_a, r: c.r }, v), callback
                );
            }
            return false;
        }

        // Invariant: h_len > self.r + c.r => h_rat > 0.0
        let h_rat = (h_len - self.r - c.r) / h_len;

        // Transform the moving capsule's interval to when the height between
        // the two capsules is zero.
        
        // Project v onto h
        let v_comp = v.dot(h) / (h_len * h_len);
        if v_comp < h_rat {
            // The capsules do not collide
            return false;
        }
        // Invariant: v_comp > h_rat, therefore v_comp > 0.0
        let coll_t = h_rat / v_comp;
        let v_travel = v * coll_t;
        // Project distance traveled before collision on self.d
        let axis_t_delta = v_travel.dot(self.d) / d_mag2;
        // Transform the times to that at the collision and determine collision
        // points
        let t_min = t_min + axis_t_delta;
        let t_max = t_max + axis_t_delta;

        // We need to check one more time to determine if the ends of the
        // capsules are colliding:
        if t_max <= 0.0 { 
            return self.collide(
                &Moving::sweep(Sphere{ c: c_a + c_d, r: c.r }, v), callback
            );
        }
        if t_min >= 1.0 {
            return self.collide(
                &Moving::sweep(Sphere { c: c_a, r: c.r }, v), callback
            );
        }

        // Invariants:
        // t_max is greater than or equal to zero
        // t_min is less than or equal to one
        // t_max is greater than t_min

        // If we have arrived here, two parallel capsules are colliding at their
        // sides. Handling this is very simple: clamp t_min and t_max to the
        // interval [0.0, 1.0] and average the two.
        let s_t = (clamp(t_min, 0.0, 1.0)
                   + clamp(t_max, 0.0, 1.0)) * 0.5;
        let o_t = (s_t - t_min) / (t_max - t_min);
        let a_c = self.a + self.d * s_t;
        let b_c = c_a + c_d * o_t + v_travel;
        let ab = b_c - a_c;
        let n = if ab.is_zero() {
            // If the distance between the two object's centers is zero, rely on
            // velocity. If velocity is zero, no contact can be generated.
            if v.is_zero() {
                return false;
            } else {
                -v.normalize()
            }
        } else {
            (b_c - a_c).normalize()
        };
        callback(Contact {
            a: a_c + n * self.r,
            b: b_c + -n * c.r,
            t: coll_t,
            n,
        });
        return true;
    }
}

//commute_contact!(Moving<Sphere>, Sphere);
//commute_contact!(Moving<Sphere>, Capsule);
//commute_contact!(Moving<Capsule>, Capsule);

impl<Recv, Arg> Collider<Contact, Arg> for Moving<Recv>
where
    Arg: Shape + Copy,
    Recv: Collider<Contact, Moving<Arg>> + Shape + Copy
{
    fn collide<F: FnMut(Contact)>(&self, other: &Arg, mut callback: F) -> bool {
        let other_moving = Moving::sweep(*other, -self.1);
        self.0.collide(&other_moving, |c: Contact| {
            let d = self.1 * c.t;
            let a = c.a + d;
            let b = c.b + d;
            callback(Contact{ a, b, ..c })
        })
    }
}

/// Any two moving object collision can be reduced to a one moving one static.
/// object collision. This is done by finding the relative velocity between the
/// two objects.
impl<K, T> Collider<Contact, Moving<T>> for Moving<K>
where
    K: Collider<Contact, Moving<T>> + Copy + Clone + Shape,
    T: Copy + Clone + Shape,
{
    fn collide<F: FnMut(Contact)>(&self, other: &Moving<T>, mut callback: F) -> bool {
        let Moving(oo, ov) = *other;
        let Moving(so, sv) = *self;
        let rel_v = ov - sv;
        let rel_obj = Moving(oo, rel_v);
        so.collide(&rel_obj, |c: Contact| {
            let a = c.a + sv * c.t;
            let b = c.b + sv * c.t;
            callback(Contact { a, b, ..c })
        })
    }
}

/// A local contact stores a contact plus the contact points relative to the
/// space of the objects colliding.
#[derive(Copy, Clone, Debug)]
pub struct LocalContact {
    pub local_a: Point3<f32>,
    pub local_b: Point3<f32>,
    pub global: Contact,
}

impl Neg for LocalContact {
    type Output = LocalContact;

    /// Negate the normal and swap contact points
    fn neg(self) -> Self {
        LocalContact {
            local_a: self.local_b,
            local_b: self.local_a,
            global: -self.global,
        }
    }
}

impl<Recv, Arg> Collider<LocalContact, Arg> for Recv
where
    Recv: Shape + Collider<Contact, Arg> + Delta,
    Arg: Shape + Delta
{
    fn collide<F: FnMut(LocalContact)>(&self, other: &Arg, mut callback: F) -> bool {
        self.collide(other, |c: Contact| {
            let a_c = self.center() + self.delta() * c.t;
            let b_c = other.center() + other.delta() * c.t;
            callback(LocalContact {
                local_a: c.a + -a_c.to_vec(),
                local_b: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}

impl<Recv, Arg> Collider<LocalContact, Moving<Arg>> for Recv
where
    Recv: Shape + Collider<Contact, Moving<Arg>>,
    Arg: Shape + Copy 
{
    fn collide<F: FnMut(LocalContact)>(&self, other: &Moving<Arg>, mut callback: F) -> bool {
        self.collide(other, |c: Contact| {
            let a_c = self.center();
            let b_c = other.as_ref().center() + other.vel() * c.t;
            callback(LocalContact {
                local_a: c.a + -a_c.to_vec(),
                local_b: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}

impl<Recv, Arg> Collider<LocalContact, Arg> for Moving<Recv>
where
    Recv: Shape + Copy,
    Arg: Shape + Collider<Contact, Moving<Recv>> + Copy
{
    fn collide<F: FnMut(LocalContact)>(&self, other: &Arg, mut callback: F) -> bool {
        other.collide(self, |c: Contact| {
            let a_c = self.as_ref().center() + self.vel() * c.t;
            let b_c = other.center();
            callback(LocalContact {
                local_a: c.b + -a_c.to_vec(),
                local_b: c.a + -b_c.to_vec(),
                global: -c
            })
        })
    }
}

impl<Recv, Arg> Collider<LocalContact, Moving<Arg>> for Moving<Recv>
where
    Recv: Shape + Collider<Contact, Moving<Arg>> + Copy,
    Arg: Shape + Copy
{
    fn collide<F: FnMut(LocalContact)>(&self, other: &Moving<Arg>, mut callback: F) -> bool {
        // Need to take a further look at this
        let Moving(oo, ov) = *other;
        let Moving(so, sv) = *self;
        let rel_v = ov - sv;
        let rel_obj = Moving(oo, rel_v);
        so.collide(&rel_obj, |c: Contact| {
            let a = c.a + sv * c.t;
            let b = c.b + sv * c.t;
            let local_a = a + -so.center().to_vec();
            let local_b = b + -(other.as_ref().center() + ov * c.t).to_vec();
            callback(LocalContact {
                local_a,
                local_b,
                global: Contact {
                    a, b,
                    ..c
                }
            })
        }) 
    }
}

#[cfg(test)]
mod tests {
    mod spheres {
        use cgmath::{Point3, Vector3};
        use geom;
        use geom::{Collider, Contact, Moving, Sphere, Rect, Triangle};

        #[test]
        fn test_moving_spheres_collision() {
            let s1 = Moving::sweep(
                Sphere {
                    c: Point3::new(-3.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let s2 = Moving::sweep(
                Sphere {
                    c: Point3::new(3.0, 0.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(-2.0, 0.0, 0.0),
            );
            
            let collision: Contact = s1.check_collision(&s2).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(-1.0, 0.0, 0.0));
            assert_eq!(collision.b, Point3::new(-1.0, 0.0, 0.0));
        }

        #[test]
        fn test_rect_collision() {
            let floor = Rect {
                c: Point3::new(0.0, 1.0, 0.0),
                u: [Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)],
                e: [3.0, 3.0],
            };
            let sphere_collide_center = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 0.0),
            );
            assert!(floor.collide(&sphere_collide_center, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            assert!(sphere_collide_center.collide(&floor, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, -1.0, 0.0));
            }));
            let sphere_collide_center_2s = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -20.0, 0.0),
            );
            assert!(floor.collide(&sphere_collide_center_2s, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.t, 0.5);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            let sphere_collide_corner = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 3.0),
            );
            assert!(floor.collide(&sphere_collide_corner, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 3.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 3.0));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            let sphere_miss_corner = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 3.00001),
            );
            assert!(!floor.collide(&sphere_miss_corner, |_:Contact| {}));
        }

        #[test]
        fn test_tri_collision() {
            let floor = Triangle {
                a: Vector3::new(1.0, 1.0, 0.0),
                b: Vector3::new(0.0, 1.0, 1.0),
                c: Vector3::new(0.0, 1.0, -1.0),
            };
            let sphere_collide_center = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 0.0),
            );
            assert!(floor.collide(&sphere_collide_center, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            let sphere_collide_corner = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 1.0),
            );
            assert!(floor.collide(&sphere_collide_corner, |c: Contact| {
                assert_relative_eq!(c.a, Point3::new(0.0, 1.0, 1.0), epsilon = geom::COLLISION_EPSILON);
                assert_relative_eq!(c.b, Point3::new(0.0, 1.0, 1.0), epsilon = geom::COLLISION_EPSILON);
                assert!((1.0 - c.t) < geom::COLLISION_EPSILON);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            let sphere_miss_corner = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 1.00001),
            );
            assert!(!floor.collide(&sphere_miss_corner, |c: Contact| { panic!("C = {:?}", c); }));
            let sphere_collide_edge = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.5, -10.0, 0.5),
            );
            assert!(floor.collide(&sphere_collide_edge, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.5, 1.0, 0.5));
                assert_eq!(c.b, Point3::new(0.5, 1.0, 0.5));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
        }
    }

    mod capsules {
        use cgmath::{Point3, Vector3};
        use geom;
        use geom::{Capsule, Sphere, Collider, Contact, Moving, Triangle, Rect};

        #[test]
        fn test_moving_sphere_collision() {
            let c = Capsule {
                a: Point3::new(4.0, 3.0, 5.5),
                d: Vector3::new(0.0, 1.0, 0.0),
                r: 2.0,
            };
            let s = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 3.0, 5.5),
                    r: 1.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let collision: Contact = c.check_collision(&s).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.0, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.0, 5.5));
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.0, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.0, 5.5));
 
        }

        #[test]
        fn test_moving_capsule_collision() {
            let s = Capsule {
                a: Point3::new(4.0, 3.0, 5.5),
                d: Vector3::new(0.0, 1.0, 0.0),
                r: 2.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 3.0, 5.5),
                    d: Vector3::new(0.0, 1.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.5, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.5, 5.5));
            let s = Capsule {
                a: Point3::new(4.0, 3.0, 5.5),
                d: Vector3::new(0.0, 1.0, 0.0),
                r: 1.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 3.0, 5.5),
                    d: Vector3::new(0.0, 1.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.a, Point3::new(3.0, 3.5, 5.5));
            assert_eq!(collision.b, Point3::new(3.0, 3.5, 5.5));
            assert_eq!(collision.t, 1.0);
            let s = Capsule {
                a: Point3::new(1.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(-2.0, 0.0, 0.0),
                    d: Vector3::new(-1.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(2.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.a, Point3::new(0.0, 0.0, 0.0));
            assert_eq!(collision.b, Point3::new(0.0, 0.0, 0.0));
            assert_eq!(collision.t, 0.5);
            let s = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 0.0, 0.0),
                    d: Vector3::new(-1.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(2.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.a, Point3::new(-1.0, 0.0, 0.0));
            assert_eq!(collision.b, Point3::new(1.0, 0.0, 0.0));
            assert_eq!(collision.t, 0.0);
            let s = Capsule {
                a: Point3::new(4.0, 3.0, 5.5),
                d: Vector3::new(0.0, 1.0, 0.0),
                r: 2.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 2.0, 5.5),
                    d: Vector3::new(0.0, 1.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.0, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.0, 5.5));
            let s = Capsule {
                a: Point3::new(4.0, 3.0, 5.5),
                d: Vector3::new(0.0, 1.0, 0.0),
                r: 2.0,
            };
            let c = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 2.5, 5.5),
                    d: Vector3::new(0.0, 1.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(1.0, 0.0, 0.0),
            );
            let collision: Contact = s.check_collision(&c).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.25, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.25, 5.5));
        }

        #[test]
        fn test_rect_collision() {
            let floor = Rect {
                c: Point3::new(0.0, 1.0, 0.0),
                u: [Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)],
                e: [3.0, 3.0],
            };
            let capsule_level_off_center = Moving::sweep(
                Capsule {
                    a: Point3::new(1.0, 13.0, 0.0),
                    d: Vector3::new(3.0, 0.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 0.0),
            );
            let mut contacts: Vec<Contact> = Vec::new();
            floor.collide(&capsule_level_off_center, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(3.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
        }

        #[test]
        fn test_tri_collision() {
            let floor = Triangle {
                a: Vector3::new(1.0, 1.0, 0.0),
                b: Vector3::new(0.0, 1.0, 1.0),
                c: Vector3::new(0.0, 1.0, -1.0),
            };
            let capsule_clip_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(0.9, 3.0, 1.0),
                    d: Vector3::new(0.0, 0.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            let mut contacts: Vec<Contact> = Vec::new();
            floor.collide(&capsule_clip_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.9, 1.0, 0.1), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(0.9, 1.0, -0.1), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_clip_off_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.9, 3.0, 0.0),
                    d: Vector3::new(0.0, 0.0, 2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_clip_off_center, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.9, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(0.9, 1.0, 0.1), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_clip_off_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.9, 3.0, 0.0),
                    d: Vector3::new(0.0, 0.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_clip_off_center, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.9, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(0.9, 1.0, -0.1), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_through_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.9, 2.0, 0.0),
                    d: Vector3::new(1.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_through_center, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 0.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.9, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_tilted_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.5, 4.0, 0.0),
                    d: Vector3::new(-1.0, -0.5, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -2.0, 0.0),
            );

            let collision: Contact = floor.check_collision(&capsule_tilted_center).unwrap();
            // Wolfram alpha gives us an answer of 0.8149827, which frankly I am
            // fine iwth considering how dubious my methodology is.
            assert_eq!(collision.t, 0.81598306);
            assert_relative_eq!(collision.a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            /*
            // Ray intersection stability must be improved!.
            let capsule_tilted_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.5, 4.0, 0.0),
                    d: Vector3::new(1.0, -1.5, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -2.0, 0.0),
            );
            let collision: Contact = floor.check_collision(&capsule_tilted_center).unwrap();
            assert_eq!(collision.t, 0.75);
            assert_relative_eq!(collision.a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            */
            let capsule_tilted_center = Moving::sweep(
                Capsule {
                    a: Point3::new(0.5, 4.0, 0.0),
                    d: Vector3::new(-1.0, -1.0, 2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -2.0, 0.0),
            );
            let collision: Contact = floor.check_collision(&capsule_tilted_center).unwrap();
            assert_relative_eq!(collision.a, Point3::new(0.0, 1.0, 1.0), epsilon =  geom::COLLISION_EPSILON);
            assert_eq!(collision.t, 0.7022774);
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 2.0, 2.0),
                    d: Vector3::new(0.0, 0.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 1.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 2);
            assert_relative_eq!(contacts[1].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 4.0, 2.0),
                    d: Vector3::new(0.0, -2.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 4.0, 0.0),
                    d: Vector3::new(0.0, 2.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 2.0, 2.0),
                    d: Vector3::new(0.0, 0.0, -4.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 1.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 2);
            assert_relative_eq!(contacts[1].a, Point3::new(0.0, 1.0, -1.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 2.0, -2.0),
                    d: Vector3::new(0.0, 0.0, 4.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, -1.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 2);
            assert_relative_eq!(contacts[1].a, Point3::new(0.0, 1.0, 1.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let floor = Triangle {
                a: Vector3::new(1.0, 1.0, 0.0),
                b: Vector3::new(0.0, 1.0, 2.0),
                c: Vector3::new(0.0, 1.0, -2.0),
            };
            let capsule_parallel_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-0.5, 2.0, 0.5),
                    d: Vector3::new(0.0, 0.0, -1.0),
                    r: 0.5,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_parallel_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.5), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 2);
            assert_relative_eq!(contacts[1].a, Point3::new(0.0, 1.0, -0.5), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
            let capsule_perp_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 2.0, 0.0),
                    d: Vector3::new(-3.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_perp_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_perp_to_edge = Moving::sweep(
                Capsule {
                    a: Point3::new(-4.0, 2.0, 0.0),
                    d: Vector3::new(3.0, 0.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_perp_to_edge, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_next_to_vert = Moving::sweep(
                Capsule {
                    a: Point3::new(2.0, 2.0, 1.0),
                    d: Vector3::new(0.0, 0.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_next_to_vert, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_next_to_vert_skewed = Moving::sweep(
                Capsule {
                    a: Point3::new(2.0, 2.0, 1.0),
                    d: Vector3::new(0.0, -1.0, -2.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_next_to_vert_skewed, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 0.5);
            assert_relative_eq!(contacts[0].a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_intersects_tri_plane = Moving::sweep(
                Capsule {
                    a: Point3::new(0.0, 4.0, 0.0),
                    d: Vector3::new(-2.0, -4.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_intersects_tri_plane, |c| contacts.push(c));
            assert_relative_eq!(contacts[0].t, 0.7639319);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
            let capsule_intersects_tri_plane = Moving::sweep(
                Capsule {
                    a: Point3::new(-1.0, 2.0, 0.0),
                    d: Vector3::new(-1.0, -2.0, 0.0),
                    r: 1.0,
                },
                Vector3::new(0.0, -1.0, 0.0),
            );
            floor.collide(&capsule_intersects_tri_plane, |c| contacts.push(c));
            assert_relative_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
        }
    }

    mod hitscans {
        use cgmath::{InnerSpace, Point3, Vector3};
        use geom;
        use geom::{Capsule, Collider, Ray};

        #[test]
        fn test_ray_collisions() {
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(1.0, -3.0, 0.0),
                d: Vector3::new(-0.25, 1.0, 0.0).normalize(),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_relative_eq!(intersection.p, Point3::new(0.5, -1.0, 0.0),
                                epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(r.p + r.d*intersection.t, Point3::new(0.5, -1.0, 0.0),
                                epsilon = geom::COLLISION_EPSILON);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(0.0, -3.0, 0.0),
                d: Vector3::new(0.25, 1.0, 0.0).normalize(),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_relative_eq!(intersection.p, Point3::new(0.5, -1.0, 0.0),
                                epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(r.p + r.d*intersection.t, Point3::new(0.5, -1.0, 0.0),
                                epsilon = geom::COLLISION_EPSILON);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(0.0, 2.0, 0.0),
                r: 2.0,
            };
            let r = Ray {
                p: Point3::new(4.0, 1.0, 0.0),
                d: Vector3::new(-1.0, 0.0, 0.0),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_eq!(intersection.p, Point3::new(2.0, 1.0, 0.0));
            assert_eq!(intersection.t, 2.0);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(3.0, 0.0, 0.0),
                d: Vector3::new(-1.0, 0.0, 0.0),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_eq!(intersection.p, Point3::new(2.0, 0.0, 0.0));
            assert_eq!(intersection.t, 1.0);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(-2.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_eq!(intersection.p, Point3::new(-1.0, 0.0, 0.0));
            assert_eq!(intersection.t, 1.0);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(-2.0, 0.5, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
            };
            let intersection = r.check_collision(&c).unwrap();
            // Funny enough, I did these calculations by myself to ensure the
            // code was correct, and it panic'd because sure enough I had done
            // the math wrong for the correct time of impact. 
            assert_relative_eq!(intersection.p, Point3::new(-0.8660254037844386, 0.5, 0.0));
            assert_relative_eq!(intersection.t, 1.13397459621556196,
                                epsilon = geom::COLLISION_EPSILON);
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(3.0, 0.5, 0.0),
                d: Vector3::new(-1.0, 0.0, 0.0),
            };
            let intersection = r.check_collision(&c).unwrap();
            assert_relative_eq!(intersection.p, Point3::new(1.8660254037844386, 0.5, 0.0));
            assert_relative_eq!(intersection.t, 1.13397459621556196,
                                epsilon = geom::COLLISION_EPSILON);
        }
    }
}
