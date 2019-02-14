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

use cgmath::{EuclideanSpace, InnerSpace, Rotation, Point3, Rotation3, Vector3, Quaternion};

use collision;

/// Maximum tolerence for error, i.e. what we consider the x86 floating
/// point epsilon.
pub const COLLISION_EPSILON: f32 = 0.000001;

/// A normal vector and a distance.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Plane {
    /// The normal vector of the plane.
    pub n: Vector3<f32>,
    /// The plane's distance from the origin.
    pub d: f32,
}

impl Plane {
    /// Create a new plane.
    pub fn new(n: Vector3<f32>, d: f32) -> Self {
        Plane {
            n,
            d
        }
    }
}

impl From<(Point3<f32>, Point3<f32>, Point3<f32>)> for Plane {
    fn from(p: (Point3<f32>, Point3<f32>, Point3<f32>)) -> Self {
        let (a, b, c) = p;
        let n = (b - a).cross(c - a).normalize();
        Plane {
            n: n,
            d: n.dot(a.to_vec()),
        }
    }
}

/// A point and a direction with infinite distance.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Ray {
    /// The origin of the ray.
    pub p: Point3<f32>,
    /// The direction of the ray. Does not need to be normalized.
    pub d: Vector3<f32>,
}

impl Ray {
    /// Create a new ray.
    pub fn new(origin: Point3<f32>, dir: Vector3<f32>) -> Self {
        Ray {
            p: origin,
            d: dir,
        }
    }

    /// Clamp a ray's distance to create a segment.
    pub fn clamp(&self, t: f32) -> Segment {
        Segment {
            a: self.p,
            b: self.p + self.d * t,
        }
    }
}

/// A point and a direction with a finite distance.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Segment {
    /// The starting point of the segment.
    pub a: Point3<f32>,
    /// The end point of the segment.
    pub b: Point3<f32>,
}

impl Segment {
    /// Create a new segment
    pub fn new(from: Point3<f32>, to: Point3<f32>) -> Self {
        Segment {
            a: from,
            b: to 
        }
    }
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Triangle {
    // TODO: There's no point in making these Vectors. Horribly inconsistent.
    /// The first point in the triangle.
    pub a: Vector3<f32>,
    /// The second point in the triangle.
    pub b: Vector3<f32>,
    /// The third point in the triangle.
    pub c: Vector3<f32>,
}

impl Triangle {
    /// Create a new triangle.
    pub fn new(a: Point3<f32>, b: Point3<f32>, c: Point3<f32>) -> Self {
        Triangle {
            a: a.to_vec(),
            b: b.to_vec(),
            c: c.to_vec(),
        }
    }

    /// Calculates the normal vector of a triangle. This is not cached.
    pub fn normal(&self) -> Vector3<f32> {
        (self.b - self.a).cross(self.c - self.a).normalize()
    }

    /// Find the barycentric co-ordinates for a given point.
    pub fn barycentric(&self, p: Point3<f32>) -> (f32, f32, f32) {
        let v0 = self.b - self.a;
        let v1 = self.c - self.a;
        let v2 = p.to_vec() - self.a;
        let d0 = v0.dot(v0);
        let d1 = v0.dot(v1);
        let d2 = v1.dot(v1);
        let d3 = v2.dot(v0);
        let d4 = v2.dot(v1);
        let denom = d0 * d2 - d1 * d1;
        let v = ( d2 * d3 - d1 * d4 ) / denom;
        let w = ( d0 * d4 - d1 * d3 ) / denom;
        ( v, w, 1.0 - v - w )
    }
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

/// Four points in space forming a 3D geometry with four faces.
pub struct Tetrahedron {
    pub a: Vector3<f32>,
    pub b: Vector3<f32>,
    pub c: Vector3<f32>,
    pub d: Vector3<f32>,
}

impl From<(Point3<f32>, Point3<f32>, Point3<f32>, Point3<f32>)> for Tetrahedron {
    fn from(p: (Point3<f32>, Point3<f32>, Point3<f32>, Point3<f32>)) -> Self {
        Tetrahedron {
            a: p.0.to_vec(),
            b: p.1.to_vec(),
            c: p.2.to_vec(),
            d: p.3.to_vec(),
        }
    }
}

/// A center point, two directions, and two half widths.
#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rectangle {
    /// The center of the rectangle.
    pub c: Point3<f32>,
    /// The directions of the rectangle. Must be normalized.
    pub u: [Vector3<f32>; 2],
    /// Half the lengths of each side of the rectangle.
    pub e: [f32; 2],
}

impl Rectangle {
    /// Create a new rectangle.
    pub fn new(center: Point3<f32>, axis: [Vector3<f32>; 2]) -> Self {
        let e = [ axis[0].magnitude(), axis[1].magnitude() ];
        let u = [ axis[0] / e[0], axis[1] / e[1] ];
        Rectangle {
            c: center,
            u,
            e
        }
    }
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
/// AABBs being closed means that a point lying on the surface of the AABB is
/// considered contained by the AABB.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct AABB {
    pub c: Point3<f32>,
    pub r: Vector3<f32>,
}

impl AABB {
    /// Create a new AABB.
    pub fn new(c: Point3<f32>, r: Vector3<f32>) -> Self {
        AABB{ c, r }
    }
}

/// An arbitrarily oriented bounding box.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OBB {
    pub c: Point3<f32>,
    pub q: Quaternion<f32>,
    pub r: Vector3<f32>,
}

impl OBB {
    /// Create a new OBB.
    pub fn new(c: Point3<f32>, r: Vector3<f32>, q: Quaternion<f32>) -> Self {
        OBB{ c, r, q }
    }
}

/// A point and a distance.
///
/// Like AABBs, spheres are closed volumes.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Sphere {
    /// The center point of the sphere.
    pub c: Point3<f32>,
    /// The radius of the sphere.
    pub r: f32,
}

impl Sphere {
    /// Create a new sphere.
    pub fn new(center: Point3<f32>, radius: f32) -> Self {
        assert!(radius > 0.0);
        Sphere {
            c: center,
            r: radius,
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Capsule {
    /// The starting point of the segment sweeping the sphere
    pub a: Point3<f32>,
    /// When added to `a` produces the end point of the segment.
    pub d: Vector3<f32>,
    /// Radius of the sphere.
    pub r: f32,
}

impl Capsule {
    /// Create a new capsule from a segment and a radius.
    pub fn new(segment: Segment, radius: f32) -> Self {
        assert!(radius > 0.0);
        Capsule {
            a: segment.a,
            d: segment.b - segment.a,
            r: radius,
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Moving<T: Shape>(pub T, pub Vector3<f32>);

impl<T: Copy + Clone + Shape> Moving<T> {
    /// Create a moving object with velocity of vel
    pub fn sweep(obj: T, vel: Vector3<f32>) -> Self {
        Moving(obj, vel)
    }

}

impl<T: Shape> Moving<T> {
    /// Return the velocity of the object.
    pub fn vel(&self) -> Vector3<f32> {
        self.1
    }
}

impl<T: Shape> AsRef<T> for Moving<T> {
    fn as_ref(&self) -> &T {
        &self.0
    }
}

impl<T: Shape> AsMut<T> for Moving<T> {
    fn as_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

/// An object with a positional derivative over a time step.
pub trait Delta {
    fn delta(&self) -> Vector3<f32>;
}

impl<T: Shape> Delta for Moving<T> {
    fn delta(&self) -> Vector3<f32> {
        self.1
    }
}


fn clamp(n: f32, min: f32, max: f32) -> f32 {
    if n < min {
        min
    } else if n > max {
        max
    } else {
        n
    }
}

pub fn closest_pts_seg(seg1: &Segment, seg2: &Segment) -> Option<(Point3<f32>, Point3<f32>)> {
    let d1 = seg1.b - seg1.a;
    let d2 = seg2.b - seg2.a;
    let a = d1.magnitude2();
    let e = d2.magnitude2();
    let r = seg1.a - seg2.a;
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
    Some((seg1.a + d1 * s, seg2.a + d2 * t))
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

    /// Returns the closest point on the shape to the given point.
    fn closest_point(&self, to: Point3<f32>) -> Point3<f32>;
}

macro_rules! impl_shape_reqs {
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

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        to + -self.n * (self.n.dot(to.to_vec()) - self.d)
    }
}

impl_shape_reqs!(Ray, p);

impl Shape for Ray {
    fn center(&self) -> Point3<f32> {
        self.p
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let p = (to - self.p).dot(self.d);
        if p < 0.0 {
            self.p
        } else {
            self.p + self.d * (p / self.d.magnitude2())
        }
    }
}

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

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let ab = self.b - self.a;
        let t = ab.dot(to - self.a);
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

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        let ap = to.to_vec() - self.a;
        let d1 = ab.dot(ap);
        let d2 = ac.dot(ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return Point3::from_vec(self.a);
        }

        let bp = to.to_vec() - self.b;
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

        let cp = to.to_vec() - self.c;
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

impl_shape_reqs!(Rectangle, c);

impl Shape for Rectangle {
    fn center(&self) -> Point3<f32> {
        self.c
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let d = to - self.c;
        let mut q = self.c;
        for i in 0..2 {
            let dist = d.dot(self.u[i]);
            q = q + self.u[i] * clamp(dist, -self.e[i], self.e[i]);
        }
        q
    }
}

impl_shape_reqs!(AABB, c);

impl Shape for AABB {    
    fn center(&self) -> Point3<f32> {
        self.c
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        Point3::new(
            clamp(to.x, self.c.x - self.r.x, self.c.x + self.r.x),
            clamp(to.y, self.c.y - self.r.y, self.c.y + self.r.y),
            clamp(to.z, self.c.z - self.r.z, self.c.z + self.r.z),
        )
    }
}

impl_shape_reqs!(OBB, c);

impl Shape for OBB {    
    fn center(&self) -> Point3<f32> {
        self.c
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let to = self.q.invert().rotate_point(to);
        self.q.rotate_point(
            Point3::new(
                clamp(to.x, self.c.x - self.r.x, self.c.x + self.r.x),
                clamp(to.y, self.c.y - self.r.y, self.c.y + self.r.y),
                clamp(to.z, self.c.z - self.r.z, self.c.z + self.r.z),
            )
        )
    }
}

impl_shape_reqs!(Sphere, c);

impl Shape for Sphere {
    fn center(&self) -> Point3<f32> {
        self.c
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let d = to - self.c;
        let rat = d.magnitude2() / (self.r * self.r);
        self.c + d * rat
    }
}

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

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let seg = Segment::from((self.a, self.a + self.d));
        let sphere = Sphere{ c: seg.closest_point(to), r: self.r };
        sphere.closest_point(to)
    }
}

/// A type that is linear in nature.
///
/// Particle types are used to represent the motion of point particles in space.
/// Point particles do not have any spatial extent, i.e. volume.
pub trait Particle: Shape + Copy {
    /// The length of the time step.
    const DT: f32;

    /// The origin of the particle.
    fn pos(&self) -> Point3<f32>;

    /// The direction of the particle.
    fn dir(&self) -> Vector3<f32>;

    /// Rotate this particle around the world space of an object.
    fn rotate_around<R: Rotation3<f32>>(&self, c: Point3<f32>, r: R) -> Self;
}

/// A ray represents a point particle traveling in a direction with infinite
/// velocity.
impl Particle for Ray {
    const DT: f32 = f32::INFINITY;

    fn pos(&self) -> Point3<f32> {
        self.p
    }

    fn dir(&self) -> Vector3<f32> {
        self.d
    }

    fn rotate_around<R: Rotation3<f32>>(&self, c: Point3<f32>, r: R) -> Self {
        let p = r.rotate_vector(self.p - c) + c.to_vec();
        let d = r.rotate_vector(self.d);
        Ray {
            p: Point3::from_vec(p),
            d
        }
    }
}


/// A segment represents a point particle traveling in a direction with finite
/// velocity.
impl Particle for Segment {
    const DT: f32 = 1.0;

    fn pos(&self) -> Point3<f32> {
        self.a
    }

    fn dir(&self) -> Vector3<f32> {
        self.b - self.a
    }
    
    fn rotate_around<R: Rotation3<f32>>(&self, c: Point3<f32>, r: R) -> Self {
        let a = r.rotate_vector(self.a - c) + c.to_vec();
        let d = r.rotate_vector(self.b - self.a);
        let b = a + d;
        Segment {
            a: Point3::from_vec(a),
            b: Point3::from_vec(b)
        }
    }

}

/// A type that is composed of vertices, edges and has a face.
///
/// The property of having a face means that a polygon type can determine if a
/// point is contained on its surface.
pub trait Polygon
    : Copy
    + Shape
    + collision::Contains<Point3<f32>>
    + Into<Plane>
{
    /// The number of vertices and edges available to query.
    const NUM_VERTICES: usize;

    /// Returns the ith vertex as a Point.
    fn vertex(&self, usize) -> Point3<f32>;

    /// Returns the ith edge of the polygon as a pair of indices.
    ///
    /// In the future this is hoped to be deprecate in favor of a constant
    /// associated array. Unfortunately it is impossible as of this moment
    /// to have arrays depend on associated constants.
    fn edge(&self, usize) -> (usize, usize);
}

impl Polygon for Triangle {
    const NUM_VERTICES: usize = 3;

    fn vertex(&self, i: usize) -> Point3<f32> {
        // I really hope this becomes a single copy instead of three, Rust
        // should know how to optimize this.
        Point3::from_vec([self.a, self.b, self.c][i])
    }

    fn edge(&self, i: usize) -> (usize, usize) {
        [(0, 1), (1, 2), (2, 0)][i]
    }
}

impl Polygon for Rectangle {
    const NUM_VERTICES: usize = 4;

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

    fn edge(&self, i: usize) -> (usize, usize) {
        [(0, 1), (1, 2), (2, 3), (3, 0)][i]
    }
}

/// A type that has volume.
///
/// Any type that has Volume can be rotated about its center.
pub trait Volumetric : Shape + Clone {
    /// Rotate the bound in place. This is useless for spheres.
    fn rotate<R: Rotation3<f32>>(self, r: R) -> Self;

    /// Rotates the object around a point.
    fn rotate_about<R: Rotation3<f32>>(mut self, r: R, p: Point3<f32>) -> Self {
        let center = self.center();
        self.set_pos(p + r.rotate_vector(center - p));
        self.rotate(r)
    }
}

impl Volumetric for AABB {
    fn rotate<R: Rotation3<f32>>(self, rot: R) -> AABB {
        let r = self.r;
        let vx = rot.rotate_vector(Vector3::new(r.x, 0.0, 0.0));
        let vy = rot.rotate_vector(Vector3::new(0.0, r.y, 0.0));
        let vz = rot.rotate_vector(Vector3::new(0.0, 0.0, r.z));
        let p1 = self.c + (vx + vy + vz);
        let p2 = self.c + (vx + vy - vz);
        let p3 = self.c + (vx - vy + vz);
        let p4 = self.c + (vx - vy - vz);
        let p5 = self.c + (-vx + vy + vz);
        let p6 = self.c + (-vx + vy - vz);
        let p7 = self.c + (-vx - vy + vz);
        let p8 = self.c + (-vx - vy - vz);
        let lower = Vector3::new(
            p1.x.min(
                p2.x
                    .min(p3.x.min(p4.x.min(p5.x.min(p6.x.min(p7.x.min(p8.x)))))),
            ),
            p1.y.min(
                p2.y
                    .min(p3.y.min(p4.y.min(p5.y.min(p6.y.min(p7.y.min(p8.y)))))),
            ),
            p1.z.min(
                p2.z
                    .min(p3.z.min(p4.z.min(p5.z.min(p6.z.min(p7.z.min(p8.z)))))),
            ),
        );
        let upper = Vector3::new(
            p1.x.max(
                p2.x
                    .max(p3.x.max(p4.x.max(p5.x.max(p6.x.max(p7.x.max(p8.x)))))),
            ),
            p1.y.max(
                p2.y
                    .max(p3.y.max(p4.y.max(p5.y.max(p6.y.max(p7.y.max(p8.y)))))),
            ),
            p1.z.max(
                p2.z
                    .max(p3.z.max(p4.z.max(p5.z.max(p6.z.max(p7.z.max(p8.z)))))),
            ),
        );
        let r = (upper - lower) / 2.0;
        let c = Point3::from_vec((upper + lower) / 2.0);
        AABB { c, r }
    }
}

impl Volumetric for OBB {
    fn rotate<R: Rotation3<f32>>(self, rot: R) -> Self {
        let q: Quaternion<f32> = rot.into();
        OBB {
            c: self.c,
            q: q * self.q,
            r: self.r,
        }
    }
}

impl Volumetric for Sphere {
    /// Rotation to a bounding sphere does nothing.
    fn rotate<R: Rotation3<f32>>(self, _: R) -> Sphere {
        self
    }
}

impl Volumetric for Capsule {
    fn rotate<R: Rotation3<f32>>(self, r: R) -> Self {
        Capsule {
            a: self.center() + r.rotate_vector(self.a - self.center()),
            d: r.rotate_vector(self.d),
            ..self
        }
    }
}

/// A type that represents a convex volume.
pub trait Convex<Support = Point3<f32>>
where
    Support: Into<Point3<f32>>
{
    /// Returns the point on the object that produces the greatest dot product
    /// with the supplied axis. The axis is expected to be normalized.
    // TODO: add an optional "initial" argument so that we can speed up the search.
    fn support(&self, axis: Vector3<f32>) -> Support;
}

impl Convex for AABB {
    fn support(&self, d: Vector3<f32>) -> Point3<f32> {
        Point3::new(
            d.x.signum() * self.r.x,
            d.y.signum() * self.r.y,
            d.z.signum() * self.r.z
        ) + self.c.to_vec()
    }
}

impl Convex for OBB {
    fn support(&self, d: Vector3<f32>) -> Point3<f32> {
        let d = self.q.invert().rotate_vector(d);
        self.q.rotate_point(
            Point3::new(
                d.x.signum() * self.r.x,
                d.y.signum() * self.r.y,
                d.z.signum() * self.r.z
            )
        ) + self.c.to_vec()
    }
}

impl Convex for Sphere {
    fn support(&self, d: Vector3<f32>) -> Point3<f32> {
        self.c + d * self.r
    }
}

impl Convex for Capsule {
    fn support(&self, d: Vector3<f32>) -> Point3<f32> {
        // I have yet to find an implementation of this that I believe to be
        // correct, so for now I'm just going to use the version I found
        // for a cylinder and add the radius to the end.
        use cgmath::Zero;
        let c = self.a + self.d * 0.5;
        let u = self.d.normalize();
        let ud = u.dot(d);
        let w = d - u * ud;
        if w.is_zero() {
            c + (self.d.magnitude() * 0.5 + self.r) * u * ud.signum()
        } else {
            c + (self.d.magnitude() * 0.5 + self.r) * u * ud.signum() + w.normalize() * self.r
        }
    }
}

/// A point that stores the local support points as well as the Minkowski
/// difference.
#[derive(Copy, Clone, Debug)]
pub struct SupportPoint {
    pub p: Point3<f32>,
    pub a: Point3<f32>,
    pub b: Point3<f32>,
}

impl Into<Point3<f32>> for SupportPoint {
    fn into(self) -> Point3<f32> {
        self.p
    }
}

impl From<Point3<f32>> for SupportPoint {
    fn from(p: Point3<f32>) -> Self {
        SupportPoint {
            p,
            a: Point3::new(0.0, 0.0, 0.0),
            b: Point3::new(0.0, 0.0, 0.0),
        }
    }
}

pub struct MinkowskiDiff<'a, 'b, S1, S2>
where
    S1: Convex + 'a,
    S2: Convex + 'b 
{
    pub s1: &'a S1,
    pub s2: &'b S2,
}

impl<'a, 'b, 'c, 'd, S1, S2> Convex for MinkowskiDiff<'a, 'b, S1, S2>
where
    S1: Convex + 'a,
    S2: Convex + 'b,
{
    fn support(&self, axis: Vector3<f32>) -> Point3<f32> {
        Point3::from_vec(self.s1.support(axis) - self.s2.support(-axis))
    }
}

impl<'a, 'b, 'c, 'd, S1, S2> Convex<SupportPoint> for MinkowskiDiff<'a, 'b, S1, S2>
where
    S1: Convex + 'a,
    S2: Convex + 'b,
{
    fn support(&self, axis: Vector3<f32>) -> SupportPoint {
        let a = self.s1.support(axis);
        let b = self.s2.support(-axis);
        let p = Point3::from_vec(a - b);
        SupportPoint {
            p,
            a,
            b
        }
    }
}

/// Computes an orthonormal basis from a given vector. This is usually used to
/// produce tangent vectors for friction contacts.
/// Code taken from http://box2d.org/2014/02/computing-a-basis/
pub fn compute_basis(n: &Vector3<f32>) -> [Vector3<f32>; 2] {
    let b = if n.x.abs() >= 0.57735 {
        Vector3::new(n.y, -n.x, 0.0)
    } else {
        Vector3::new(0.0, n.z, -n.y)
    }.normalize();
    [b, n.cross(b)]
}

#[cfg(test)]
mod tests {
    mod triangle {
        use cgmath::{Point3, EuclideanSpace, InnerSpace};
        use geom::{Triangle, Shape, COLLISION_EPSILON};

        #[test]
        fn test_tri_closest_pt() {
            let tri = Triangle::new(
                Point3::new(2.0, 3.5, 0.0),
                Point3::new(-2.0, -1.5, 0.0),
                Point3::new(2.0, -1.5, 0.0)
            );
            assert!(tri.closest_point(Point3::new(0.0, 0.0, 0.0)).to_vec().magnitude2() < COLLISION_EPSILON);
        }
    }

    mod convex {
        use cgmath::{Vector3, Point3};
        use geom::*;

        #[test]
        fn test_capsule_support_fn() {
            let capsule = Capsule::new(Segment::new(Point3::new(2.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)), 1.0);
            assert_eq!(capsule.support(Vector3::new(0.0, 1.0, 0.0)), Point3::new(5.0, 1.0, 0.0));
            assert_eq!(capsule.support(Vector3::new(-1.0, 0.0, 0.0)), Point3::new(1.0, 0.0, 0.0));
        }
    }
}
