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

use std::ops::{Add, Sub, Mul, Div};
use cgmath::{EuclideanSpace, InnerSpace, Point3, Vector3, Zero};

use geom::*;
use collision::*;

/// A type that can overlap, contain, and be combined with one another.
///
/// Bounds can also be scaled and extended by scalar values to produce proxy
/// volumes. They can also be rotated and have a surface area.
/// 
/// Bounds are Shapes and thus have a center and can be displaced.
pub trait Bound
    : Copy
    + Add<Vector3<f32>, Output = Self>
    + Sub<Vector3<f32>, Output = Self>
    + Mul<f32, Output = Self>           // Scalar multiply
    + Div<f32, Output = Self>           // Scalar divide
    + Add<f32, Output = Self>           // Scalar extend
    + Sub<f32, Output = Self>           // Scalar shrink
    + Volumetric
    + Overlaps<Self>
    + Contains<Self>
{
    /// Produce a bound that encloses the two arguments.
    fn combine(&Self, &Self) -> Self;

    /// Return a measure of the area of the object.
    fn surface_area(&self) -> f32;
}

/// A type that can be decomposed into a bound. 
pub trait BoundedBy<B: Bound> {
    fn bounds(&self) -> B;
}

/// All geometries that satisfy Bound are bounded by themselves.
impl<B: Bound> BoundedBy<B> for B {
    #[inline(always)]
    fn bounds(&self) -> B {
        *self
    }
}

impl<B: Bound, T: Copy + Clone + Shape +BoundedBy<B>> BoundedBy<B> for Moving<T> {
    /// The bounds for a swept object is the bounds extended in the direction
    /// and magnitude of the velocity.
    fn bounds(&self) -> B {
        let s_bounds = self.0.bounds();
        let e_bounds = s_bounds + self.1;
        B::combine(&s_bounds, &e_bounds)
    }
}

////////////////////////////////////////////////////////////////////////////////
// AABBs satisfy Bound

impl Mul<f32> for AABB {
    type Output = Self;

    /// Expand AABB
    fn mul(self, s: f32) -> AABB {
        AABB{ r: self.r * s, ..self }
    }
}

impl Div<f32> for AABB {
    type Output = Self;

    /// Compress AABB
    fn div(self, s: f32) -> Self {
        AABB{ r: self.r / s, ..self }
    }
}

impl Add<f32> for AABB {
    type Output = Self;

    /// Extend AABB
    fn add(self, s: f32) -> AABB {
        AABB{ r: self.r + Vector3::new(s, s, s), ..self }
    }
}

impl Sub<f32> for AABB {
    type Output = Self;

    /// Shrink AABB
    fn sub(self, s: f32) -> Self {
        AABB{ r: self.r - Vector3::new(s, s, s), ..self }
    }
}

impl Bound for AABB {
    /// The returned AABB is the smallest volume possible that encloses both
    /// arguments. At least I think, I Can't remember if that claim is true.
    /// I'll have to check sometime.
    fn combine(a: &AABB, b: &AABB) -> AABB {
        let lower = Vector3::new(
            (a.c.x - a.r.x).min(b.c.x - b.r.x),
            (a.c.y - a.r.y).min(b.c.y - b.r.y),
            (a.c.z - a.r.z).min(b.c.z - b.r.z),
        );
        let upper = Vector3::new(
            (a.c.x + a.r.x).max(b.c.x + b.r.x),
            (a.c.y + a.r.y).max(b.c.y + b.r.y),
            (a.c.z + a.r.z).max(b.c.z + b.r.z),
        );
        let r = (upper - lower) / 2.0;
        assert!(r.x >= 0.0);
        assert!(r.y >= 0.0);
        assert!(r.z >= 0.0);
        let c = Point3::from_vec((upper + lower) / 2.0);
        AABB { c, r }
    }

    fn surface_area(&self) -> f32 {
        self.r.x * self.r.y + self.r.y * self.r.z + self.r.z * self.r.x
    }
}

impl BoundedBy<AABB> for Triangle {
    fn bounds(&self) -> AABB {
        let c = (self.a + self.b + self.c) / 3.0;
        let d0 = (self.a.x - c.x)
            .abs()
            .max((self.b.x - c.x).abs().max((self.c.x - c.x).abs()));
        let d1 = (self.a.y - c.y)
            .abs()
            .max((self.b.y - c.y).abs().max((self.c.y - c.y).abs()));
        let d2 = (self.a.z - c.z)
            .abs()
            .max((self.b.z - c.z).abs().max((self.c.z - c.z).abs()));
        AABB {
            c: Point3::from_vec(c),
            r: Vector3::new(d0, d1, d2),
        }
    }
}

impl BoundedBy<AABB> for Rectangle {
    fn bounds(&self) -> AABB {
        let p1 = self.c + self.u[0] * self.e[0];
        let p2 = self.c + self.u[1] * self.e[1];
        let d0 = (p1.x - self.c.x).abs().max((p2.x - self.c.x).abs());
        let d1 = (p1.y - self.c.y).abs().max((p2.y - self.c.y).abs());
        let d2 = (p1.z - self.c.z).abs().max((p2.z - self.c.z).abs());
        AABB {
            c: self.c,
            r: Vector3::new(d0, d1, d2),
        }
    }
}

impl BoundedBy<AABB> for Sphere {
    fn bounds(&self) -> AABB {
        AABB {
            c: self.c,
            r: Vector3::new(self.r, self.r, self.r),
        }
    }
}

impl BoundedBy<AABB> for Capsule {
    fn bounds(&self) -> AABB {
        // Include length of d with radius to cover all rotations
        let r = self.r + self.d.magnitude() * 0.5;
        AABB {
            c: self.a + self.d*0.5,
            r: Vector3::new(r, r, r),
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
// Spheres satisfy Bound

impl Mul<f32> for Sphere {
    type Output = Self;

    fn mul(self, s: f32) -> Self {
        Sphere{ r: s * self.r, ..self }
    }
}

impl Div<f32> for Sphere {
    type Output = Self;

    fn div(self, s: f32) -> Self {
        Sphere{ r: s / self.r, ..self }
    }
}

impl Add<f32> for Sphere {
    type Output = Self;

    fn add(self, s: f32) -> Self {
        Sphere{ r: s + self.r, ..self }
    }
}

impl Sub<f32> for Sphere {
    type Output = Self;

    fn sub(self, s: f32) -> Self {
        Sphere{ r: self.r - s, ..self }
    }
}

impl Bound for Sphere {
    fn combine(a: &Sphere, b: &Sphere) -> Sphere {
        let d = b.c - a.c;
        let r = b.r - a.r;
        if r * r >= d.magnitude2() {
            if a.r >= b.r {
                *a
            } else {
                *b
            }
        } else {
            let dist = d.magnitude();
            let r = (dist + a.r + b.r) * 0.5;
            Sphere {
                c: a.c + if dist > COLLISION_EPSILON {
                    ((r - a.r) / dist) * d
                } else {
                    Vector3::zero()
                },
                r,
            }
        }
    }

    #[inline(always)]
    fn surface_area(&self) -> f32 {
        self.r * self.r
    }
}

impl BoundedBy<Sphere> for Triangle {
    fn bounds(&self) -> Sphere {
        let c = (self.a + self.b + self.c) / 3.0;
        let r = (self.a - c)
            .magnitude2()
            .max((self.b - c).magnitude2().max((self.c - c).magnitude2()))
            .sqrt();
        Sphere {
            c: Point3::from_vec(c),
            r,
        }
    }
}

impl BoundedBy<Sphere> for Rectangle {
    fn bounds(&self) -> Sphere {
        Sphere {
            c: self.c,
            r: (self.e[0] + self.e[1]).sqrt(),
        }
    }
}

/// Taking the bounding sphere of an AABB creates a sphere that contains all of
/// the extreme points of the AABB. Thus, a bounding sphere for an AABB created
/// from another bounding sphere will be larger in size than the original
/// sphere.
impl BoundedBy<Sphere> for AABB {
    fn bounds(&self) -> Sphere {
        Sphere {
            c: self.c,
            r: self.r.magnitude(),
        }
    }
}

impl BoundedBy<Sphere> for Capsule {
    fn bounds(&self) -> Sphere {
        // Include length of d with radius to cover all rotations
        let r = self.r + self.d.magnitude() * 0.5;
        Sphere {
            c: self.a + self.d*0.5,
            r,
        }
    }
}

#[cfg(test)]
mod tests {
    mod bounds {
        use cgmath::{Point3, Vector3};
        use bounds::Bound;
        use geom::{Sphere, AABB};
        use collision::{Contains, Overlaps};

        #[test]
        fn test_aabb() {
            let bound1 = AABB {
                c: Point3::new(0.0, 0.0, 0.0),
                r: Vector3::new(1.0, 1.0, 1.0),
            };
            let bound2 = AABB {
                c: Point3::new(0.0, 2.0, 0.0),
                r: Vector3::new(1.0, 1.0, 1.0),
            };
            let bound3 = AABB {
                c: Point3::new(0.0, 3.0, 0.0),
                r: Vector3::new(1.0, 1.0, 1.0),
            };
            let combined = AABB::combine(&bound1, &bound2);
            assert!(bound1.overlaps(&bound2));
            assert!(!bound1.overlaps(&bound3));
            assert!(!bound1.contains(&bound2));
            assert!(combined.contains(&bound1));
            assert!(combined.contains(&bound2));
            assert!(!combined.contains(&bound3));
        }

        #[test]
        fn test_sphere() {
            let bound1 = Sphere {
                c: Point3::new(0.0, 0.0, 0.0),
                r: 1.0,
            };
            let bound2 = Sphere {
                c: Point3::new(0.0, 2.0, 0.0),
                r: 1.0,
            };
            let bound3 = Sphere {
                c: Point3::new(0.0, 3.0, 0.0),
                r: 1.0,
            };
            let combined = Sphere::combine(&bound1, &bound2);
            assert!(bound1.overlaps(&bound2));
            assert!(!bound1.overlaps(&bound3));
            assert!(!bound1.contains(&bound2));
            assert!(combined.contains(&bound1));
            assert!(combined.contains(&bound2));
            assert!(!combined.contains(&bound3));
        }

        #[test]
        fn test_mixed() {
            use bounds::BoundedBy;

            let bound1 = Sphere {
                c: Point3::new(0.0, 0.0, 0.0),
                r: 1.0,
            };
            let bound2 = AABB {
                c: Point3::new(0.0, 2.0, 0.0),
                r: Vector3::new(1.0, 1.0, 1.0),
            };
            let bound3 = Sphere {
                c: Point3::new(0.0, 3.0, 0.0),
                r: 1.0,
            };
            let combined_sphere = Sphere::combine(&bound1, &bound2.bounds());
            let combined_aabb = AABB::combine(&bound1.bounds(), &bound2);
            let combined_unknown = Bound::combine(&bound1, &bound2.bounds());
            assert!(bound1.overlaps(&bound2));
            assert!(!bound1.overlaps(&bound3));
            let bound2_bound: Sphere = bound2.bounds();
            assert!(!bound1.contains(&bound2_bound));
            assert!(combined_sphere.contains(&bound1));
            assert!(combined_sphere.contains(&bound2_bound));
            assert!(!combined_sphere.contains(&bound3));
            let bound1_bound: AABB = bound1.bounds();
            assert!(combined_aabb.contains(&bound1_bound));
            assert!(combined_aabb.contains(&bound2));
            let bound3_bound: AABB = bound3.bounds();
            assert!(!combined_aabb.contains(&bound3_bound));
            assert!(combined_unknown.contains(&bound1));
            assert!(combined_unknown.contains(&bound2_bound));
            assert!(!combined_unknown.contains(&bound3));
        }
    }
}
