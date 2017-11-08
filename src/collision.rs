use std::f32;
use std::ops::Neg;

use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Quaternion,  Rotation,
             Vector3, Zero};

use bitset::FixedSizeBitSet;
use geom::*;

/// A type that can overlap another.
///
/// Overlaps is the most simple form of discrete collision detection and is
/// supported by every Bound on itself and almost every Shape.
/// 
/// No Moving object satisfies Overlaps.
pub trait Overlaps<RHS> {
    /// Returns true if the two objects overlap and false otherwise.
    fn overlaps(&self, rhs: &RHS) -> bool;
}

impl Overlaps<AABB> for AABB {
    fn overlaps(&self, rhs: &AABB) -> bool {
        // TODO: Add SIMD version of this
        (self.c.x - rhs.c.x).abs() <= (self.r.x + rhs.r.x)
            && (self.c.y - rhs.c.y).abs() <= (self.r.y + rhs.r.y)
            && (self.c.z - rhs.c.z).abs() <= (self.r.z + rhs.r.z)
    }
}

impl Overlaps<Sphere> for AABB {
    #[inline(always)]
    fn overlaps(&self, rhs: &Sphere) -> bool {
        rhs.overlaps(self)
    }
}

impl Overlaps<AABB> for Sphere {
    fn overlaps(&self, rhs: &AABB) -> bool {
        // TODO: Add SIMD version of this
        let mut d = 0.0;
        for i in 0..3 {
            let e = self.c[i] - (rhs.c[i] - rhs.r[i]);
            if e < 0.0 {
                if e < -self.r {
                    return false;
                }
                d += e * e;
            } else {
                // This feels... wrong given my AABB representation.
                let e = self.c[i] - (rhs.c[i] + rhs.r[i]);
                if e > 0.0 {
                    if e > self.r {
                        return false;
                    }
                    d += e * e;
                }
            }
        }
        d <= self.r * self.r
    }
}

impl Overlaps<Sphere> for Sphere {
    fn overlaps(&self, rhs: &Sphere) -> bool {
        let r = self.r + rhs.r;
        (rhs.c - self.c).magnitude2() <= r * r
    }
}

/// A type that can completely subsume another.
///
/// Contains is another common form of discrete collision detection. Most often
/// used internally to determine if an object contains a point.
pub trait Contains<RHS> {
    /// Returns true if the current object contains the argument.
    fn contains(&self, rhs: &RHS) -> bool;
}

impl Contains<Point3<f32>> for Plane {
    fn contains(&self, p: &Point3<f32>) -> bool {
        relative_eq!(self.n.dot(p.to_vec()), self.d, epsilon = COLLISION_EPSILON)
    }
}

impl Contains<Point3<f32>> for Triangle {
    fn contains(&self, p: &Point3<f32>) -> bool {
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
        u >= 0.0 && v >= 0.0 && (u + v) < 1.0
    }
}

impl Contains<Point3<f32>> for  Rectangle {
    fn contains(&self, p: &Point3<f32>) -> bool {
        // This might be sped up by simply finding the closest point and checking
        // if it is equal to p. I'm not sure though.
        let p = p.to_vec();
        let n = self.u[0].cross(self.u[1]);
        relative_eq!(p.dot(n), n.dot(self.c.to_vec()), epsilon = COLLISION_EPSILON)
            && p.dot(self.u[0]).abs() <= self.e[0]
            && p.dot(self.u[1]).abs() <= self.e[1]
    }
}

impl Contains<Point3<f32>> for AABB {
    fn contains(&self, p: &Point3<f32>) -> bool {
        (self.c.x - p.x).abs() <= self.r.x
            && (self.c.y - p.y).abs() <= self.r.y
            && (self.c.z - p.z).abs() <= self.r.z
    }
}

impl Contains<Point3<f32>> for Sphere {
    fn contains(&self, p: &Point3<f32>) -> bool {
         (p - self.c).magnitude2() <= self.r*self.r
    }
}


impl Contains<AABB> for AABB {
    fn contains(&self, rhs: &AABB) -> bool {
        let rhs_max = rhs.c + rhs.r;
        let rhs_min = rhs.c + -rhs.r;
        self.contains(&rhs_max) && self.contains(&rhs_min)
    }
}

/// An important consequence of Spheres being closed: two spheres with the same
/// center and radius contain each other. 
impl Contains<Sphere> for Sphere {
    fn contains(&self, rhs: &Sphere) -> bool {
        if self.r < rhs.r {
            return false;
        } 
        let r = self.r - rhs.r;
        (rhs.c - self.c).magnitude2() <= r * r
    }
}

/// A collision between a non-volumetric object and a volumetric object.
#[derive(Copy, Clone, Debug)]
pub struct Intersection {
    /// The point of intersection.
    pub p: Point3<f32>,
    /// The time of intersection relative to the velocity of the non-volumetric
    /// object.
    pub t: f32,
}

/// A type that can collide with a volumetric object and produce a single point
/// of contact.
///
/// Intersects for particle types implements collisions against volumetric and
/// planar geometries. Intersects is not a commutative operator. 
pub trait Intersects<RHS> : Particle {
    /// Returns an Intersection if one exists.
    fn intersection(&self, rhs: &RHS) -> Option<Intersection>;
}

impl<P: Particle> Intersects<Plane> for P {
    fn intersection(&self, p: &Plane) -> Option<Intersection> {
        let denom = p.n.dot(self.dir());
        if denom == 0.0 {
            return None;
        }
        let t = (p.d - p.n.dot(self.pos().to_vec())) / denom;
        if t <= 0.0 || t > P::DT {
            return None;
        }
        Intersection {
            p: self.pos() + self.dir() * t,
            t,
        }.into()
    }
}

impl<Part, Poly> Intersects<Poly> for Part
where
    Part: Particle,
    Poly: Polygon
{
    fn intersection(&self, poly: &Poly) -> Option<Intersection> {
        let p: Plane = (*poly).into();
        if let Some(inter) = self.intersection(&p) {
            if poly.contains(&inter.p) {
                return inter.into();
            }
        }
        None
    }
}

impl<P: Particle> Intersects<AABB> for P {
    fn intersection(&self, a: &AABB) -> Option<Intersection> {
        let (mut t_min, mut t_max): (f32, f32) = (0.0, f32::INFINITY);
        let p = self.pos();
        let d = self.dir();
        for dim in 0..3 {
            if d[dim].abs() < COLLISION_EPSILON {
                if (p[dim] - a.c[dim]).abs() > a.r[dim] {
                    return None;
                }
            } else {
                let ood = 1.0 / d[dim];
                let t1 = (a.c[dim] - a.r[dim] - p[dim]) * ood;
                let t2 = (a.c[dim] + a.r[dim] - p[dim]) * ood;
                if t1 > t2 {
                    t_min = t_min.max(t2);
                    t_max = t_max.min(t1);
                } else {
                    t_min = t_min.max(t1);
                    t_max = t_max.min(t2);
                }
                if t_min > t_max {
                    return None;
                }
            }
        }
        if t_min > P::DT {
            return None;
        }
        Intersection {
            p: p + d * t_min,
            t: t_min,
        }.into()
    }
}

impl<P: Particle> Intersects<Sphere> for P {
    fn intersection(&self, s: &Sphere) -> Option<Intersection> {
        let p = self.pos();
        let d = self.dir();
        let m = p - s.c;
        let a = d.magnitude2();
        let b = m.dot(d);
        let c = m.magnitude2() - s.r * s.r;
        if c > 0.0 && b > 0.0 {
            return None;
        }
        let discr = b * b - a * c;
        if discr < 0.0 {
            return None;
        }
        let t = ((-b - discr.sqrt()) / a).max(0.0);
        if t > P::DT {
            return None;
        }
        Intersection {
            p: p + t * d,
            t,
        }.into()
    }
}

impl<P: Particle> Intersects<Capsule> for P {
    fn intersection(&self, cap: &Capsule) -> Option<Intersection> {
        // This code is horrible and needs to be completely redone.
        let p = self.pos();
        let d = self.dir();
        let m = p - cap.a;
        let md = m.dot(cap.d);
        let nd = d.dot(cap.d);
        let dd = cap.d.dot(cap.d);
        let nn = d.magnitude2();
        let mn = m.dot(d);
        let a = dd * nn - nd * nd;
        let k = m.magnitude2() - cap.r * cap.r;
        if a.abs() < COLLISION_EPSILON {
            let (b, c) = if md < 0.0 {
                (mn, k)
            } else if md > dd {
                let m2 = p - (cap.a + cap.d);
                (m2.dot(d), m2.magnitude2() - cap.r * cap.r)
            } else {
                // Already colliding
                return None;
            };
            if c > 0.0 && b > 0.0 {
                 return None;
             }
            let discr = b * b - nn*c;
            if discr < 0.0 {
                return None;
            }
            let t = ((-b - discr.sqrt()) / nn).max(0.0);
            if t > P::DT {
                return None;
            }
            return Intersection {
                p: p + t * d,
                t,
            }.into()
        }
        let c = dd * k - md * md;
        let b = dd * mn - nd * md;
        let discr = b * b - a * c;
        if discr < 0.0 {
            return None;
        }
        let t = (-b - discr.sqrt()) / a;
        if t < 0.0 {
            // Intersection is behind ray.
            return None;
        }
        // I would like to avoid taking a square root twice, but for now that's
        // the price we pay.
        let t = if md + t * nd < 0.0 {
            if mn > 0.0 && k > 0.0 {
                return None;
            }
            let discr = mn * mn - nn*k;
            if discr < 0.0 {
                return None;
            }
            ((-mn - discr.sqrt()) / nn).max(0.0)
        } else if md + t * nd > dd {
            let m2 = p - (cap.a + cap.d);
            let b = m2.dot(d);
            let c = m2.magnitude2() - cap.r * cap.r;
            if c > 0.0 && b > 0.0 {
                return None;
            }
            let discr = b * b - nn*c;
            if discr < 0.0 {
                return None;
            }
            ((-b - discr.sqrt()) / nn).max(0.0)
        } else {
            t
        };
        if t > P::DT {
            return None;
        }
        Intersection {
                p: p + t * d,
                t,
        }.into()
    }
}

impl<P: Particle> Intersects<Moving<Sphere>> for P {
    fn intersection(&self, s: &Moving<Sphere>) -> Option<Intersection> {
        // ray intersection with moving spheres is exactly the same problem as
        // capsule intersection
        self.intersection(
            &Capsule {
                a: s.0.c,
                d: s.1,
                r: s.0.r,
            }
        )
    }
}


/// A point of contact between two objects occuring during a timestep.
///
/// Contact models a contact point generated by at least one moving volumetric
/// geometry and some other geometry.
#[derive(Copy, Clone, Debug)]
pub struct Contact {
    /// Contact point at time of collision for collider in global coordinates.
    pub a: Point3<f32>,
    /// Contact point at time of collision for collidee in global coordinates.
    pub b: Point3<f32>,
    /// Collision normal on the surface of the collider.
    pub n: Vector3<f32>,
    /// Time of impact. This is guaranteed to be in the interal [0, 1]. Contacts
    /// with a time of 0 can be considered resting contacts or a contact from
    /// the previous frame.
    pub t: f32,
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

/// A type that can produce a point of contact with another.
///
/// Contacts models a hybrid discrete/continuous approach to collision detection
/// that uses continuous detection to find precise contact points and collision
/// normals over discrete timesteps.
///
/// Either the type implementing Contacts or the type paramater must be
/// Volumetric. Both can be as well. In addition, at least one of the types must
/// be Moving in some way.
///
/// A Contact collision can potentially produce multiple contacts, in the case of
/// Capsule/Polygon collision and Compound collisions. This is handled by passing
/// a closure to the collision handler.
pub trait Contacts<RHS> {
    /// Calls the closure for each contact found. Returns true if any contact was
    /// found.
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, callback: F) -> bool;

    /// Returns the last contact found, if one exists.
    fn last_contact(&self, rhs: &RHS) -> Option<Contact> {
        let mut contact = None;
        self.contacts(rhs, |c|{ contact = Some(c); });
        contact
    }
}

macro_rules! commute_contacts {
    (
        $recv:ty, $arg:ty
    ) => {
        impl Contacts<$arg> for $recv {
            fn contacts<F: FnMut(Contact)>(&self, rhs: &$arg, mut callback: F) -> bool {
                rhs.contacts(self, |c: Contact|{ callback(-c) })
            }
        }
    };
}

/*
impl Contacts<> for  {
    fn contacts<F: FnMut(Contact)>(&self, : &, mut callback: F) -> bool {

    }
}
*/

impl Contacts<Moving<Sphere>> for Plane {
    fn contacts<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
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

impl Contacts<Moving<Capsule>> for Plane {
    fn contacts<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
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
        self.contacts(&ms, callback)
    }
}

commute_contacts!{ Moving<Sphere>, Plane }
commute_contacts!{ Moving<Capsule>, Plane }

impl<Poly: Polygon> Contacts<Moving<Sphere>> for Poly {
    fn contacts<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
        let &Moving(s, v) = sphere;
        let mut collision = false;
        let p: Plane = (*self).into();
        p.contacts(sphere, |contact: Contact| {
            // If the point lies on the face we know automatically that it is
            // a valid collision.
            if self.contains(&contact.a) {
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
            for edge_i in 0..Poly::NUM_VERTICES{
                let (a, b) = self.edge(edge_i);
                let v1 = self.vertex(a);
                let v2 = self.vertex(b);
                let c = Capsule{ a: v1, d: v2 - v1, r: s.r };
                if let Some(i) = ray.intersection(&c) {
                    if i.t <= 1.0 && i.t < first_t {
                        first_t = i.t;
                        tri_p = Segment::from((v1, v2)).min_dist(&i.p);
                    }
                }
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

commute_contacts!{ Moving<Sphere>, Triangle }
commute_contacts!{ Moving<Sphere>, Rectangle }
commute_contacts!{ Moving<Capsule>, Triangle }
commute_contacts!{ Moving<Capsule>, Rectangle }

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
impl<Poly: Polygon> Contacts<Moving<Capsule>> for Poly  {
    fn contacts<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
        let &Moving(c, v) = capsule;
        let p: Plane = (*self).into();
        // Check if the capsule is already colliding.
        let denom = p.n.dot(c.d.normalize());
        if denom.abs() > COLLISION_EPSILON {
            let t = (p.d - p.n.dot(c.a.to_vec())) / denom;
            if t <= 1.0 && t >= 0.0 {
                // Already colliding the plane
                let q = c.a + c.d * t;
                if self.contains(&q) {
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
            p.last_contact(&start_sphere)
            .map_or_else(
                || { 
                    p.last_contact(&end_sphere)
                        .map(|c1: Contact| { (c1, -c.d) })
                },
                |c1: Contact| {
                    p.last_contact(&end_sphere)
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
            if self.contains(&contact.a) {
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
                for edge_i in 0..Poly::NUM_VERTICES {
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
                for edge_i in 0..Poly::NUM_VERTICES {
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
        }
        // Intersect the Minkowski sum of the triangle and the capsule with
        // the ray originating at the capsule's origin.
        // Unfortunately, we have to use an vector here because rust does not
        // support stack allocation if we want to support any number of
        // edges. In order to subvert that as much as possilbe we're using
        // a pool Block64 as a bitset.
        // Thus, any face that wishes to compute here may only have 64
        // vertices (for now). Hopefully this will change when generic constants
        // are added.
        // Honestly, that is a rather reasonable requirement. I'm certainly
        // never going to use more than four.
        if Poly::NUM_VERTICES > 64 {
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
        for edge_i in 0..Poly::NUM_VERTICES {
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
            if let Some(inter) = ray.intersection(&edge_sum) {
                if inter.t > best_par.0.min(1.0) {
                    continue;
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

            } else if let Some(inter) = ray.intersection(&Capsule{ a: edge_a, d: -c.d, r: c.r }) {
                if inter.t > best_par.0.min(1.0) {
                    continue;
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
            }
        }
        // Perform edge collisions.
        let mut best_sum = (
            f32::INFINITY,
            Point3::new(0.0,0.0,0.0),
        );

        for edge_i in 0..Poly::NUM_VERTICES {
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
            p.contacts(&Moving::sweep(s, v), |contact: Contact| {
                if best_sum.0 > contact.t
                    && (tris[0].contains(&contact.a)
                        || tris[1].contains(&contact.b))
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
                    if let Some(inter) = ray.intersection(&bottom_edge) {
                        if inter.t <= 1.0 && inter.t <= best_sum.0 {
                            let q: Point3<f32> = Segment::from((edge_a, edge_b))
                                .min_dist(&inter.p);
                            best_sum = (inter.t, q);
                        }
                    }
                    let top_edge = Capsule {
                        a: edge_a + -c.d,
                        d: edge_b - edge_a,
                        r: c.r
                    };
                    if let Some(inter) = ray.intersection(&top_edge) {
                        if inter.t <= 1.0 && inter.t <= best_sum.0 {
                            let plane_p = inter.p + c.d;
                            let q: Point3<f32> = Segment::from((edge_a, edge_b))
                                .min_dist(&plane_p);
                            best_sum = (inter.t, q);
                        }
                    }
                    // Hopefully this gets unrolled.
                    for &(vert, is_parallel) in [
                        (edge_a, a_on_parallel_edge),
                        (edge_b, b_on_parallel_edge),
                    ].iter() {
                        if is_parallel {
                            continue;
                        }
                        let cap = Capsule{ a: vert, d: -c.d, r: c.r };
                        if let Some(inter) = ray.intersection(&cap) {
                            if inter.t <= 1.0 && inter.t <= best_sum.0 {
                                best_sum = (inter.t, vert);
                            }
                        }
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


impl Contacts<Moving<Sphere>> for Sphere {
    fn contacts<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
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
        if let Some(intersect) = ray.intersection(&Sphere{ c: s.c, r }) {
            if intersect.t <= 1.0 {
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
                });
                return true;
            }
        }
        false
    }
}

commute_contacts!{ Sphere, Moving<Capsule> }

impl Contacts<Moving<Sphere>> for Capsule {
    fn contacts<F: FnMut(Contact)>(&self, sphere: &Moving<Sphere>, mut callback: F) -> bool {
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
        if let Some(intersect) = ray.intersection(&Capsule{ r: s.r + self.r, ..*self }) {
            if intersect.t <= 1.0 {
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
                });
                return true;
            }
        }
        false
    }
}

impl Contacts<Moving<Capsule>> for Capsule {
    fn contacts<F: FnMut(Contact)>(&self, capsule: &Moving<Capsule>, mut callback: F) -> bool {
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
            return ss.contacts(capsule, callback);
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
                return self.contacts(
                    &Moving::sweep(Sphere{ c: c_a + c_d, r: c.r }, v), callback
                );
            }
            if t_min >= 1.0 {
                return self.contacts(
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
            return self.contacts(
                &Moving::sweep(Sphere{ c: c_a + c_d, r: c.r }, v), callback
            );
        }
        if t_min >= 1.0 {
            return self.contacts(
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

impl<Recv, Arg> Contacts<Arg> for Moving<Recv>
where
    Arg: Shape + Copy,
    Recv: Contacts<Moving<Arg>> + Shape + Copy
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &Arg, mut callback: F) -> bool {
        let rhs_moving = Moving::sweep(*rhs, -self.1);
        self.0.contacts(&rhs_moving, |c| {
            let d = self.1 * c.t;
            let a = c.a + d;
            let b = c.b + d;
            callback(Contact{ a, b, ..c });
        })
    }
}

/// Any two moving object collision can be reduced to a one moving one static
/// object collision. This is done by finding the relative velocity between the
/// two objects.
impl<Recv, Arg> Contacts<Moving<Arg>> for Moving<Recv>
where
    Recv: Contacts<Moving<Arg>> + Shape + Copy,
    Arg: Shape + Copy
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &Moving<Arg>, mut callback: F) -> bool {
        let Moving(geom_a, v_a) = *self;
        let Moving(geom_b, v_b) = *rhs;
        geom_a.contacts(&Moving::sweep(geom_b, v_b - v_a), |c| {
            let a = c.a + v_a * c.t;
            let b = c.b + v_a * c.t;
            callback(Contact{ a, b, ..c })
        })
    }
}

/// A point of contact between two objects that includes the contact points for
/// each object in terms of the object's center.
///
/// A LocalContact is derived from a regular Contact and includes all of the same
/// information including the contact points in the local coordinates of the
/// objects they belong to.
#[derive(Copy, Clone, Debug)]
pub struct LocalContact {
    /// Global contact point `a` relative to the center of object a at the time of
    /// collision.
    pub local_a: Point3<f32>,
    /// Global contact point `b` relative to the center of object b at the time of
    /// collision.
    pub local_b: Point3<f32>,
    /// Contact the LocalContact was derived from.
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

/// A type that can produce a point of contact with another and convert it to
/// local coordinates.
///
/// Because the points in a Contact are global at the time of collision, they
/// are not immediately useful to physics resolution. Upon being translated to
/// local coordinates, which refer to the points as if the origin was the center
/// of the object they belong to.
pub trait LocalContacts<RHS> {
    /// Calls the closure for each contact found. Returns true if any contact was
    /// found.
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &RHS, callback: F) -> bool;
    
    /// Returns the last contact found, if one exists.
    fn last_local_contact(&self, rhs: &RHS) -> Option<LocalContact> {
        let mut contact = None;
        self.local_contacts(rhs, |c|{ contact = Some(c); });
        contact
    }
}

impl<Recv, Arg> LocalContacts<Arg> for Recv
where
    Recv: Contacts<Arg> + Shape + Delta,
    Arg: Shape + Delta
{
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &Arg, mut callback: F) -> bool {
        self.contacts(rhs, |c| {
            let a_c = self.center() + self.delta() * c.t;
            let b_c = rhs.center() + rhs.delta() * c.t;
            callback(LocalContact {
                local_a: c.a + -a_c.to_vec(),
                local_b: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}

impl<Recv, Arg> LocalContacts<Moving<Arg>> for Recv
where
    Recv: Contacts<Moving<Arg>> + Shape,
    Arg: Shape + Copy 
{
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &Moving<Arg>, mut callback: F) -> bool {
        self.contacts(rhs, |c| {
            let a_c = self.center();
            let b_c = rhs.as_ref().center() + rhs.vel() * c.t;
            callback(LocalContact {
                local_a: c.a + -a_c.to_vec(),
                local_b: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}

impl<Recv, Arg> LocalContacts<Arg> for Moving<Recv>
where
    Recv: Shape + Copy,
    Arg: Contacts<Moving<Recv>> + Shape + Copy
{
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &Arg, mut callback: F) -> bool {
        rhs.contacts(self, |c| {
            let a_c = self.as_ref().center() + self.vel() * c.t;
            let b_c = rhs.center();
            callback(LocalContact {
                local_a: c.b + -a_c.to_vec(),
                local_b: c.a + -b_c.to_vec(),
                global: -c
            })
        })
    }
}

impl<Recv, Arg> LocalContacts<Moving<Arg>> for Moving<Recv>
where
    Recv: Contacts<Moving<Arg>> + Shape + Copy,
    Arg: Shape + Copy
{
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &Moving<Arg>, mut callback: F) -> bool {
        // Need to take a further look at this
        let Moving(geom_a, v_a) = *self;
        let Moving(geom_b, v_b) = *rhs;
        geom_a.contacts(&Moving::sweep(geom_b, v_b - v_a), |c| {
            let a = c.a + v_a * c.t;
            let b = c.b + v_a * c.t;
            let local_a = a + -geom_a.center().to_vec();
            let local_b = b + -(rhs.as_ref().center() + v_b * c.t).to_vec();
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
    mod intersections {
        use cgmath::{InnerSpace, Point3, Vector3};
        use geom;
        use geom::{Capsule, Ray};
        use collision::{Intersects};

        #[test]
        fn test_ray_intersections() {
            let c = Capsule {
                a: Point3::new(0.0, 0.0, 0.0),
                d: Vector3::new(1.0, 0.0, 0.0),
                r: 1.0,
            };
            let r = Ray {
                p: Point3::new(1.0, -3.0, 0.0),
                d: Vector3::new(-0.25, 1.0, 0.0).normalize(),
            };
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
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
            let intersection = r.intersection(&c).unwrap();
            assert_relative_eq!(intersection.p, Point3::new(1.8660254037844386, 0.5, 0.0));
            assert_relative_eq!(intersection.t, 1.13397459621556196,
                                epsilon = geom::COLLISION_EPSILON);
        }
    }

    mod spheres {
        use cgmath::{Point3, Vector3};
        use geom;
        use geom::{Moving, Sphere, Rect, Triangle};
        use collision::*;

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
            
            let collision: Contact = s1.last_contact(&s2).unwrap();
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
            assert!(floor.contacts(&sphere_collide_center, |c: Contact| {
                assert_eq!(c.a, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.b, Point3::new(0.0, 1.0, 0.0));
                assert_eq!(c.t, 1.0);
                assert_eq!(c.n, Vector3::new(0.0, 1.0, 0.0));
            }));
            assert!(sphere_collide_center.contacts(&floor, |c: Contact| {
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
            assert!(floor.contacts(&sphere_collide_center_2s, |c: Contact| {
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
            assert!(floor.contacts(&sphere_collide_corner, |c: Contact| {
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
            assert!(!floor.contacts(&sphere_miss_corner, |_:Contact| {}));
        }

        #[test]
        fn test_tri_collision() {
            let floor = Triangle {
                a: Vector3::new(1.0, 1.0, 0.0),
                c: Vector3::new(0.0, 1.0, 1.0),
                b: Vector3::new(0.0, 1.0, -1.0),
            };
            let sphere_collide_center = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.0, -10.0, 0.0),
            );
            assert!(floor.contacts(&sphere_collide_center, |c: Contact| {
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
            assert!(floor.contacts(&sphere_collide_corner, |c: Contact| {
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
            assert!(!floor.contacts(&sphere_miss_corner, |c: Contact| { panic!("C = {:?}", c); }));
            let sphere_collide_edge = Moving::sweep(
                Sphere {
                    c: Point3::new(0.0, 13.0, 0.0),
                    r: 2.0,
                },
                Vector3::new(0.5, -10.0, 0.5),
            );
            assert!(floor.contacts(&sphere_collide_edge, |c: Contact| {
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
        use geom::{Capsule, Sphere, Moving, Triangle, Rect};
        use collision::{Contacts, Contact};

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
            let collision: Contact = c.last_contact(&s).unwrap();
            assert_eq!(collision.t, 1.0);
            assert_eq!(collision.a, Point3::new(2.0, 3.0, 5.5));
            assert_eq!(collision.b, Point3::new(2.0, 3.0, 5.5));
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            let collision: Contact = s.last_contact(&c).unwrap();
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
            floor.contacts(&capsule_level_off_center, |c| contacts.push(c));
            assert_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(1.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_relative_eq!(contacts[1].a, Point3::new(3.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            contacts.clear();
        }
        
        #[test]
        fn test_tri_collision() {
            let floor = Triangle {
                a: Vector3::new(1.0, 1.0, 0.0),
                c: Vector3::new(0.0, 1.0, 1.0),
                b: Vector3::new(0.0, 1.0, -1.0),
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
            floor.contacts(&capsule_clip_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_clip_off_center, |c| contacts.push(c));
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
            floor.contacts(&capsule_clip_off_center, |c| contacts.push(c));
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
            floor.contacts(&capsule_through_center, |c| contacts.push(c));
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

            let collision: Contact = floor.last_contact(&capsule_tilted_center).unwrap();
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
            let collision: Contact = floor.last_contact(&capsule_tilted_center).unwrap();
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_parallel_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_perp_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_perp_to_edge, |c| contacts.push(c));
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
            floor.contacts(&capsule_next_to_vert, |c| contacts.push(c));
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
            floor.contacts(&capsule_next_to_vert_skewed, |c| contacts.push(c));
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
            floor.contacts(&capsule_intersects_tri_plane, |c| contacts.push(c));
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
            floor.contacts(&capsule_intersects_tri_plane, |c| contacts.push(c));
            assert_relative_eq!(contacts[0].t, 1.0);
            assert_relative_eq!(contacts[0].a, Point3::new(0.0, 1.0, 0.0), epsilon = geom::COLLISION_EPSILON);
            assert_eq!(contacts.len(), 1);
            contacts.clear();
        }
    }        
}
        
