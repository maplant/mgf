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
use std::ops::{AddAssign, SubAssign};
use cgmath::{EuclideanSpace, InnerSpace, SquareMatrix, Rotation,
             Matrix, Matrix3, Point3, Quaternion, Vector3, One, Zero};

use bounds::*;
use compound::*;
use collision::*;
use geom::*;
use manifold::*;

/// An type that has a moment of inertia.
pub trait Inertia {
    fn tensor(&self, m: f32) -> Matrix3<f32>;
}

impl Inertia for Sphere {
    fn tensor(&self, m: f32) -> Matrix3<f32> {
        let i = 0.4 * m * self.r * self.r;
        Matrix3::new(i, 0.0, 0.0, 0.0, i, 0.0, 0.0, 0.0, i)
    }
}

impl Inertia for Capsule {
    fn tensor(&self, m: f32) -> Matrix3<f32> {
        // Compute the inertia tensor for a capsule aligned vertically along the
        // y axis
        let h = self.d.magnitude();
        let r = self.r;
        // Distribute the mass between the hemispheres and cylinder by volume.
        let mh = m * 2.0 * r / (4.0 * r + 3.0 * h);
        let mc = m * h / (4.0 / 3.0 * r + h);
        // Tensor for cylinder
        let ic_x = 1.0 / 12.0 * mc * (3.0 * r * r + h * h);
        let ic_y = 0.5 * mc * r * r;
        let ic_z = ic_x;
        // Tensor for two hemispheres
        let is_x = mh * (3.0 * r + 2.0 * h) / 4.0 * h;
        let is_y = 4.0 / 5.0 * mh * r * r;
        let is_z = is_x;
        // Resulting tensor:
        let (i_x, i_y, i_z) = (ic_x + is_x, ic_y + is_y, ic_z + is_z);
        // Find the rotation of the capsule
        let dst = self.d;
        let src = Vector3::new(0.0, 1.0, 0.0) * h;
        let rot = Matrix3::<f32>::from(Quaternion::from_arc(src, dst, None));
        rot * Matrix3::new(i_x, 0.0, 0.0, 0.0, i_y, 0.0, 0.0, 0.0, i_z) * rot.transpose()
    }
}

impl Inertia for Component {
    fn tensor(&self, m: f32) -> Matrix3<f32> {
        match self {
            &Component::Sphere(s) => s.tensor(m),
            &Component::Capsule(c) => c.tensor(m),
        }
    }
}

/// A generic physical body that has a mass, a volume, and experiences linear and
/// rotational movement.
///
/// A RigidBody is technically a shape, although its geomtries are all considered
/// to be in motion.
#[derive(Clone)]
pub struct RigidBody {
    pub restitution: f32,
    pub friction: f32,
    pub total_inv_mass: f32,
    pub inv_moment_body: Matrix3<f32>,
    pub inv_moment: Matrix3<f32>,
    pub v: Vector3<f32>,
    pub omega: Vector3<f32>,
    pub v_step: Vector3<f32>,
    pub linear_m: Vector3<f32>,
    pub angular_m: Vector3<f32>,
    pub force: Vector3<f32>,
    pub torque: Vector3<f32>,
    // Position and rotation are stored in collider.
    pub collider: Compound,
    pub component_masses: Vec<f32>,
}

impl RigidBody {
    /// Construct a new RigidBody from a Vec of Components and masses.
    pub fn new(
        restitution: f32,
        friction: f32,
        world_force: Vector3<f32>,
        mut shapes: Vec<Component>,
        masses: Vec<f32>
    ) -> Self {
        // Calculate the center of mass of the RigidBody.
        let mut total_mass = 0.0;
        let mut center = Point3::new(0.0, 0.0, 0.0);
        for (shape, mass) in shapes.iter().zip(masses.iter()) {
            center += shape.center().to_vec() * *mass;
            total_mass += *mass;
        }
        // Displace the shapes by the center of mass.
        let inv_mass = 1.0 / total_mass;
        center *= inv_mass;
        for shape in shapes.iter_mut() {
            *shape -= center.to_vec();
        }
        // Now we can properly calculate the tensor of the body
        let tensor = shapes.iter().zip(masses.iter())
            .fold(Matrix3::zero(), |tensor_sum, (shape, mass)| {
                let disp = shape.center().to_vec();
                let outer_prod = Matrix3::from_cols(
                    disp * disp.x,
                    disp * disp.y,
                    disp * disp.z
                );
                tensor_sum + shape.tensor(*mass)
                    + *mass * (Matrix3::one() * disp.magnitude2() - outer_prod)
            });
        let inv_moment = tensor.invert().unwrap();
        RigidBody {
            restitution,
            friction,
            total_inv_mass: inv_mass,
            inv_moment_body: inv_moment,
            inv_moment,
            v: Vector3::zero(),
            omega: Vector3::zero(),
            v_step: Vector3::zero(),
            linear_m: Vector3::zero(),
            angular_m: Vector3::zero(),
            force: world_force / inv_mass,
            torque: Vector3::zero(),
            collider: Compound::new(shapes),
            component_masses: masses,
        }
    }
}    

impl AddAssign<Vector3<f32>> for RigidBody {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.collider.disp += v
    }
}

impl SubAssign<Vector3<f32>> for RigidBody {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.collider.disp -= v
    }
}

impl Shape for RigidBody {
    /// The center returned by a RigidBody is the center of mass for that
    /// object. 
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.collider.disp)
    }
}

impl Delta for RigidBody {
    fn delta(&self) -> Vector3<f32> {
        self.v_step
    }
}

impl BoundedBy<AABB> for RigidBody {
    fn bounds(&self) -> AABB {
        let b1 = self.collider.bounds();
        let b2 = b1 + self.v_step;
        Bound::combine(&b1, &b2)
    }
}

impl BoundedBy<Sphere> for RigidBody {
    fn bounds(&self) -> Sphere {
        let b1 = self.collider.bounds();
        let b2 = b1 + self.v_step;
        Bound::combine(&b1, &b2)
    }
}

impl<RHS> Contacts<RHS> for RigidBody
where
    RHS: Contacts<Moving<Component>> + BoundedBy<AABB>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        let conj_rot = self.collider.rot.conjugate();
        let mut rhs_bounds = rhs.bounds().rotate(conj_rot);
        let rhs_center = rhs_bounds.center();
        let bounds_disp = conj_rot.rotate_point(rhs_center + -self.collider.disp) + self.collider.disp;
        rhs_bounds.set_pos(bounds_disp);
        let rhs_bounds: AABB = Moving::sweep(rhs_bounds, -self.v_step).bounds();
        let mut collided = false;
        self.collider.bvh.query(&rhs_bounds, |&comp| {
            let shape = Moving::sweep(
                comp.rotate(self.collider.rot) + self.collider.disp,
                self.v_step
            );
            rhs.contacts(&shape, |c| { collided = true; callback(-c) });
        });
        collided
    }
}

/// A physical object that cannot move.
///
/// A static body is an immovable object that can be constructed from any
/// shape, the only other parameter required is a value for friction.
/// 
/// A StaticBody is not a Shape and cannot be moved. StaticBodies contain
/// a reference to the Shape they mimic and thus should not live very long.
pub struct StaticBody<'a, S: Shape + 'a> {
    pub friction: f32,
    pub shape: &'a S,
}

impl<'a, S: Shape + 'a> StaticBody<'a, S> {
    /// Construct a new StaticBody from a Shape.
    pub fn new(friction: f32, shape: &'a S) -> Self {
        StaticBody {
            friction,
            shape
        }
    }
}

impl<'a, S, B> BoundedBy<B> for StaticBody<'a, S>
where
    B: Bound,
    S: Shape + BoundedBy<B>
{
    fn bounds(&self) -> B {
        self.shape.bounds()
    }
}

impl<'a, P, S> Intersects<StaticBody<'a, S>> for P
where
    P: Intersects<S>,
    S: Shape
{
    fn intersection(&self, rhs: &StaticBody<'a, S>) -> Option<Intersection> {
        self.intersection(rhs.shape)
    }
}

impl<'a, S, RHS> Contacts<RHS> for StaticBody<'a, S>
where
    S: Contacts<RHS> + Shape
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, callback: F) -> bool {
        self.shape.contacts(rhs, callback)
    }
}

impl<'a, S, Recv> LocalContacts<StaticBody<'a, S>> for Recv
where
    S: Shape,
    Recv: Contacts<S> + Shape + Delta
{
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &StaticBody<'a, S>, mut callback: F) -> bool {
        self.contacts(rhs.shape, |c| {
            let a_c = self.center() + self.delta() * c.t;
            let b_c = rhs.shape.center();
            callback(LocalContact {
                local_b: c.a + -a_c.to_vec(),
                local_a: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}

/// A description of the physical state of an object.
#[derive(Copy, Clone)]
pub struct PhysicsState {
    /// Restitution is a measure of how much kinetic energy is retained in a
    /// collision. 100% of kinetic energy retention corresponds to a coefficient
    /// of one.
    pub restitution: f32,
    /// We simplify friction in this case to simply be a ratio of the normal force
    /// applied tangentially to an object during collision.
    pub friction: f32,
    /// We only ever need inverse mass for calcuations, plus it gives a neat
    /// advantage that we can represent immovable objects with an infinite mass,
    /// or an inverse mass of zero.
    pub inv_mass: f32,
    /// The inverse moment of inertia tensor of the object correctly oriented.
    pub inv_moment: Matrix3<f32>,
    /// The position of the object at the end of the timestep.
    pub x: Point3<f32>,
    /// The linear velocity of the object
    pub v: Vector3<f32>,
    /// The angular velocity of the object
    pub omega: Vector3<f32>,
}

/// A type that exhibits physical properties.
///
/// A PhysicsObject's primary function is to return a PhysicsState to be used
/// during collision resolution. Beyond that it has various methods to be updated
/// or choose to ignore such updates (for example, calling `apply_impulse` on a
/// StaticBody is a no-op).
pub trait PhysicsObject {
    /// Integrate the object over the timestep
    fn integrate(&mut self, dt: f32);

    /// Return the physics state of the object
    fn state(&self) -> PhysicsState;

    /// Apply linear and angular impulse to the object
    fn apply_impulse(&mut self, linear: Vector3<f32>, angular: Vector3<f32>);

    /// Update the linear and angular positional derivatives
    fn update_dx(&mut self);
}

impl PhysicsObject for RigidBody {
    fn integrate(&mut self, dt: f32) {
        self.collider += self.v_step;
        self.collider.rot = (self.collider.rot
                             + Quaternion::from_sv(0.0, self.omega)
                             * 0.5 * self.collider.rot).normalize();
        let r = Matrix3::from(self.collider.rot);
        self.inv_moment = r * self.inv_moment_body * r.transpose();
        self.linear_m += self.force * dt;
        self.angular_m += self.torque * dt;
        self.update_dx();
        self.v_step = self.v * dt;
    }

    fn state(&self) -> PhysicsState {
        PhysicsState {
            restitution: self.restitution,
            friction: self.friction,
            inv_mass: self.total_inv_mass,
            inv_moment: self.inv_moment,
            x: self.collider.center() + self.v_step,
            v: self.v,
            omega: self.omega,
        }
    }

    fn apply_impulse(&mut self, linear: Vector3<f32>, angular: Vector3<f32>) {
        self.linear_m += linear;
        self.angular_m += angular;
    }

    fn update_dx(&mut self) {
        self.v = self.linear_m * self.total_inv_mass;
        self.omega = self.inv_moment * self.angular_m;
    }
}

impl<'a, S: Shape + 'a> PhysicsObject for StaticBody<'a, S> {
    #[inline(always)]
    fn integrate(&mut self, _dt: f32) {
        // Do nothing
    }

    fn state(&self) -> PhysicsState {
        PhysicsState {
            restitution: 0.0,
            friction: self.friction,
            inv_mass: 0.0,  // Infinite mass
            inv_moment: Matrix3::zero(),
            x: self.shape.center(),
            v: Vector3::zero(),
            omega: Vector3::zero(),
        }
    }
    
    #[inline(always)]
    fn apply_impulse(&mut self, _linear: Vector3<f32>, _angular: Vector3<f32>) {
        // Do nothing
    }

    #[inline(always)]
    fn update_dx(&mut self) {
        // Do nothing
    }
}

/// A set of constants necessary when performing collision resolution.
pub trait PhysicsConfig {
    const PENETRATION_SLOP: f32;
    const BAUMGARTE: f32; 
}

/// The suggested set of parameters to use when resolving collisions.
pub struct DefaultPhysConfig {}

impl PhysicsConfig for DefaultPhysConfig {
    const PENETRATION_SLOP: f32 = 0.005;
    const BAUMGARTE: f32 = 0.2;
}

impl PhysicsState {
    /// Perform collision resolution on two physics states from a contact and a
    /// timestep.
    pub fn resolve_contact<Config, ObjA, ObjB>(
        obj_a: &mut ObjA,
        obj_b: &mut ObjB,
        contact: &LocalContact,
        dt: f32
    ) where Config: PhysicsConfig,
            ObjA: PhysicsObject,
            ObjB: PhysicsObject
    {
        let (state_a, state_b) = (obj_a.state(), obj_b.state());
        let (xa, va, oa) = (state_a.x, state_a.v, state_a.omega);
        let (xb, vb, ob) = (state_b.x, state_b.v, state_b.omega);
        let ca = contact.local_a.to_vec() + xa.to_vec();
        let cb = contact.local_b.to_vec() + xb.to_vec();
        let ra = ca - xa.to_vec();
        let rb = cb - xb.to_vec();
        let ra_cn = ra.cross(contact.global.n);
        let rb_cn = rb.cross(contact.global.n);

        // Penetration is defined as the distance between the two contact points
        // dotted with the normal vector.
        let pen = (cb - ca).dot(contact.global.n);

        // The value of restitutiuon we want to use is the maximum of the two
        // object's restitution.
        let max_restitution = state_a.restitution.max(state_b.restitution);

        let dv = vb + ob.cross(rb) - va - oa.cross(ra);
        let rel_v = dv.dot(contact.global.n);

        let bias = Config::BAUMGARTE / dt * if pen > 0.0 {
            0.0
        } else {
            -pen + Config::PENETRATION_SLOP
        } + if rel_v < -1.0 {
            -max_restitution * rel_v
        } else {
            0.0
        };

        let normal_mass = 1.0 /
            (state_a.inv_mass + ra_cn.dot(state_a.inv_moment * ra_cn)
             + state_b.inv_mass + rb_cn.dot(state_b.inv_moment * rb_cn));

        // Calculate and clamp impulse due to collision resolution
        let lambda = (normal_mass * (-rel_v + bias)).max(0.0);
        let impulse = contact.global.n * lambda;

        // Apply impulse
        obj_a.apply_impulse(-impulse, -state_a.inv_moment * ra.cross(impulse));
        obj_b.apply_impulse( impulse,  state_b.inv_moment * rb.cross(impulse));

        // Reduce the normal impulse by some constant and calculate friction.
        let normal_impulse = lambda / 5.0;
        let mix_friction = (state_a.friction * state_b.friction).sqrt();
        let tangents = contact.global.compute_basis();
        let mut impulse = Vector3::zero();
        for tangent in tangents.iter() {
            let ra_ct = ra.cross(*tangent);
            let rb_ct = rb.cross(*tangent);
            let tangent_mass = 1.0 /
                (state_a.inv_mass + ra_ct.dot(state_a.inv_moment * ra_ct)
                 + state_b.inv_mass + rb_ct.dot(state_b.inv_moment * rb_ct));
            let lambda = -dv.dot(*tangent) * tangent_mass;
            let max_lambda = mix_friction * normal_impulse;
            let clamped = lambda.max(-max_lambda).min(max_lambda);
            impulse += clamped * tangent;
        }

        // Apply friction impulse
        obj_a.apply_impulse(-impulse, -state_a.inv_moment * ra.cross(impulse));
        obj_b.apply_impulse( impulse,  state_b.inv_moment * rb.cross(impulse));
    }

    /// Perform collision resolution from a manifold.
    /// This is the preferred method for resolving two bodies.
    pub fn resolve_manifold<Config, ObjA, ObjB>(
        obj_a: &mut ObjA,
        obj_b: &mut ObjB,
        manifold: &Manifold,
        dt: f32
    ) where Config: PhysicsConfig,
            ObjA: PhysicsObject,
            ObjB: PhysicsObject
    {
        for contact in &manifold.contacts {
            PhysicsState::resolve_contact::<Config, ObjA, ObjB>(obj_a, obj_b, contact, dt);
        }
        obj_a.update_dx();
        obj_b.update_dx();
    }
}

#[cfg(test)]
mod tests {
    mod physics {
        use cgmath::{Vector3, Point3, One, Zero};

        use geom::*;
        use physics::*;

        #[test]
        fn test_rigid_body() {
            let body_sample = RigidBody::new(1.0, 0.0, Vector3::zero(),
                                             vec![ Component::from(Sphere{ c: Point3::new(-5.0, 0.0, 0.0),
                                                                           r: 1.0 }),
                                                   Component::from(Sphere{ c: Point3::new(5.0, 0.0, 0.0),
                                                                           r: 1.0 }) ],
                                             vec![ 1.0, 1.0 ]);
            let mut body1 = body_sample.clone();
            let mut body2 = body1.clone();

            body1.set_pos(Point3::new(-10.0, 0.0, 0.0));
            body1.apply_impulse(Vector3::new(5.0, 0.0, 0.0) * 2.0, Vector3::zero());
            body1.integrate(1.0);
            body2.set_pos(Point3::new(10.0, 0.0, 0.0));
            body2.apply_impulse(Vector3::new(-5.0, 0.0, 0.0) * 2.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c| {
                manifold.push(c);
            }));
            
            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(-5.2005005, 0.0, 0.0));
            assert_eq!(body2.v, Vector3::new(5.2005005, 0.0, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            // Collisions with no penetration should produce no extra impulse besides slop:
            let mut body1 = body_sample.clone();
            let mut body2 = body_sample.clone();

            body1.set_pos(Point3::new(-10.0, 0.0, 0.0));
            body1.apply_impulse(Vector3::new(4.0, 0.0, 0.0) * 2.0, Vector3::zero());
            body1.integrate(1.0);
            body2.set_pos(Point3::new(10.0, 0.0, 0.0));
            body2.apply_impulse(Vector3::new(-4.0, 0.0, 0.0) * 2.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c| {
                manifold.push(c);
            }));
            
            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(-4.0004997, 0.0, 0.0));
            assert_eq!(body2.v, Vector3::new(4.0004997, 0.0, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            // Collisions with no penetration should produce no extra impulse besides slop:
            let mut body1 = body_sample.clone();
            let mut body2 = body_sample.clone();

            body1 += Vector3::new(0.0, 5.0, 0.0);
            body1.apply_impulse(Vector3::new(0.0, -4.0, 0.0) * 2.0, Vector3::zero());
            body1.integrate(1.0);
            body2 -= Vector3::new(0.0, 5.0, 0.0);
            body2.apply_impulse(Vector3::new(0.0, 4.0, 0.0) * 2.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c| {
                manifold.push(c);
            }));

            assert_eq!(manifold.contacts.len(), 2);
            
            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(0.0, 4.0639954, 0.0));
            assert_eq!(body2.v, Vector3::new(0.0, -4.0639954, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            // Test against a capsule RigidBody
            let capsule_sample = RigidBody::new(1.0, 0.0, Vector3::zero(),
                                             vec![ Component::from(Capsule{ a: Point3::new(-5.0, 0.0, 0.0),
                                                                            d: Vector3::new(10.0, 0.0, 0.0),
                                                                            r: 1.0 }) ],
                                             vec![ 2.0 ]);

            let mut body1 = capsule_sample.clone();
            let mut body2 = body_sample.clone();

            body1 += Vector3::new(0.0, 5.0, 0.0);
            body1.apply_impulse(Vector3::new(0.0, -4.0, 0.0) * 2.0, Vector3::zero());
            body1.integrate(1.0);
            body2 -= Vector3::new(0.0, 5.0, 0.0);
            body2.apply_impulse(Vector3::new(0.0, 4.0, 0.0) * 2.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c| {
                manifold.push(c);
            }));

            assert_eq!(manifold.contacts.len(), 2);
            
            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(0.0, 2.0769472, 0.0));
            assert_eq!(body2.v, Vector3::new(0.0, -2.0769472, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            // This collision should be exactly the same as the two spheres collision, except
            // with only one contact point.
            let mut body1 = capsule_sample.clone();
            let mut body2 = capsule_sample.clone();

            body1 += Vector3::new(0.0, 5.0, 0.0);
            body1.apply_impulse(Vector3::new(0.0, -4.0, 0.0) * 2.0, Vector3::zero());
            body1.integrate(1.0);
            body2 -= Vector3::new(0.0, 5.0, 0.0);
            body2.apply_impulse(Vector3::new(0.0, 4.0, 0.0) * 2.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c| {
                manifold.push(c);
            }));

            assert_eq!(manifold.contacts.len(), 1);
            
            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(0.0, 4.0004997, 0.0));
            assert_eq!(body2.v, Vector3::new(0.0, -4.0004997, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            // Test a raft-like object composed of four capsules.
            let raft = RigidBody::new(
                1.0,
                0.0,
                Vector3::zero(),
                vec![
                    Component::from(Capsule{
                        a: Point3::new(0.0, 0.0, 0.0),
                        d: Vector3::new(2.0, 0.0, 0.0),
                        r: 1.0,
                    }),
                    Component::from(Capsule{
                        a: Point3::new(2.0, 0.0, 0.0),
                        d: Vector3::new(0.0, 0.0, 2.0),
                        r: 1.0,
                    }),
                    Component::from(Capsule{
                        a: Point3::new(2.0, 0.0, 2.0),
                        d: Vector3::new(-2.0, 0.0, 0.0),
                        r: 1.0,
                    }),
                    Component::from(Capsule{
                        a: Point3::new(0.0, 0.0, 2.0),
                        d: Vector3::new(0.0, 0.0, -2.0),
                        r: 1.0,
                    })
                ],
                vec![ 1.0, 1.0, 1.0, 1.0 ]
            );

            let mut body1 = raft.clone();
            let mut body2 = raft.clone();

            body1 += Vector3::new(0.0, 5.0, 0.0);
            body1.apply_impulse(Vector3::new(0.0, -4.0, 0.0) * 4.0, Vector3::zero());
            body1.integrate(1.0);
            body2 -= Vector3::new(0.0, 5.0, 0.0);
            body2.apply_impulse(Vector3::new(0.0, 4.0, 0.0) * 4.0, Vector3::zero());
            body2.integrate(1.0);

            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c: LocalContact| {
                manifold.push(c);
            }));

            assert_eq!(manifold.contacts.len(), 4);

            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);
            assert_eq!(body1.v, Vector3::new(0.0, 8.338712, 0.0));
            assert_eq!(body2.v, Vector3::new(0.0, -8.338712, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
            assert_eq!(body2.collider.rot, Quaternion::one());

            let static_rect = Rect {
                c: Point3::new(0.0, 0.0, 0.0),
                u: [ Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0) ],
                e: [ 5.0, 5.0 ],
            };

            let mut body1 = raft.clone();
            let mut body2 = StaticBody::new(0.0, &static_rect);

            
            body1 += Vector3::new(0.0, 5.0, 0.0);
            body1.apply_impulse(Vector3::new(0.0, -4.0, 0.0) * 4.0, Vector3::zero());
            body1.integrate(1.0);

            assert_eq!(body1.v, Vector3::new(0.0, -4.0, 0.0));
 
            let mut manifold = Manifold::new();
            assert!(body1.local_contacts(&body2, |c: LocalContact| {
                manifold.push(c);
            }));

            assert_eq!(manifold.contacts.len(), 4);

            PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(&mut body1, &mut body2, &manifold, 1.0);

            // Unfortunately an excess of an impulse comes from too many contact points.
            assert_eq!(body1.v, Vector3::new(0.0, 8.337941, 0.0));
            assert_eq!(body1.collider.rot, Quaternion::one());
        }
    }
}
