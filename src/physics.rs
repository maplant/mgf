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

/// An type that has a moment of inertia tensor.
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
/// A CompoundDynamicBody is technically a shape, although its geometries are all considered
/// to be in motion.
#[derive(Clone)]
pub struct CompoundDynamicBody {
    pub restitution: f32,
    pub friction: f32,
    pub total_inv_mass: f32,
    pub inv_moment_body: Matrix3<f32>,
    pub inv_moment: Matrix3<f32>,
    pub v: Vector3<f32>,
    pub omega: Vector3<f32>,
    pub v_step: Vector3<f32>,
    pub force: Vector3<f32>,
    pub torque: Vector3<f32>,
    // Position and rotation are stored in collider.
    pub collider: Compound,
    pub component_masses: Vec<f32>,
}

impl CompoundDynamicBody {
    /// Construct a new CompoundDynamicBody from a Vec of Components and masses.
    pub fn new(
        restitution: f32,
        friction: f32,
        world_force: Vector3<f32>,
        mut shapes: Vec<Component>,
        masses: Vec<f32>
    ) -> Self {
        // Calculate the center of mass of the CompoundDynamicBody.
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
        CompoundDynamicBody {
            restitution,
            friction,
            total_inv_mass: inv_mass,
            inv_moment_body: inv_moment,
            inv_moment,
            v: Vector3::zero(),
            omega: Vector3::zero(),
            v_step: Vector3::zero(),
            force: world_force / inv_mass,
            torque: Vector3::zero(),
            collider: Compound::new(shapes),
            component_masses: masses,
        }
    }
}    

impl AddAssign<Vector3<f32>> for CompoundDynamicBody {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.collider.disp += v
    }
}

impl SubAssign<Vector3<f32>> for CompoundDynamicBody {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.collider.disp -= v
    }
}

impl Shape for CompoundDynamicBody {
    /// The center returned by a CompoundDynamicBody is the center of mass for that
    /// object. 
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.collider.disp)
    }
}

impl Delta for CompoundDynamicBody {
    fn delta(&self) -> Vector3<f32> {
        self.v_step
    }
}

impl BoundedBy<AABB> for CompoundDynamicBody {
    fn bounds(&self) -> AABB {
        let b1 = self.collider.bounds();
        let b2 = b1 + self.v_step;
        Bound::combine(&b1, &b2)
    }
}

impl BoundedBy<Sphere> for CompoundDynamicBody {
    fn bounds(&self) -> Sphere {
        let b1 = self.collider.bounds();
        let b2 = b1 + self.v_step;
        Bound::combine(&b1, &b2)
    }
}

impl<RHS> Contacts<RHS> for CompoundDynamicBody
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
                comp.rotate_about(self.collider.rot,
                                  Point3::new(0.0, 0.0, 0.0)) + self.collider.disp,
                self.v_step
            );
            rhs.contacts(&shape, |c| { collided = true; callback(-c) });
        });
        collided
    }
}

/// A dynamic body based on a single component.
#[derive(Clone)]
pub struct SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    pub restitution: f32,
    pub friction: f32,
    pub inv_mass: f32,
    pub inv_moment_body: Matrix3<f32>,
    pub inv_moment: Matrix3<f32>,
    pub v: Vector3<f32>,
    pub omega: Vector3<f32>,
    pub q: Quaternion<f32>,
    pub force: Vector3<f32>,
    pub torque: Vector3<f32>,
    pub collider: Moving<S>,
}


impl<S> SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    /// Construct a new SimpleDynamicBody.
    pub fn new(
        restitution: f32,
        friction: f32,
        world_force: Vector3<f32>,
        shape: S,
        mass: f32
    ) -> Self {
        let inv_mass = 1.0 / mass;
        let inv_moment = shape.tensor(mass).invert().unwrap();
        SimpleDynamicBody {
            restitution,
            friction,
            inv_mass,
            inv_moment_body: inv_moment,
            inv_moment,
            v: Vector3::zero(),
            omega: Vector3::zero(),
            q: Quaternion::one(),
            force: world_force / inv_mass,
            torque: Vector3::zero(),
            collider: Moving::sweep(shape, Vector3::zero()),
        }
    }
}    

impl<S> AddAssign<Vector3<f32>> for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    fn add_assign(&mut self, v: Vector3<f32>) {
        *self.collider.as_mut() += v
    }
}

impl<S> SubAssign<Vector3<f32>> for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    fn sub_assign(&mut self, v: Vector3<f32>) {
        *self.collider.as_mut() -= v
    }
}

impl<S> Shape for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    fn center(&self) -> Point3<f32> {
        self.collider.as_ref().center()
    }
}


impl<S> Delta for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    fn delta(&self) -> Vector3<f32> {
        self.collider.delta()
    }
}

impl<S, B> BoundedBy<B> for SimpleDynamicBody<S>
where
    S: BoundedBy<B> + Volumetric + Inertia + Copy,
    B: Bound
{
    fn bounds(&self) -> B {
        self.collider.bounds()
    }
}

impl<S, RHS> Contacts<RHS> for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy, 
    RHS: Contacts<Moving<S>>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        rhs.contacts(&self.collider, |c|callback(-c))
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
                local_a: c.a + -a_c.to_vec(),
                local_b: c.b + -b_c.to_vec(),
                global: c
            })
        })
    }
}


/// A description of the positional derivative of an object.
#[derive(Copy, Clone)]
pub struct Velocity {
    pub linear: Vector3<f32>,
    pub angular: Vector3<f32>,
}

/// A type that exhibits physical properties.
pub trait PhysicsObject {
    /// Integrate the object over the time step.
    fn integrate(&mut self, dt: f32);

    fn get_dx(&self) -> Velocity;

    fn set_dx(&mut self, v: Velocity);

    /// Return the position of the object at the end of the time step.
    fn pos(&self) -> Point3<f32>;

    fn inv_mass(&self) -> f32;

    fn inv_moment(&self) -> Matrix3<f32>;

    fn restitution(&self) -> f32;

    fn friction(&self) -> f32;
}

impl PhysicsObject for CompoundDynamicBody {
    fn integrate(&mut self, dt: f32) {
        self.collider += self.v_step;
        self.collider.rot = (self.collider.rot
                             + Quaternion::from_sv(0.0, self.omega)
                             * 0.5 * self.collider.rot).normalize();
        let r = Matrix3::from(self.collider.rot);
        self.inv_moment = r * self.inv_moment_body * r.transpose();
        self.v += self.force * self.total_inv_mass * dt;
        self.omega += self.inv_moment * self.torque * dt;
        self.v_step = self.v * dt;

        self.v *= 1.0 / (1.0 + dt * 0.00001);
        self.omega *= 1.0 / (1.0 + dt * 0.1);
    }

    fn get_dx(&self) -> Velocity {
        Velocity {
            linear: self.v,
            angular: self.omega
        }
    }

    fn set_dx(&mut self, v: Velocity) {
        self.v = v.linear;
        self.omega = v.angular;
    }

    fn pos(&self) -> Point3<f32> {
        self.collider.center() + self.v_step
    }

    fn inv_mass(&self) -> f32 {
        self.total_inv_mass
    }

    fn inv_moment(&self) -> Matrix3<f32> {
        self.inv_moment
    }

    fn restitution(&self) -> f32 {
        self.restitution
    }

    fn friction(&self) -> f32 {
        self.friction
    }
}

impl<S> PhysicsObject for SimpleDynamicBody<S>
where
    S: Volumetric + Inertia + Copy
{
    fn integrate(&mut self, dt: f32) {
        *self.collider.as_mut() += self.collider.1;
        self.q = (self.q
                  + Quaternion::from_sv(0.0, self.omega)
                  * 0.5 * self.q).normalize();
        let r = Matrix3::from(self.q);
        self.inv_moment = r * self.inv_moment_body * r.transpose();
        self.v += self.inv_mass * self.force * dt;
        self.omega += self.inv_moment * self.torque * dt;
        self.collider.1 = self.v * dt;

        self.v *= 1.0 / (1.0 + dt * 0.00001);
        self.omega *= 1.0 / (1.0 + dt * 0.1);
    }

    fn get_dx(&self) -> Velocity {
        Velocity {
            linear: self.v,
            angular: self.omega
        }
    }

    fn set_dx(&mut self, v: Velocity) {
        self.v = v.linear;
        self.omega = v.angular;
    }

    fn pos(&self) -> Point3<f32> {
        self.collider.as_ref().center() + self.collider.1
    }

    fn inv_mass(&self) -> f32 {
        self.inv_mass
    }

    fn inv_moment(&self) -> Matrix3<f32> {
        self.inv_moment
    }

    fn restitution(&self) -> f32 {
        self.restitution
    }

    fn friction(&self) -> f32 {
        self.friction
    }
}

impl<'a, S: Shape + 'a> PhysicsObject for StaticBody<'a, S> {
    #[inline(always)]
    fn integrate(&mut self, _dt: f32) {}

    #[inline(always)]
    fn get_dx(&self) -> Velocity {
        Velocity {
            linear: Vector3::zero(),
            angular: Vector3::zero(),
        }
    }

    #[inline(always)]
    fn set_dx(&mut self, _v: Velocity) {}

    fn pos(&self) -> Point3<f32> { self.shape.center() }

    fn inv_mass(&self) -> f32 {
        0.0
    }

    fn inv_moment(&self) -> Matrix3<f32> {
        Matrix3::zero()
    }

    fn restitution(&self) -> f32 {
        0.0
    }

    fn friction(&self) -> f32 {
        0.0
    }
}
