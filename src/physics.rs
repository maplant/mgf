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
use std::slice::Iter;
use cgmath::{InnerSpace, SquareMatrix, EuclideanSpace, Matrix, Matrix3, Point3,
             Quaternion, Vector3, Zero, One};

use crate::compound::*;
use crate::geom::*;
use crate::solver::*;

/// Any type that has a moment of inertia tensor.
pub trait Inertia {
    fn tensor(&self, m: f32) -> Matrix3<f32>;
}

impl Inertia for Sphere {
    fn tensor(&self, m: f32) -> Matrix3<f32> {
        let i = 0.4 * m * self.r * self.r;
        let i = Matrix3::new(
            i, 0.0, 0.0,
            0.0, i, 0.0,
            0.0, 0.0, i
        );
        let disp = self.c.to_vec();
        let outer = Matrix3::from_cols(
            disp * disp.x,
            disp * disp.y,
            disp * disp.z
        );
        i + m * (Matrix3::one() * disp.dot(disp) - outer)
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
        let i = rot * Matrix3::new(
            i_x, 0.0, 0.0,
            0.0, i_y, 0.0,
            0.0, 0.0, i_z
        ) * rot.transpose();
        let disp = self.center().to_vec();
        let outer = Matrix3::from_cols(
            disp * disp.x,
            disp * disp.y,
            disp * disp.z
        );
        i + m * (Matrix3::one() * disp.dot(disp) - outer)
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

impl Inertia for OBB {
    fn tensor(&self, m: f32) -> Matrix3<f32> {
        // Thank you wikipedia
        let (x, y, z) = (
            self.r.x * 2.0,
            self.r.y * 2.0,
            self.r.z * 2.0
        );
        let i_x = 1.0 / 12.0 * m * (y * y + z * z);
        let i_y = 1.0 / 12.0 * m * (x * x + z * z);
        let i_z = 1.0 / 12.0 * m * (x * x + y * y);
        let rot = Matrix3::<f32>::from(self.q);
        let i = rot * Matrix3::new(
            i_x, 0.0, 0.0,
            0.0, i_y, 0.0,
            0.0, 0.0, i_z
        ) * rot.transpose();
        let disp = self.center().to_vec();
        let outer = Matrix3::from_cols(
            disp * disp.x,
            disp * disp.y,
            disp * disp.z
        );
        i + m * (Matrix3::one() * disp.dot(disp) - outer)
    }
}

/// A description of the physical state of an object minus linear and angular
/// velocity.
pub struct RigidBodyInfo {
    pub x: Point3<f32>,
    pub restitution: f32,
    pub friction: f32,
    pub inv_mass: f32,
    pub inv_moment: Matrix3<f32>,
}

/// A description of the positional derivative of an object.
#[derive(Copy, Clone)]
pub struct Velocity {
    pub linear: Vector3<f32>,
    pub angular: Vector3<f32>,
}

/// A vector of rigid bodies. 
#[derive(Clone)]
pub struct RigidBodyVec {
    pub x: Vec<Point3<f32>>,
    pub q: Vec<Quaternion<f32>>,
    v: Vec<Vector3<f32>>,
    omega: Vec<Vector3<f32>>,
    force: Vec<Vector3<f32>>,
    torque: Vec<Vector3<f32>>,
    restitution: Vec<f32>,
    friction: Vec<f32>,
    inv_mass: Vec<f32>,
    inv_moment_body: Vec<Matrix3<f32>>,
    inv_moment: Vec<Matrix3<f32>>,
    constructor: Vec<ComponentConstructor>,
    pub collider: Vec<Moving<Component>>,
}

/// A reference to an element of a RigidBodyVec.
#[derive(Copy, Clone)]
pub enum RigidBodyRef {
    Dynamic(usize),
    Static{ center: Point3<f32>, friction: f32 },
}

impl From<usize> for RigidBodyRef {
    fn from(i: usize) -> RigidBodyRef {
        RigidBodyRef::Dynamic(i)
    }
}

impl Into<usize> for RigidBodyRef {
    fn into(self) -> usize {
        match self {
            RigidBodyRef::Dynamic(i) => i,
            _ => panic!("not stored"),
        }
    }
}

impl RigidBodyVec {
    /// Create an empty RigidBodyVec.
    pub fn new() -> Self {
        RigidBodyVec {
            x: Vec::new(),
            q: Vec::new(),
            v: Vec::new(),
            omega: Vec::new(),
            force: Vec::new(),
            torque: Vec::new(),
            restitution: Vec::new(),
            friction: Vec::new(),
            inv_mass: Vec::new(),
            inv_moment_body: Vec::new(),
            inv_moment: Vec::new(),
            constructor: Vec::new(),
            collider: Vec::new(),
        }
    }

    /// Add a body to the rigid body. 
    pub fn add_body(&mut self, collider: Component, mass: f32, restitution: f32, friction: f32, world_force: Vector3<f32>) -> RigidBodyRef {
        let id = self.x.len();
        let (x, q, constructor) = collider.deconstruct();
        self.x.push(x);
        self.q.push(q);
        self.v.push(Vector3::zero());
        self.omega.push(Vector3::zero());
        self.force.push(world_force * mass);
        self.torque.push(Vector3::zero());
        self.restitution.push(restitution);
        self.friction.push(friction);
        self.inv_mass.push(1.0 / mass);
        let inv_moment = (collider - x.to_vec()).tensor(mass).invert().unwrap();
        self.inv_moment_body.push(inv_moment);
        self.inv_moment.push(inv_moment);
        self.constructor.push(constructor);
        self.collider.push(Moving::sweep(collider, Vector3::zero()));
        RigidBodyRef::Dynamic(id)
    }

    /// Calculate the rotation, velocity, tensor, and collider for each rigid
    /// body.
    pub fn integrate(&mut self, dt: f32) {
        unsafe { // ONLY for get_unchecked
            // Update rotation:
            for (i, q) in self.q.iter_mut().enumerate() {
                *q = (*q + Quaternion::from_sv(0.0, *self.omega.get_unchecked(i) * dt)
                     * 0.5 * *q).normalize();
            }
            // Update moment:
            for (i, moment) in self.inv_moment.iter_mut().enumerate() {
                let r = Matrix3::from(*self.q.get_unchecked(i));
                *moment = r * *self.inv_moment_body.get_unchecked(i) * r.transpose();
            }
            // Update linear velocity:
            for (i, v) in self.v.iter_mut().enumerate() {
                *v += *self.force.get_unchecked(i) * *self.inv_mass.get_unchecked(i) * dt;
            }
            // Update angular velocity:
            for (i, omega) in self.omega.iter_mut().enumerate() {
                *omega += *self.inv_moment.get_unchecked(i) * self.torque.get_unchecked(i) * dt;
            }
            // Update colliders:
            for (i, collider) in self.collider.iter_mut().enumerate() {
                *collider = Moving::sweep(
                    self.constructor.get_unchecked(i).construct(
                        *self.x.get_unchecked(i),
                        *self.q.get_unchecked(i)
                    ),
                    *self.v.get_unchecked(i) * dt
                );
            }
        }
    }

    /// Return an iterator for the colliders of the rigid body.
    pub fn colliders(&self) -> Iter<Moving<Component>> {
        self.collider.iter()
    }

    /// Finish all motion that was initiated at the beginning of the frame by
    /// integrate.
    pub fn complete_motion(&mut self) {
        unsafe {
            // Update position from last time step:
            for (i, x) in self.x.iter_mut().enumerate() {
                *x += self.collider.get_unchecked(i).delta();
            }
        }
    }
}

impl ConstrainedSet<RigidBodyRef, Velocity, RigidBodyInfo> for RigidBodyVec {
    fn get(&self, i: RigidBodyRef) -> (Velocity, RigidBodyInfo) {
        match i {
            RigidBodyRef::Dynamic(i) =>
                (
                    Velocity {
                        linear: self.v[i],
                        angular: self.omega[i],
                    },
                    RigidBodyInfo {
                        x: self.x[i] + self.collider[i].delta(),
                        restitution: self.restitution[i],
                        friction: self.friction[i],
                        inv_mass: self.inv_mass[i],
                        inv_moment: self.inv_moment[i],
                    }
                ),
            RigidBodyRef::Static{ center, friction } =>
                (
                    Velocity {
                        linear: Vector3::zero(),
                        angular: Vector3::zero(), 
                    },
                    RigidBodyInfo {
                        x: center,
                        restitution: 0.0,
                        friction, 
                        inv_mass: 0.0, 
                        inv_moment: Matrix3::zero(),
                    }
                ),
        }
    }

    fn set(&mut self, i: RigidBodyRef, v: Velocity) {
        match i {
            RigidBodyRef::Dynamic(i) => {
                self.v[i] = v.linear;
                self.omega[i] = v.angular;
            },
            RigidBodyRef::Static{ .. } => (),
        }
    }
}

#[cfg(test)]
mod tests {
    mod inertia {
        #[test]
        fn test_tensors() {
            use cgmath::{Point3, Matrix3};
            use crate::geom::Sphere;
            use crate::physics::Inertia;

            let s = Sphere{ c: Point3::new(0.0, 0.0, 0.0), r: 1.0 };
            assert_eq!(
                s.tensor(1.0),
                Matrix3::new(
                    0.4, 0.0, 0.0,
                    0.0, 0.4, 0.0,
                    0.0, 0.0, 0.4
                )
            );
        }
    }
}
            
