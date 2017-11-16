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
use std::marker::PhantomData;
use cgmath::{EuclideanSpace, InnerSpace};

use smallvec::SmallVec;

use manifold::*;
use physics::*;

pub trait ConstrainedSet<Index, Constrained, Extra>
where
    Index: Copy
{
    fn get(&self, Index) -> (Constrained, Extra);
    fn set(&mut self, Index, Constrained);
}

pub trait Constraint {
    type Index: Copy;
    type Constrained;
    type Extra;

    /// Solve the constraint.
    fn solve<T: ConstrainedSet<Self::Index, Self::Constrained, Self::Extra>>(&mut self, &mut T);
}

pub struct Solver<C: Constraint> {
    constraints: SmallVec<[C; 10]>,
}

impl<C: Constraint> Solver<C> {
    pub fn new() -> Self {
        Solver {
            constraints: SmallVec::new(),
        }
    }

    pub fn add_constraint(&mut self, constraint: C) {
        self.constraints.push(constraint);
    }

    pub fn solve<T: ConstrainedSet<C::Index, C::Constrained, C::Extra>>(&mut self, cs: &mut T, iters: usize) {
        for _ in 0..iters {
            for constraint in self.constraints.iter_mut() {
                constraint.solve(cs);
            }
        }
    }
}

pub struct ContactConstraint<Index, Params = DefaultContactConstraintParams>
where
    Index: Copy,
    Params: ContactConstraintParams
{
    obj_a: Index,
    obj_b: Index,
    manifold: Manifold,
    friction: f32,
    states: SmallVec<[ContactState; 4]>,
    params: PhantomData<Params>,
}

impl<Index, Params> ContactConstraint<Index, Params>
where
    Index: Copy,
    Params: ContactConstraintParams
{
    pub fn new<T: ConstrainedSet<Index, Velocity, RigidBodyInfo>>(pool: &T, obj_a: Index, obj_b: Index, manifold: Manifold, dt: f32) -> Self {
        let (
            Velocity { linear: va, angular: oa },
            RigidBodyInfo {
                x: xa,
                restitution: rest_a,
                friction: fric_a,
                inv_mass: inv_mass_a,
                inv_moment: inv_moment_a
            }
        ) = pool.get(obj_a);

        let (
            Velocity { linear: vb, angular: ob },
            RigidBodyInfo {
                x: xb,
                restitution: rest_b,
                friction: fric_b,
                inv_mass: inv_mass_b,
                inv_moment: inv_moment_b
            }
        ) = pool.get(obj_b);

        // Mix restitution and friction values:
        let restitution = rest_a.max(rest_b);
        let friction = (fric_a * fric_b).sqrt();

        // Calculate contact states for each contact
        let mut states = SmallVec::with_capacity(manifold.contacts.len());
        for &(local_a, local_b) in manifold.contacts.iter() {
            let ra = local_a.to_vec();
            let rb = local_b.to_vec();
            let ca = ra + xa.to_vec();
            let cb = rb + xb.to_vec();
            let ra_cn = ra.cross(manifold.normal);
            let rb_cn = rb.cross(manifold.normal);

            // Penetration is defined as the distance between the two contact points
            // dotted with the normal vector.
            let pen = (cb - ca).dot(manifold.normal);

            let dv = vb + ob.cross(rb) - va - oa.cross(ra);
            let rel_v = dv.dot(manifold.normal);

            let bias = -Params::BAUMGARTE / dt * if pen > 0.0 {
                0.0
            } else {
                pen + Params::PENETRATION_SLOP
            } + if rel_v < -1.0 {
                -restitution * rel_v
            } else {
                0.0
            };

            let normal_mass = 1.0 /
                (inv_mass_a + ra_cn.dot(inv_moment_a * ra_cn)
                 + inv_mass_b + rb_cn.dot(inv_moment_b * rb_cn));

            let tangent_mass = [
                {
                    let ra_ct = ra.cross(manifold.tangent_vector[0]);
                    let rb_ct = rb.cross(manifold.tangent_vector[0]);
                    1.0 / (inv_mass_a + ra_ct.dot(inv_moment_a * ra_ct)
                           + inv_mass_b + rb_ct.dot(inv_moment_b * rb_ct))
                },
                {
                    let ra_ct = ra.cross(manifold.tangent_vector[1]);
                    let rb_ct = rb.cross(manifold.tangent_vector[1]);
                    1.0 / (inv_mass_a + ra_ct.dot(inv_moment_a * ra_ct)
                           + inv_mass_b + rb_ct.dot(inv_moment_b * rb_ct))
                }
            ];
            
            states.push(ContactState {
                bias,
                normal_mass,
                normal_impulse: 0.0,
                tangent_mass,
                tangent_impulse: [ 0.0, 0.0 ]
            })
        }

        ContactConstraint {
            obj_a,
            obj_b,
            manifold,
            friction,
            states,
            params: PhantomData,
        }
    }
}

impl<Index, Params> Constraint for ContactConstraint<Index, Params>
where
    Index: Copy,
    Params: ContactConstraintParams
{
    type Index = Index;
    type Constrained = Velocity;
    type Extra = RigidBodyInfo;

    fn solve<T: ConstrainedSet<Index, Velocity, RigidBodyInfo>>(&mut self, pool: &mut T) {
        let (
            Velocity{ linear: mut va, angular: mut oa },
            RigidBodyInfo{ inv_mass: inv_mass_a, inv_moment: inv_moment_a, .. }
        ) = pool.get(self.obj_a);

        let (
            Velocity{ linear: mut vb, angular: mut ob },
            RigidBodyInfo{ inv_mass: inv_mass_b, inv_moment: inv_moment_b, .. }
        ) = pool.get(self.obj_b);

        for (i, ref mut contact_state) in self.states.iter_mut().enumerate() {
            let (local_a, local_b) = self.manifold.contacts[i];
            let (ra, rb) = (local_a.to_vec(), local_b.to_vec());
            let dv = vb + ob.cross(rb) - va - oa.cross(ra);

            // Calculate friction impulse
            for i in 0..2 {
                let lambda = -dv.dot(self.manifold.tangent_vector[i]) * contact_state.tangent_mass[i];
                let max_lambda = self.friction * contact_state.normal_impulse;
                let prev_impulse = contact_state.tangent_impulse[i];
                contact_state.tangent_impulse[i] =
                    clamp(-max_lambda, max_lambda, prev_impulse + lambda);
                let impulse = self.manifold.tangent_vector[i] * lambda;
                va -= impulse * inv_mass_a;
                oa -= inv_moment_a * ra.cross(impulse);
                vb += impulse * inv_mass_b;
                ob += inv_moment_b * rb.cross(impulse);
            }

            let dv = vb + ob.cross(rb) - va - oa.cross(ra);
            // Calculate normal impulse
            let vn = dv.dot(self.manifold.normal);
            let lambda = contact_state.normal_mass * (-vn + contact_state.bias);
            let prev_impulse = contact_state.normal_impulse;
            contact_state.normal_impulse = (prev_impulse + lambda).max(0.0);
            let lambda = contact_state.normal_impulse - prev_impulse;

            // Apply normal impulse
            let impulse = self.manifold.normal * lambda;
            va -= impulse * inv_mass_a;
            oa -= inv_moment_a * ra.cross(impulse);
            vb += impulse * inv_mass_b;
            ob += inv_moment_b * rb.cross(impulse);
            
        }
        pool.set(self.obj_a, Velocity{ linear: va, angular: oa });
        pool.set(self.obj_b, Velocity{ linear: vb, angular: ob });
    }

}

struct ContactState {
    bias: f32,
    normal_mass: f32,
    normal_impulse: f32,
    tangent_mass: [f32; 2],
    tangent_impulse: [f32; 2],
}


/// A type that describes parameters used when solving contact constraints.
pub trait ContactConstraintParams {
    const PENETRATION_SLOP: f32;
    const BAUMGARTE: f32;
}

/// The suggested set of parameters to use when resolving collisions.
pub struct DefaultContactConstraintParams {}

/// A contact constraint solver for PhysicsObjects.
///
/// This is an unsafe interface in a lot of instances.
impl ContactConstraintParams for DefaultContactConstraintParams {
    const PENETRATION_SLOP: f32 = 0.05;
    const BAUMGARTE: f32 = 0.2;
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
