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

struct ContactState {
    bias: f32,
    normal_mass: f32,
    normal_impulse: f32,
    tangent_mass: [f32; 2],
    tangent_impulse: [f32; 2],
}

/// A constraint created by two objects coming in contact.
pub struct ContactConstraint<'a> {
    /// Mixed restitution of the two objects
    pub friction: f32,
    pub obj_a: &'a mut PhysicsObject,
    pub obj_b: &'a mut PhysicsObject,
    manifold: Manifold,
    states: SmallVec<[ContactState; 4]>,
}

/// A type that describes parameters used when solving contact constraints.
pub trait ContactSolverParams {
    const PENETRATION_SLOP: f32;
    const BAUMGARTE: f32;
}

/// The suggested set of parameters to use when resolving collisions.
pub struct DefaultContactSolverParams {}

/// A contact constraint solver for PhysicsObjects.
///
/// This is an unsafe interface in a lot of instances.
impl ContactSolverParams for DefaultContactSolverParams {
    const PENETRATION_SLOP: f32 = 0.05;
    const BAUMGARTE: f32 = 0.2;
}

/// A PGE based constraint solver for contacts. 
pub struct ContactSolver<'a, Params = DefaultContactSolverParams>
where
    Params: ContactSolverParams
{
    constraints: SmallVec<[ContactConstraint<'a>; 20]>,
    params: PhantomData<Params>,
}

impl<'a, Params: ContactSolverParams> ContactSolver<'a, Params> {
    /// Construct a new empty solver
    pub fn new() -> Self {
        ContactSolver {
            constraints: SmallVec::new(),
            params: PhantomData,
        }
    }

    /// Construct a new empty solver with a given capacity
    pub fn with_capacity(cap: usize) -> Self {
        ContactSolver {
            constraints: SmallVec::with_capacity(cap),
            params: PhantomData
        }
    }

    /// Add a constraint to the solver and pre-solve the constraint.
    ///
    /// Unfortunately this interface requires taking two mutable references, so
    /// it's hard to work with this interface safely. 
    pub fn add_constraint<'b, 'c>(
        &mut self,
        obj_a: &'b mut PhysicsObject,
        obj_b: &'c mut PhysicsObject,
        manifold: Manifold,
        dt: f32
    ) where
        'b: 'a,
        'c: 'a,
    {
        if manifold.contacts.len() == 0 {
            // Nothing to do.
            return;
        }

        // Pre-solve the constraint:
        let (xa, xb) = (obj_a.pos(), obj_b.pos());
        let Velocity{ linear: va, angular: oa } = obj_a.get_dx();
        let Velocity{ linear: vb, angular: ob } = obj_b.get_dx();

        let inv_mass_a = obj_a.inv_mass();
        let inv_moment_a = obj_a.inv_moment();
        let inv_mass_b = obj_b.inv_mass();
        let inv_moment_b = obj_b.inv_moment();

        // Mix restitution and friction values:
        let restitution = obj_a.restitution().max(obj_b.restitution());
        let friction = (obj_a.friction() * obj_b.friction()).sqrt();

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
        
        self.constraints.push(ContactConstraint {
            friction,
            obj_a,
            obj_b,
            manifold,
            states,
        });
    }

    /// Solves the added constraints and updates the inserted physics objects
    /// with the correct new velocities.
    pub fn solve(&mut self, iters: usize) {
        for _ in 0..iters {
            for &mut ContactConstraint{
                friction,
                ref mut obj_a, ref mut obj_b,
                ref manifold,
                ref mut states,
            } in self.constraints.iter_mut() {
                let Velocity{ linear: mut va, angular: mut oa } = obj_a.get_dx();
                let Velocity{ linear: mut vb, angular: mut ob } = obj_b.get_dx();
                let inv_mass_a = obj_a.inv_mass();
                let inv_moment_a = obj_a.inv_moment();
                let inv_mass_b = obj_b.inv_mass();
                let inv_moment_b = obj_b.inv_moment();
                for (i, ref mut contact_state) in states.iter_mut().enumerate() {
                    let (local_a, local_b) = manifold.contacts[i];
                    let (ra, rb) = (local_a.to_vec(), local_b.to_vec());
                    let dv = vb + ob.cross(rb) - va - oa.cross(ra);

                    // Calculate friction impulse
                    for i in 0..2 {
                        let lambda = -dv.dot(manifold.tangent_vector[i]) * contact_state.tangent_mass[i];
                        let max_lambda = friction * contact_state.normal_impulse;
                        let prev_impulse = contact_state.tangent_impulse[i];
                        contact_state.tangent_impulse[i] =
                            clamp(-max_lambda, max_lambda, prev_impulse + lambda);
                        let impulse = manifold.tangent_vector[i] * lambda;
                        va -= impulse * inv_mass_a;
                        oa -= inv_moment_a * ra.cross(impulse);
                        vb += impulse * inv_mass_b;
                        ob += inv_moment_b * rb.cross(impulse);
                    }

                    let dv = vb + ob.cross(rb) - va - oa.cross(ra);
                    // Calculate normal impulse
                    let vn = dv.dot(manifold.normal);
                    let lambda = contact_state.normal_mass * (-vn + contact_state.bias);
                    let prev_impulse = contact_state.normal_impulse;
                    contact_state.normal_impulse = (prev_impulse + lambda).max(0.0);
                    let lambda = contact_state.normal_impulse - prev_impulse;

                    // Apply normal impulse
                    let impulse = manifold.normal * lambda;
                    va -= impulse * inv_mass_a;
                    oa -= inv_moment_a * ra.cross(impulse);
                    vb += impulse * inv_mass_b;
                    ob += inv_moment_b * rb.cross(impulse);
                   
                }
                obj_a.set_dx(Velocity{ linear: va, angular: oa });
                obj_b.set_dx(Velocity{ linear: vb, angular: ob });
            }
        }
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
