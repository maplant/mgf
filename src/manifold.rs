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

use cgmath::prelude::*;
use cgmath::{Vector3, Point3};

use smallvec::SmallVec;
use collision::*;
use geom::*;

/// A type that supplies constants to a ContactPruner.
pub trait PruningParams {
    /// Minimum distance between two points required to not reject a point.
    const PERSISTENT_THRESHOLD_SQ: f32;

    // const MAX_NUMBER_OF_CONTACTS: usize
}

/// Default parameters supplied for pruning
pub struct DefaultPruningParams {}

impl PruningParams for DefaultPruningParams {
    const PERSISTENT_THRESHOLD_SQ: f32 = 0.5;
}

/// Structure for pruning unnecessary contact points.
pub struct ContactPruner<Params = DefaultPruningParams>
where
    Params: PruningParams
{
    min_col_time: f32,
    contacts: SmallVec<[LocalContact; 4]>,
    params: PhantomData<Params>,
}

impl<Params: PruningParams> ContactPruner<Params> {
    pub fn new() -> Self {
        ContactPruner {
            min_col_time: f32::INFINITY,
            contacts: SmallVec::new(),
            params: PhantomData,
        }
    }

    pub fn with_capacity(cap: usize) -> Self {
        ContactPruner {
            min_col_time: f32::INFINITY,
            contacts: SmallVec::with_capacity(cap),
            params: PhantomData,
        }
    }

    /// Determines if the contact given is far away from all of the current
    /// contact points stored in the Manifold, and if it is pushes it into the vec.
    /// If two contact points are too close to each other, the furthest from the
    /// center of the geometries is chosen. 
    pub fn push(&mut self, new_contact: LocalContact) {
        if new_contact.global.t < self.min_col_time - COLLISION_EPSILON {
            // If the collision occurs earlier then the collisions in this manifold,
            // it is the preferred collision.
            self.contacts.clear();
            self.contacts.push(new_contact);
            self.min_col_time = new_contact.global.t;
            return;
        } else if new_contact.global.t > self.min_col_time + COLLISION_EPSILON {
            return;
        }
        for old_contact in self.contacts.iter_mut() {
            let ra = new_contact.global.a - old_contact.global.a;
            let rb = new_contact.global.b - old_contact.global.b;
            if ra.magnitude2() <= Params::PERSISTENT_THRESHOLD_SQ ||
                rb.magnitude2() <= Params::PERSISTENT_THRESHOLD_SQ
            {
                // Use the contact that appears to increase the distance between
                // the contact points and the object's center of mass.
                let prev_dist = old_contact.local_a.to_vec().magnitude2() +
                    old_contact.local_b.to_vec().magnitude2();
                let new_dist = new_contact.local_a.to_vec().magnitude2() +
                    new_contact.local_b.to_vec().magnitude2();
                if prev_dist < new_dist {
                    *old_contact = new_contact;
                }
                return;
            }
        }
        self.contacts.push(new_contact);
    }

    pub fn clear(&mut self) {
        self.min_col_time = f32::INFINITY;
        self.contacts.clear();
    }
}

/// A set of contacts between two objects.  
#[derive(Clone, Debug)]
pub struct Manifold {
    pub time: f32,
    pub normal: Vector3<f32>,
    pub tangent_vector: [Vector3<f32>; 2],
    /// List of the local contact points.
    pub contacts: SmallVec<[(Point3<f32>, Point3<f32>); 4]>,
}

impl From<LocalContact> for Manifold {
    fn from(lc: LocalContact) -> Self {
        Manifold {
            time: lc.global.t,
            normal: lc.global.n,
            tangent_vector: compute_basis(&lc.global.n),
            contacts: SmallVec::from_vec(vec![(lc.local_a, lc.local_b)]),
        }
    }
}

impl<P: PruningParams> From<ContactPruner<P>> for Manifold {
    fn from(pruner: ContactPruner<P>) -> Self {
        let mut contacts: SmallVec<[(Point3<f32>, Point3<f32>); 4]>
            = SmallVec::with_capacity(pruner.contacts.len());
        let avg_normal = pruner.contacts.iter()
            .fold(Vector3::zero(),
                  |sum, lc| {
                      contacts.push((lc.local_a, lc.local_b));
                      sum + lc.global.n
                  }) / (pruner.contacts.len() as f32);
        Manifold {
            time: pruner.min_col_time,
            normal: avg_normal,
            tangent_vector: compute_basis(&avg_normal),
            contacts
        }
    }
}


impl Manifold {
    pub fn new() -> Self {
        Manifold {
            time: 0.0,
            normal: Vector3::zero(),
            tangent_vector: [Vector3::zero(), Vector3::zero()],
            contacts: SmallVec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.contacts.len()
    }
}

