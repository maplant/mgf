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
use cgmath::{EuclideanSpace, InnerSpace};

use smallvec::SmallVec;
use collision::*;
use geom::*;

/// Minimum distance between two points required to not reject a point.
pub const PERSISTENT_THRESHOLD_SQ: f32 = 1.0;

/// Structure for pruning unnecessary contact points.
#[derive(Clone, Debug)]
pub struct Manifold {
    pub min_coll_time: f32,
    pub contacts: SmallVec<[LocalContact; 4]>,
}

impl Manifold {
    pub fn new() -> Self {
        Manifold {
            min_coll_time: f32::INFINITY,
            contacts: SmallVec::new(),
        }
    }

    pub fn with_capacity(cap: usize) -> Self {
        Manifold {
            min_coll_time: f32::INFINITY,
            contacts: SmallVec::with_capacity(cap),
        }
    }

    /// Determines if the contact given is far away from all of the current
    /// contact points stored in the Manifold, and if it is pushes it into the vec.
    /// If two contact points are too close to each other, the furthest from the
    /// center of the geometries is chosen. 
    pub fn push(&mut self, new_contact: LocalContact) {
        if new_contact.global.t < self.min_coll_time - COLLISION_EPSILON {
            // If the collision occurs earlier then the collisions in this manifold,
            // it is the preferred collision.
            self.contacts.clear();
            self.contacts.push(new_contact);
            self.min_coll_time = new_contact.global.t;
            return;
        } else if new_contact.global.t > self.min_coll_time + COLLISION_EPSILON {
            return;
        }
        for old_contact in self.contacts.iter_mut() {
            let ra = new_contact.global.a - old_contact.global.a;
            let rb = new_contact.global.b - old_contact.global.b;
            if ra.magnitude2() <= PERSISTENT_THRESHOLD_SQ ||
                rb.magnitude2() <= PERSISTENT_THRESHOLD_SQ
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
}
