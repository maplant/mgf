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

//! A low-level 3D collision and physics library intended for use in 3D video
//! game development.
//!
//! # Collision detection overview
//!
//! An object implements collision detection routines through the various traits
//! provided by MGF. Each trait provides a disparate form of collision detection:
//!
//! - `Overlaps`: Specifies a collision where the only information required is
//!   whether or not two objects overlap. 
//! - `Contains`: Specifies a collision where one object completely contains
//!   another. For bounding volumes.
//! - `Intersection`: For linear geometries, describes a collision at some time
//!   during the motion of a particle (represented with a Ray or Segment) with
//!   only one point of contact.
//! - `Contacts`: For moving volumetric geometries, describes a collision with
//!   possibly multiple contacts occuring at some time during the motion of the
//!   geometries.
//! - `LocalContacts`: Similar to `Contacts`, but returns contacts with local
//!    contact points as well as globals.
//!
//! In addition to these traits, MGF supports broad-phase colision detection
//! through a bounding volume hierarch (`BVH`)
//!
//! # Physics Overview
//!
//! Physics is implemented through the following objects:
//!
//! - `RigidBodyVec`, a vector of rigid bodies.
//! - `Solver`, a generic constraint-based solver.
//! - `ContactConstraint`, a constraint modeling the collision of two objects.
//!
//! After contact points are accumulated, they are pruned via a `ContactPruner` and put
//! through a constraint solver.
//!
//! ```rust
//! use mgf::cgmath::prelude::*;
//! use mgf::cgmath::*;
//! use mgf::*;
//!
//! const TIMESTEP: f32 = 1.0;
//! const RESTITUTION: f32 = 0.3;
//! const FRICTION: f32 = 0.5;
//! const MASS: f32 = 1.0;
//! const SOLVER_ITERS: usize = 20;
//!
//! let gravity = Vector3::new(0.0, -9.8, 0.0);
//!
//! let mut rigid_bodies = RigidBodyVec::new();
//! let mut sphere = Component::from(Sphere{ c: Point3::new(0.0, 0.0, 0.0), r: 1.0 });
//!
//! sphere.set_pos(Point3::new(-5.0, 0.0, 0.0));
//! let body_a = rigid_bodies.add_body(sphere, MASS, RESTITUTION, FRICTION, gravity);
//!
//! sphere.set_pos(Point3::new(5.0, 0.0, 0.0));
//! let body_b = rigid_bodies.add_body(sphere, MASS, RESTITUTION, FRICTION, gravity);
//!
//! // Set the velocity for the objects.
//! rigid_bodies.set(body_a, Velocity{ linear: Vector3::new(0.0,  4.0, 0.0), angular: Vector3::zero() });
//! rigid_bodies.set(body_b, Velocity{ linear: Vector3::new(0.0, -4.0, 0.0), angular: Vector3::zero() });
//!
//! // Integrate the objects, creating colliders for each.
//! rigid_bodies.integrate(TIMESTEP);
//!
//! // Collide the objects:
//! let mut pruner: ContactPruner = ContactPruner::new();
//! let body_ai: usize = body_a.into();
//! let body_bi: usize = body_b.into();
//! rigid_bodies.collider[body_ai].local_contacts(
//!     &rigid_bodies.collider[body_bi],
//!     | lc | {
//!         pruner.push(lc);
//!     }
//! );
//!
//! // Create the constraint solver:
//! let mut solver = Solver::<ContactConstraint<_>>::new();
//!
//! // Create a manifold to pass to the contact solver as a constraint:
//!
//! // Obviously two spheres will only produce one contact. To avoid using a
//! // pruner in this case a Manifold may be constructed from a single LocalContact:
//! // let manifold = Manifold::from(local_contact);
//!
//! let manifold = Manifold::from(pruner);
//! solver.add_constraint(
//!     ContactConstraint::new(
//!         &rigid_bodies,
//!         body_a, body_b,
//!         manifold,
//!         TIMESTEP
//!     )
//! );
//!
//! // Solve the collision:
//! solver.solve(&mut rigid_bodies, SOLVER_ITERS);
//! ```
#[macro_use]
pub extern crate cgmath;
extern crate smallvec;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;

mod bvh;
pub use bvh::*;

mod bounds;
pub use bounds::*;

pub mod bitset;

mod compound;
pub use compound::*;

mod collision;
pub use collision::*;

mod solver;
pub use solver::*;

mod geom;
pub use geom::*;

mod manifold;
pub use manifold::*;

mod mesh;
pub use mesh::*;

mod physics;
pub use physics::*;

mod pool;
pub use pool::*;

mod simplex;
pub use simplex::*;
