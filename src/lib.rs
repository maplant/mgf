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
//! - `SimpleDynamicBody`, a dynamic rigid body with a single geometry as a volume.
//! - `CompoundDynamicBody`, a dynamic rigid body with a `Compound` geometry.
//! - `StaticBody`, a static body with infinite mass (i.e. terrain) that can be any shape.
//!
//! After contact points are accumulated, they are pruned via a `ContactPruner` and put
//! through a constraint solver (called `ContactSolver`). The constraint solver is a
//! slightly unsafe interface that requires passing a mutable pointer for each object in
//! a contact. It will operate fine if you unsafely pass the same mutable pointer multiple
//! times for different contacts. In fact, that's very typical. In the future this interface
//! may be made safer but for now, be safe. Here is an example of a typical use case:
//!
//! ```rust
//! use mgf::cgmath::prelude::*;
//! use mgf::cgmath::*;
//! use mgf::{ContactSolver, ContactPruner, LocalContacts, Manifold, 
//!           PhysicsObject, SimpleDynamicBody, Shape, Sphere, Velocity};
//!
//! const TIMESTEP: f32 = 1.0;
//! const RESTITUTION: f32 = 0.3;
//! const FRICTION: f32 = 0.5;
//!
//! let gravity = Vector3::new(0.0, -9.8, 0.0);
//! let mass = 1.0;
//! let mut body_a = SimpleDynamicBody::new(
//!     RESTITUTION, FRICTION, gravity,
//!     Sphere {
//!         c: Point3::new(0.0, 0.0, 0.0),
//!         r: 1.0,
//!     },
//!     mass
//! );
//! // It is generally faster to clone a physics object than it is to construct
//! // a new one.
//! let mut body_b = body_a.clone();
//! 
//! body_a.set_pos(Point3::new(-5.0, 0.0, 0.0));
//! body_b.set_pos(Point3::new( 5.0, 0.0, 0.0));
//!
//! // Set the linear velocity of each of the objects.
//! body_a.set_dx(Velocity{ linear: Vector3::new(0.0, 4.0, 0.0), angular: Vector3::zero() });
//! body_b.set_dx(Velocity{ linear: Vector3::new(0.0, -4.0, 0.0), angular: Vector3::zero() });
//!
//! // ball_a and ball_b are set to collide at the end of one unit of time.
//! // Integrate them to set their motion and hitboxes.
//! body_a.integrate(TIMESTEP);
//! body_b.integrate(TIMESTEP);
//!
//! // Collide the objects:
//! let mut pruner: ContactPruner = ContactPruner::new();
//! body_a.local_contacts(&body_b, |contact|{ pruner.push(contact); });
//!
//! // Put the contacts in the constraint solver.
//! // This is where the interface gets hairy, as adding a constraint to the
//! // solver requires two mutable references. Often times physics objects will
//! // be stored in a Vec or a Pool, so obtaining mutable references for two
//! // objects is difficult. Unsafe is often necessary.
//! const SOLVER_ITERS: usize = 20;
//! let mut solver: ContactSolver = ContactSolver::new();
//! solver.add_constraint(&mut body_a, &mut body_b, Manifold::from(pruner), TIMESTEP);
//! solver.solve(SOLVER_ITERS);
//! ```

#[macro_use]
pub extern crate cgmath;
extern crate smallvec;

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

