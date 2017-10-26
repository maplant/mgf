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
//! # Physics resolution overview
//!
//! Types that satisfy the `PhysicsObject` trait can resolve collisions in
//! accordance with the laws of classical mechanics. To do this it is
//! recommended to generate a `Manifold` type and pass it to the resolution
//! function. This method handles any geometry and number of contacts:
//!
//! ```rust
//! use mgf::cgmath::{Vector3, Point3, Zero};
//! use mgf::{Sphere, RigidBody, LocalContacts, Manifold, PhysicsState,
//!           DefaultPhysConfig, Component, PhysicsObject, Shape};
//!
//! const RESTITUTION: f32 = 1.0;
//! const FRICTION: f32 = 0.0;
//!
//! let gravity = Vector3::new(0.0, -9.8, 0.0);
//!
//! let body_geom = vec![
//!     Component::from(Sphere{ c: Point3::new(0.0, 0.0, 0.0), r: 1.0 })
//! ];
//! let body_masses = vec![ 1.0 ];
//! 
//! let mut rigid_body_a = RigidBody::new(RESTITUTION, FRICTION, gravity,
//!                                       body_geom, body_masses);
//! // Since constructing a rigid body calculates the center of mass and the
//! // inertia tensor for an object, it is more efficient to clone a rigid
//! // body into a new one that shares the same geometry.
//! let mut rigid_body_b = rigid_body_a.clone();
//!
//! // Set the initial positions of the bodies
//! rigid_body_a.set_pos(Point3::new(0.0, -5.0, 0.0));
//! rigid_body_b.set_pos(Point3::new(0.0,  5.0, 0.0));
//!
//! // Apply an impulse to the bodies to set them in motion.
//! rigid_body_a.apply_impulse(
//!     Vector3::new(5.0, 0.0, 0.0),    // Linear impulse
//!     Vector3::zero()                 // Angular impulse
//! );
//! rigid_body_b.apply_impulse(Vector3::new(-5.0, 0.0, 0.0), Vector3::zero());
//!
//! // Integrate the bodies once per frame to produce the correct motion.
//! rigid_body_a.integrate(1.0);
//! rigid_body_b.integrate(1.0);
//!
//! let mut manifold = Manifold::new();
//!
//! rigid_body_a.local_contacts(&rigid_body_b, |c|{ manifold.push(c); });
//!
//! PhysicsState::resolve_manifold::<DefaultPhysConfig, _, _>(
//!     &mut rigid_body_a,
//!     &mut rigid_body_b,
//!     &manifold,
//!     1.0     // The same timestep passed to integrate.
//! );
//! ```
//!
//! `DefaultPhysConfig` is a struct that passes various constants to the
//! resolver. It does not incur any runtime penalty.

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

