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

