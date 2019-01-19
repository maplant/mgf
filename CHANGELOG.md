# Change Log

All notable changes to this project will be documented in this file, following
the format defined at [keepachangelog.com](http://keepachangelog.com/).
This project adheres to [Semantic Versioning](http://semver.org/). 

## [v1.2.2] - 2019-01-18

- Implemented `Convex` for `Capsule`. I have no idea if this is correct at all,
  hopefully will add some unit tests soon.
- Update cgmath version to 0.17.

## [v1.2.0] - 2018-09-15

Discrete collision detection has been added, although at the moment it is 
limited to convex object collisions. That means stationary objects that only 
implement convex cannot collide with certain primitives, such as triangles and 
planes.

- Added a type parameter to `Convex` in order to allow for points that include extra information.
- Added `SupportPoint` for storing local support points during GJK.
- Added `rotate_around` method to trait `Particle`.
- Added `OBB` object to represent abitrarily oriented bounding boxes.
- Expanded `Simplex` to take a type parameter indicating the desired point type.
- Added `new` support functions to geometries. 
- Fixed a major bug in the `Triangle` implementation of `closest_point`.
- Added convenience methods to `Pool`.

## [v1.1.0] - 2018-04-02

A lot of work went into adding support for discrete collision detection. It is 
not feature complete and not very well tested but the basic algorithm works fine.

For now, only separation can be queried (and thus if two objects are penetrating). 
In the future contact generation will be added.

- Added `Penetrates` trait to query the separation between two shapes.
- Added `Convex` trait to query the support function of a shape.
- Added `MinkowskiDiff` object to construct the Minkowski difference of two shapes.
- Added `Simplex` struct to facilitate in GJK

## [v1.0.0] - 2017-11-20

Although there have been breaking changes between minor versions in the past,
not having a clear definition of how to use version numbers was a detriment. 
Thus, because there are breakind changes between this version and the last, 
version 1.0.0 marks (among other things) the first time we took semantic 
versioning seriously.

### Added

- Added `RigidBodyVec` and associated `RigidBodyRef` as a faster way to contain and update rigid bodies.
- Added `RigidBodyInfo` as a more data-driven alternative to `PhysicsObject`
- Added `ConstrainedSet` trait as a way to look up and set constraint variables safely.
- Added `Constraint` trait as a generic constraint type.
- Added `Solver` as a generic constraint solver.
- Added `ComponentConstructor` as a way to dynamically create `Components` at runtime. 

### Changed

- `ContactConstraint` has been modified to fit the greater constraint solver API changes.
- `Mesh` now satisfies `Volumetric` and thus can be rotated.
- `Mesh` no longer contains per-vertex normal because it was misleading for it
  to contain information it didn't use in collisions. 
- `Volumetric` trait methods now take a borrowed self instead of a self reference to allow 
  for `Mesh` to propely implement rotation.
- Renamed `ContactSolverParams` to `ContactConstraintParams`.
- Renamed `DefaultContactSolverParams` to `DefaultContactConstraintParams`.

### Removed 

- `SimpleDynamicBody`
- `CompoundDynamicBody`
- `StaticBody`
- `PhysicsObject`
- `Vertex`
- `ContactSolver`
