# Change Log

All notable changes to this project will be documented in this file, following
the format defined at [keepachangelog.com](http://keepachangelog.com/).
This project adheres to [Semantic Versioning](http://semver.org/).

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
