# MGF: Matt's Game Framework

[![Documentation](https://docs.rs/mgf/badge.svg)](https://docs.rs/mgf)
[![Version](https://img.shields.io/crates/v/mgf.svg)](https://crates.io/crates/mgf)
[![License](https://img.shields.io/crates/l/mgf.svg)](https://github.com/DataAnalysisCosby/mgf/blob/master/LICENSE)
[![Downloads](https://img.shields.io/crates/d/mgf.svg)](https://crates.io/crates/mgf)

MGF is a collision detection and physics library for use in 3D video games.

MGF is intended to very light weight and uses cgmath as a math backend.

The library provides various features such as:

- structures to define shapes: `Ray`, `Segment`, `AABB`, `Rectangle`, `Triangle`, `Sphere`, `Capsule`
- structures to define aggregate shapes: `Mesh`, `Compound`
- discrete collision detection: `Overlaps`, `Contains`
- continuous collision detection: `Intersection`, `Contact`, `LocalContact`
- a bounding volume hierarchy: `BVH`
- rigid body physics: `SimpleDynamicBody`, `CompoundDynamicBody`, `StaticBody`, `ContactSolver`
- dynamic containers: `Pool`

MGF is very much in its infancy and is therefore not feature complete. If you
notice any errors or poorly thought out interfaces be sure to let me know.

## 3D only

For the time being MGF is solely designed to handle 3D video games. This 
reflects my own use of MGF. If there is enough demand for MGF to support 2D
games, it may in the future.

## Differences between other collision detection/physics libraries

Most collision detection and physics libraries provide similar to functionality
to MGF, what differs is how they are implemented. Most libraries use the GJK
algorithm to perform every type of collision, including continuous ones. MGF focuses
on providing fast and accurate exact collision detection for moving geometries
commonly found in games, such as `Spheres` and `Capsules`, and does not and cannot
provide fast and accurate moving mesh on mesh collisions.

One interesting result of the implementation is that the direction of the normal
of a polygon determines the "side" it is facing, in a similar manner to how back-face
culling is performed. This allows us to determine inter-object penetration much more
accurately than we would be able to otherwise. The normal force essentially acts
as the direction of propulsion for the object during physics resolution. 

## Examples

You can find working examples in the mgf_demo folder. Demos require gfx and
gfx_glutin to display visuals. Be sure to build the demos in release mode to get
adequate performance:

```
cargo build --release
```

[Here is a link to a video of the demo in action.](https://www.youtube.com/watch?v=bPMm2_ttSq8)

## Contributing

MGF is welcome to anyone's contribution, and any part of the interface is open to 
discussion. Although if you are going to contribute new features please do 
be sure to include unit tests. 
