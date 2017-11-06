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
- rigid body physics: `SimpleDynamicBody`, `CompoundDynamicBody`, `StaticBody`
- dynamic containers: `Pool`

MGF is very much in its infancy and is therefore not feature complete. If you
notice any errors or poorly thought out interfaces be sure to let me know.

## 3D only

For the time being MGF is solely designed to handle 3D video games. This 
reflects my own use of MGF. If there is enough demand for MGF to support 2D
games, it may in the future.

## Examples

You can find working examples in the mgf_demo folder. Demos require gfx and
gfx_glutin to display visuals. Be sure to build the demos in release mode to get
adequate performance.

```
cargo build --release
```

to build 

## Contributing

MGF is welcome to anyone's contribution, and any part of the interface is open to 
discussion. Although if you are going to contribute new features please do 
be sure to include unit tests. 
