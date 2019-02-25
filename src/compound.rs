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
use std::vec::Vec;
use std::ops::{Add, AddAssign, Sub, SubAssign};
use cgmath::prelude::*;
use cgmath::{EuclideanSpace, Rotation, Rotation3, Vector3, Point3, Quaternion, One, Zero};

use smallvec::SmallVec;

use crate::bvh::*;
use crate::bounds::*;
use crate::collision::*;
use crate::geom::*;

/// A component is a generic volume that can either be a Sphere or Capsule at
/// runtime. Anything that can collide with a Sphere and a Capsule can collide
/// with a component. 
#[derive(Copy, Clone, Debug)]
pub enum Component {
    Sphere(Sphere),
    Capsule(Capsule),
    // More shapes to come...
}

impl Component {
    /// Create a ComponentConstructor and return the position and rotation of
    /// the Component.
    pub fn deconstruct(&self) -> (Point3<f32>, Quaternion<f32>, ComponentConstructor)  {
        match self {
            &Component::Sphere(Sphere{ r, c }) =>
                (c, Quaternion::one(), ComponentConstructor::Sphere{ r }),
            &Component::Capsule(Capsule{ r, a, d }) => {
                let h = d.magnitude();
                let rot = Quaternion::from_arc(Vector3::new(0.0, 1.0, 0.0) * h, d, None);
                (a + d * 0.5, rot, ComponentConstructor::Capsule{ r, half_h: h * 0.5 })
            },
        }
    }
}

impl Volumetric for Component {
    fn rotate<R: Rotation3<f32>>(self, r: R) -> Self {
        match self {
            Component::Sphere(s) => Component::Sphere(s.rotate(r)),
            Component::Capsule(c) => Component::Capsule(c.rotate(r)),
        }
    }
}

impl From<Sphere> for Component {
    fn from(s: Sphere) -> Self {
        Component::Sphere(s)
    }
}

impl From<Capsule> for Component {
    fn from(c: Capsule) -> Self {
        Component::Capsule(c)
    }
}

impl Add<Vector3<f32>> for Component {
    type Output = Self;

    fn add(self, v: Vector3<f32>) -> Self {
        match self {
            Component::Sphere(s) => Component::Sphere(s + v),
            Component::Capsule(c) => Component::Capsule(c + v),
        }
    }
}
        
impl Sub<Vector3<f32>> for Component {
    type Output = Self;

    fn sub(self, v: Vector3<f32>) -> Self {
        match self {
            Component::Sphere(s) => Component::Sphere(s - v),
            Component::Capsule(c) => Component::Capsule(c - v),
        }
    }
}

impl AddAssign<Vector3<f32>> for Component {
    fn add_assign(&mut self, v: Vector3<f32>) {
        match self {
            &mut Component::Sphere(ref mut s) => *s += v,
            &mut Component::Capsule(ref mut c) => *c += v,
        }
    }
}

impl SubAssign<Vector3<f32>> for Component {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        match self {
            &mut Component::Sphere(ref mut s) => *s -= v,
            &mut Component::Capsule(ref mut c) => *c -= v,
        }
    }
}

impl Shape for Component {
    fn center(&self) -> Point3<f32> {
        match self {
            &Component::Sphere(s) => s.center(),
            &Component::Capsule(c) => c.center(),
        }
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        match self {
            &Component::Sphere(s) => s.closest_point(to),
            &Component::Capsule(c) => c.closest_point(to),
        }
    }
}

impl BoundedBy<AABB> for Component {
    fn bounds(&self) -> AABB {
        match self {
            &Component::Sphere(s) => s.bounds(),
            &Component::Capsule(c) => c.bounds(),
        }
    }
}

impl BoundedBy<Sphere> for Component {
    fn bounds(&self) -> Sphere {
       match self {
            &Component::Sphere(s) => s.bounds(),
            &Component::Capsule(c) => c.bounds(),
        }
    }
}

impl<P: Particle> Intersects<Component> for P {
    fn intersection(&self, rhs: &Component) -> Option<Intersection> {
        match rhs {
            &Component::Sphere(ref s) => self.intersection(s),
            &Component::Capsule(ref c) => self.intersection(c),
        }
    }
}

macro_rules! impl_component_collision {
    (
        $recv:ty
    ) => {
        impl Contacts<Moving<Component>> for $recv {
            fn contacts<F: FnMut(Contact)>(&self, rhs: &Moving<Component>, callback: F) -> bool {
                match rhs.0 {
                    Component::Sphere(s) => self.contacts(&Moving::sweep(s, rhs.1), callback),
                    Component::Capsule(c) => self.contacts(&Moving::sweep(c, rhs.1), callback),
                }
            }
        }
    };
}

impl_component_collision!{ Plane }
impl_component_collision!{ Triangle }
impl_component_collision!{ Rectangle }
impl_component_collision!{ Sphere }
impl_component_collision!{ Capsule }

impl<RHS> Contacts<RHS> for Moving<Component>
where
    RHS: Contacts<Moving<Sphere>> + Contacts<Moving<Capsule>>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        match self.0 {
            Component::Sphere(s) => rhs.contacts(&Moving::sweep(s, self.1), |c|callback(-c)),
            Component::Capsule(c) => rhs.contacts(&Moving::sweep(c, self.1), |c|callback(-c)),
        }
    }
}

impl LocalContacts<Moving<Component>> for Moving<Component> {
    fn local_contacts<F: FnMut(LocalContact)>(&self, rhs: &Moving<Component>, mut callback: F) -> bool {
        self.contacts(
            rhs,
            | c | {
                callback(
                    LocalContact {
                        local_a: c.a + -(self.0.center() + self.1 * c.t).to_vec(),
                        local_b: c.b + -(rhs.0.center() + rhs.1 * c.t).to_vec(),
                        global: c,
                    }
                );
            }
        )
    }
}

/// A description of a Component minus rotation and position.
#[derive(Copy, Clone, Debug)]
pub enum ComponentConstructor {
    Sphere{ r: f32 },
    Capsule{ r: f32, half_h: f32 },
    // More shapes to come...
}

impl ComponentConstructor {
    /// Create a component from a component constructor.
    pub fn construct<R: Rotation3<f32>>(&self, p: Point3<f32>, rot: R) -> Component {
        match self {
            &ComponentConstructor::Sphere{ r } => Component::Sphere(Sphere{ r, c: p }),
            &ComponentConstructor::Capsule{ r, half_h } => {
                let d = rot.rotate_vector(Vector3::new(0.0, 1.0, 0.0) * half_h); 
                Component::Capsule(Capsule{ r: r, a: p + -d, d: d * 2.0 })
            },
        }
    }
}

/// An aggregate structure of Spheres and Capsules. Has a position and rotation.
#[derive(Clone)]
pub struct Compound {
    /// The displacement of the object.
    pub disp: Vector3<f32>,
    /// The rotation of the object. Assumed to be normalized.
    pub rot: Quaternion<f32>,
    /// Indices of the geometries composing the compound in the BVH.
    /// One-to-one with the constructing vector.
    pub shapes: SmallVec<[usize; 1]>,
    /// BVH storing the components to improve collision efficiency.
    pub bvh: BVH<AABB, Component>,
}

impl Compound {
    pub fn new(components: Vec<Component>) -> Self {
        let mut bvh: BVH<AABB, Component> = BVH::new();
        let mut shapes: SmallVec<[usize; 1]> = SmallVec::with_capacity(components.len());
        for component in components.iter() {
            shapes.push(bvh.insert(component, *component));
        }
        Compound {
            disp: Vector3::zero(),
            rot: Quaternion::one(),
            shapes,
            bvh,
        }
    }
}

impl AddAssign<Vector3<f32>> for Compound {
    fn add_assign(&mut self, v: Vector3<f32>) {
        self.disp += v
    }
}

impl SubAssign<Vector3<f32>> for Compound {
    fn sub_assign(&mut self, v: Vector3<f32>) {
        self.disp -= v
    }
}

impl BoundedBy<AABB> for Compound {
    fn bounds(&self) -> AABB {
        self.bvh[self.bvh.root()].rotate(self.rot) + self.disp
    }
}

impl BoundedBy<Sphere> for Compound {
    fn bounds(&self) -> Sphere {
        let s: Sphere = self.bvh[self.bvh.root()].bounds();
        s + self.disp
    }
}

impl Shape for Compound {
    /// The point returned by a compound shape is the displacement of the
    /// object and not the center of mass. It is impossible for Compound
    /// to calculate the center of mass given it has no information regarding
    /// the mass of individual components.
    fn center(&self) -> Point3<f32> {
        Point3::from_vec(self.disp)
    }

    fn closest_point(&self, to: Point3<f32>) -> Point3<f32> {
        let mut best_p = Point3::new(0.0, 0.0, 0.0);
        let mut best_dist: f32 = f32::INFINITY;
        for shape in self.shapes.iter() {
            let new_p = self.bvh.get_leaf(*shape).closest_point(to);
            let new_dist = (to - new_p).magnitude2();
            if new_dist < best_dist {
                best_p = new_p;
                best_dist = new_dist;
            }
        }
        best_p
    }
}

impl<P: Particle> Intersects<Compound> for P {
    fn intersection(&self, rhs: &Compound) -> Option<Intersection> {
        let conj_rot = rhs.rot.conjugate();
        let p = conj_rot.rotate_point(self.pos() + -rhs.disp) + rhs.disp;
        let d = conj_rot.rotate_vector(self.dir());
        let r = Ray{ p, d };
        let mut result: Option<Intersection> = None;
        rhs.bvh.raytrace(&r, |&comp, inter| {
            if inter.t > P::DT {
                return;
            }
            let shape = comp.rotate(rhs.rot) + rhs.disp;
            if let Some(inter) = self.intersection(&shape) {
                if let Some(res) = result {
                    if inter.t > res.t {
                        return;
                    }
                }
                result = Some(inter)
            }
        });
        result
    }
}

impl<RHS> Contacts<RHS> for Compound
where
    RHS: Contacts<Component> + BoundedBy<AABB>
{
    fn contacts<F: FnMut(Contact)>(&self, rhs: &RHS, mut callback: F) -> bool {
        // We assume that self.rot is normalized.
        let conj_rot = self.rot.conjugate();
        let mut rhs_bounds = rhs.bounds().rotate(conj_rot);
        let rhs_center = rhs_bounds.center();
        let bounds_disp = conj_rot.rotate_point(rhs_center + -self.disp) + self.disp;
        rhs_bounds.set_pos(bounds_disp);
        let mut collided = false;
        self.bvh.query(&rhs_bounds, |&comp| {
            let shape = comp.rotate_about(self.rot, Point3::new(0.0, 0.0, 0.0)) + self.disp;
            rhs.contacts(&shape, |c| { collided = true; callback(-c) });
        });
        collided
    }
}

#[cfg(test)]
mod tests {
    mod compound {
        use cgmath::InnerSpace;
        use crate::compound::*;
        use crate::collision::Contacts;

        #[test]
        fn test_compound() {
            let components = vec![
                Component::Sphere(Sphere{ c: Point3::new(-5.0, 0.0, 0.0), r: 1.0 }),
                Component::Sphere(Sphere{ c: Point3::new(5.0, 0.0, 0.0), r: 1.0 }),
            ];
            let mut compound = Compound::new(components);
            let test_sphere = Moving::sweep(Sphere{ c: Point3::new(0.0, 8.0, 0.0), r: 1.0 },
                                            Vector3::new(0.0, -1.5, 0.0));
            assert!(!compound.contacts(&test_sphere, |c: Contact| { panic!("c = {:?}", c); }));
            // rotate compounds
            compound.rot = Quaternion::from_arc(Vector3::new(1.0, 0.0, 0.0),
                                                Vector3::new(0.0, 1.0, 0.0),
                                                None).normalize();
            let contact: Contact = compound.last_contact(&test_sphere).unwrap();
            assert_relative_eq!(contact.t, 0.6666663, epsilon = COLLISION_EPSILON);
            assert_relative_eq!(contact.a, Point3::new(0.0, 6.0, 0.0), epsilon = COLLISION_EPSILON);

            let static_rect = Rect {
                c: Point3::new(0.0, -2.0, 0.0),
                u: [ Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0) ],
                e: [ 6.0, 6.0 ],
            };

            compound.rot = Quaternion::one();

            let _contact: Contact = compound.last_contact(&Moving::sweep(static_rect, Vector3::new(0.0, 3.0, 0.0))).unwrap();
        }
    }
}
