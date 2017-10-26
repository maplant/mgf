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

use std::vec::Vec;
use std::ops::{Add, AddAssign, Sub, SubAssign};
use cgmath::{EuclideanSpace, Rotation, Vector3, Point3, Quaternion, One, Zero};

use bvh::*;
use bounds::*;
use geom::*;

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
    /// Rotates the component abound the origin.
    pub fn rotate(self, r: Quaternion<f32>) -> Self {
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
/*
macro_rules! impl_component_collision {
    (
        $recv:ty
    ) => {
        impl Collider<Contact, Component> for Moving<$recv> {
            fn collide<F: FnMut(Contact)>(&self, rhs: &Component, callback: F) -> bool {
                match rhs {
                    &Component::Sphere(ref s) => self.collide(s, callback),
                    &Component::Capsule(ref c) => self.collide(c, callback),
                }
            }
        }
    };
}

impl_component_collision!(Plane);
impl_component_collision!(Triangle);
impl_component_collision!(Rectangle);
impl_component_collision!(Sphere);
impl_component_collision!(Capsule);
*/
impl Collider<Contact, Moving<Component>> for Rectangle {
    fn collide<F: FnMut(Contact)>(&self, rhs: &Moving<Component>, callback: F) -> bool {
        match rhs.0 {
            Component::Sphere(s) => self.collide(&Moving::sweep(s, rhs.1), callback),
            Component::Capsule(c) => self.collide(&Moving::sweep(c, rhs.1), callback),
        }
    }
}

/*
impl<T> Collider<Contact, T> for Component
where
    T: Collider<Contact, Sphere> + Collider<Contact, Capsule>
{
    fn collide<F: FnMut(Contact)>(&self, rhs: &T, mut callback: F) -> bool {
        match self {
            &Component::Sphere(ref s) => rhs.collide(s, |c| callback(-c)),
            &Component::Capsule(ref c) => rhs.collide(c, |c| callback(-c)),
        }
    }
}
 */   

/*
impl<T> Collider<Contact, Moving<T>> for Component
where
    T: Collider<Contact, Moving<Sphere>> + Collider<Contact, Moving<Capsule>> + Shape + Copy
{
    fn collide<F: FnMut(Contact)>(&self, rhs: &Moving<T>, mut callback: F) -> bool {
        match self {
            &Component::Sphere(s) => {
                rhs.as_ref().collide(&Moving::sweep(s, -rhs.1), |c: Contact| callback(-c))
            },
            &Component::Capsule(c) => {
                rhs.as_ref().collide(&Moving::sweep(c, -rhs.1), |c: Contact| callback(-c))
            },
        }
    }
}
*/

impl<T> Collider<Contact, T> for Moving<Component>
where
    T: Collider<Contact, Moving<Sphere>> + Collider<Contact, Moving<Capsule>> 
{
    fn collide<F: FnMut(Contact)>(&self, rhs: &T, mut callback: F) -> bool {
        match self.0 {
            Component::Sphere(s) => rhs.collide(&Moving::sweep(s, self.1), |c|callback(-c)),
            Component::Capsule(c) => rhs.collide(&Moving::sweep(c, self.1), |c|callback(-c)),
        }
    }
}

impl Collider<Contact, Component> for Moving<Sphere> {
    fn collide<F: FnMut(Contact)>(&self, rhs: &Component, callback: F) -> bool {
        match rhs {
            &Component::Sphere(ref s) => self.collide(s, callback),
            &Component::Capsule(ref c) => self.collide(c, callback),
        }
    }
}

impl Collider<Contact, Component> for Moving<Capsule> {
    fn collide<F: FnMut(Contact)>(&self, rhs: &Component, callback: F) -> bool {
        match rhs {
            &Component::Sphere(ref s) => self.collide(s, callback),
            &Component::Capsule(ref c) => self.collide(c, callback),
        }
    }
}

#[derive(Clone)]
pub struct Compound {
    pub disp: Vector3<f32>,
    /// Rotation, assumed to be normalized.
    pub rot: Quaternion<f32>,
    pub shapes: Vec<Component>,
    pub bvh: BVH<AABB, usize>,
}

impl Compound {
    pub fn new(shapes: Vec<Component>) -> Self {
        let mut bvh: BVH<AABB, usize> = BVH::new();
        for (i, component) in shapes.iter().enumerate() {
            bvh.insert(component, i);
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
}

impl<T> Collider<Contact, T> for Compound
where
    T: Collider<Contact, Component> + BoundedBy<AABB>
{
    fn collide<F: FnMut(Contact)>(&self, rhs: &T, mut callback: F) -> bool {
        // We assume that self.rot is normalized.
        let conj_rot = self.rot.conjugate();
        let mut rhs_bounds = rhs.bounds().rotate(conj_rot);
        let rhs_center = rhs_bounds.center();
        let bounds_disp = conj_rot.rotate_point(rhs_center + -self.disp) + self.disp;
        rhs_bounds.set_pos(bounds_disp);
        let mut collided = false;
        self.bvh.collide(&rhs_bounds, |comp_i| {
            let shape = self.shapes[comp_i].rotate(self.rot) + self.disp;
            rhs.collide(&shape, |c| { collided = true; callback(-c) });
        });
        collided
    }
}

#[cfg(test)]
mod tests {
    mod compound {
        use cgmath::InnerSpace;
        use compound::*;

        #[test]
        fn test_compound() {
            let components = vec![
                Component::Sphere(Sphere{ c: Point3::new(-5.0, 0.0, 0.0), r: 1.0 }),
                Component::Sphere(Sphere{ c: Point3::new(5.0, 0.0, 0.0), r: 1.0 }),
            ];
            let mut compound = Compound::new(components);
            let test_sphere = Moving::sweep(Sphere{ c: Point3::new(0.0, 8.0, 0.0), r: 1.0 },
                                            Vector3::new(0.0, -1.5, 0.0));
            assert!(!compound.collide(&test_sphere, |c: Contact| { panic!("c = {:?}", c); }));
            // rotate compounds
            compound.rot = Quaternion::from_arc(Vector3::new(1.0, 0.0, 0.0),
                                                Vector3::new(0.0, 1.0, 0.0),
                                                None).normalize();
            let contact: Contact = compound.check_collision(&test_sphere).unwrap();
            assert_relative_eq!(contact.t, 0.6666663, epsilon = COLLISION_EPSILON);
            assert_relative_eq!(contact.a, Point3::new(0.0, 6.0, 0.0), epsilon = COLLISION_EPSILON);

            let static_rect = Rect {
                c: Point3::new(0.0, -2.0, 0.0),
                u: [ Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0) ],
                e: [ 6.0, 6.0 ],
            };

            let _contact: Contact = static_rect.check_collision(&Moving::sweep(compound.shapes[0], Vector3::new(0.0, -3.0, 0.0))).unwrap();
        }
    }
}
