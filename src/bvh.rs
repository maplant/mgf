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

use std::cmp;
use std::ops::Index;

use smallvec::SmallVec;

use crate::bounds::{Bound, BoundedBy};
use crate::pool::Pool;

use crate::collision::{Intersects, Intersection};

use serde::{Serialize, Deserialize};

/// A Bounding Volume Hierarchy.
#[derive(Serialize, Deserialize)]
pub struct BVH<B: Bound, V> {
    root: usize,
    pool: Pool<BVHNode<B, V>>,
}

#[derive(Serialize, Deserialize)]
struct BVHNode<B: Bound, V> {
    height: i32,
    parent: usize,
    bounds: B,
    node_type: BVHNodeType<V>,
}

#[derive(Serialize, Deserialize)]
enum BVHNodeType<V> {
    Leaf(V),
    Parent(usize, usize),
}

impl<B, V> Clone for BVH<B, V>
where
    B: Bound,
    V: Clone,
{
    fn clone(&self) -> Self {
        BVH {
            root: self.root,
            pool: self.pool.clone(),
        }
    }
}

impl<B, V> Clone for BVHNode<B, V>
where
    B: Bound,
    V: Clone,
{
    fn clone(&self) -> Self {
        BVHNode {
            height: self.height,
            parent: self.parent,
            bounds: self.bounds,
            node_type: self.node_type.clone(),
        }
    }
}

impl<V: Clone> Clone for BVHNodeType<V> {
    fn clone(&self) -> Self {
        match self {
            &BVHNodeType::Leaf(ref v) => BVHNodeType::Leaf(v.clone()),
            &BVHNodeType::Parent(lc, rc) => BVHNodeType::Parent(lc, rc),
        }
    }
}

impl<B: Bound, V> BVH<B, V> {
    /// Creates an empty BVH.
    pub fn new() -> Self {
        BVH {
            root: 0,
            pool: Pool::new(),
        }
    }

    /// Creates a BVH with a preallocated array of cap.
    pub fn with_capacity(cap: usize) -> Self {
        BVH {
            root: 0,
            pool: Pool::with_capacity(cap),
        }
    }

    /// Determines if the BVH is empty.
    pub fn empty(&self) -> bool {
        self.pool.empty()
    }

    /// Removes all entries from the BVH.
    pub fn clear(&mut self) {
        self.root = 0;
        self.pool.clear();
    }

    fn insert_node(&mut self, bounds: B, node_type: BVHNodeType<V>) -> usize {
        self.pool.push(BVHNode {
            bounds,
            node_type,
            height: -1,
            parent: 0,
        })
    }

    /// Inserts an item into the BVH, rebalancing if necessary. All IDs returned
    /// prior to insert remain valid afterward.
    pub fn insert<K: BoundedBy<B>>(&mut self, key: &K, val: V) -> usize {
        let bounds = key.bounds();
        let leaf = self.insert_node(bounds, BVHNodeType::Leaf(val));
        if self.pool.len() == 1 {
            self.root = leaf;
            return leaf;
        }
        let mut best = self.root;
        loop {
            if let BVHNodeType::Parent(child1, child2) = self.pool[best].node_type {
                let curr_bounds = self.pool[best].bounds;
                let area = curr_bounds.surface_area();
                let combined_bounds = B::combine(&curr_bounds, &bounds);
                let combined_area = combined_bounds.surface_area();
                let no_descent_cost = combined_area * 2.0;
                let inheritance_cost = (combined_area - area) * 2.0;

                let child_cost = |child: usize| -> f32 {
                    if let BVHNodeType::Parent(_, _) = self.pool[child].node_type {
                        let old_area = self.pool[child].bounds.surface_area();
                        let new_area = B::combine(&bounds, &self.pool[child].bounds).surface_area();
                        new_area - old_area + inheritance_cost
                    } else {
                        B::combine(&bounds, &self.pool[child].bounds).surface_area() +
                            inheritance_cost
                    }
                };

                let child1_cost = child_cost(child1);
                let child2_cost = child_cost(child2);

                // Descend according to minimum cost
                if no_descent_cost < child1_cost && no_descent_cost < child2_cost {
                    break;
                }

                best = if child1_cost < child2_cost {
                    child1
                } else {
                    child2
                };
            } else {
                break;
            }
        }

        // Create a new parent
        let old_parent = self.pool[best].parent;
        //        let new_parent = self.peek_id();
        let best_bounds = self.pool[best].bounds;
        let new_parent = self.insert_node(
            B::combine(&bounds, &best_bounds),
            BVHNodeType::Parent(best, leaf),
        );
        self.pool[new_parent].parent = old_parent;
        self.pool[new_parent].height = self.pool[best].height + 1;

        if best != self.root {
            if let BVHNodeType::Parent(child1, child2) = self.pool[old_parent].node_type {
                self.pool[old_parent].node_type = if child1 == best {
                    BVHNodeType::Parent(new_parent, child2)
                } else {
                    BVHNodeType::Parent(child1, new_parent)
                };
            }
        } else {
            self.root = new_parent;
        }
        //        self.pool[new_parent].node_type = BVHNodeType::Parent(best, leaf);
        self.pool[best].parent = new_parent;
        self.pool[leaf].parent = new_parent;

        // Walk up the tree fixing the heights and bounds.
        let mut i = self.pool[leaf].parent;
        loop {
            i = self.balance(i);

            // Gauranteed to be a parent.
            if let BVHNodeType::Parent(child1, child2) = self.pool[i].node_type {
                self.pool[i].height =
                    1 + cmp::max(self.pool[child1].height, self.pool[child2].height);
                self.pool[i].bounds =
                    B::combine(&self.pool[child1].bounds, &self.pool[child2].bounds);
                if i == self.root {
                    break;
                };
            }

            i = self.pool[i].parent;
        }

        leaf
    }

    /// Removes a leaf node from the BVH.
    pub fn remove(&mut self, leaf: usize) {
        let parent = self.pool[leaf].parent;
        self.pool.remove(leaf);
        if leaf == self.root {
            self.root = 0;
            return;
        }
        if let BVHNodeType::Parent(child1, child2) = self.pool[parent].node_type {
            let sibling = if child1 == leaf { child2 } else { child1 };
            if self.root != parent {
                let grand_parent = self.pool[parent].parent;
                if let BVHNodeType::Parent(child1, child2) = self.pool[grand_parent].node_type {
                    self.pool[grand_parent].node_type = if child1 == parent {
                        BVHNodeType::Parent(sibling, child2)
                    } else {
                        BVHNodeType::Parent(child1, sibling)
                    }
                }
                self.pool[sibling].parent = grand_parent;
                self.pool.remove(parent);
                let mut i = grand_parent;
                loop {
                    i = self.balance(i);
                    if let BVHNodeType::Parent(child1, child2) = self.pool[i].node_type {
                        self.pool[i].bounds =
                            B::combine(&self.pool[child1].bounds, &self.pool[child2].bounds);
                        self.pool[i].height =
                            1 + cmp::max(self.pool[child1].height, self.pool[child2].height);
                        if self.root == i {
                            break;
                        }

                        i = self.pool[i].parent;
                    }
                }
            } else {
                self.root = sibling;
                self.pool.remove(parent);
            }
        }
    }

    /// Returns the index of the root node.
    pub fn root(&self) -> usize {
        if self.empty() {
            panic!("BVH is empty, there is no root node");
        }
        self.root
    }

    pub fn get_leaf(&self, i: usize) -> &V {
        if let &BVHNodeType::Leaf(ref leaf) = &self.pool[i].node_type {
            leaf
        } else {
            panic!("node at index {} is not a leaf", i);
        }
    }

    /// Finds each entry in the BVH that has a bound that overlaps the bound of
    /// the passed object. Performs a depth first search for all objects.
    ///
    /// callback is called for each entry that overlaps arg. A reference to the
    /// value stored at the leaf is passed.
    pub fn query<Arg, F>(&self, arg: &Arg, mut callback: F)
    where
        Arg: BoundedBy<B>,
        F: FnMut(&V)
    {
        if self.empty() {
            return;
        }
        let arg_bounds = arg.bounds();
        // 64 entries should be enough, since 2^64 - 1 is the maximum number of
        // items a Vec can store.
        let mut stack = SmallVec::<[usize; 64]>::new();
        stack.push(self.root);
        while let Some(top) = stack.pop() {
            if arg_bounds.overlaps(&self.pool[top].bounds) {
                match self.pool[top].node_type {
                    BVHNodeType::Leaf(ref val) => {
                        callback(val);
                    },

                    BVHNodeType::Parent(lchild, rchild) => {
                        stack.push(lchild);
                        stack.push(rchild);
                    }
                }
            }
        }
    }

    /// Finds each entry in the BVH that has a bound that overlaps the bound of
    /// the passed object. Performs a depth first search for all objects.
    ///
    /// callback is called for each entry that overlaps arg. A mutable reference
    /// to the value stored at the leaf is passed.
    pub fn query_mut<Arg, F>(&mut self, arg: &Arg, mut callback: F)
    where
        Arg: BoundedBy<B>,
        F: FnMut(&mut V)
    {
        if self.empty() {
            return;
        }
        let arg_bounds = arg.bounds();
        let mut stack = SmallVec::<[usize; 64]>::new();
        stack.push(self.root);
        while let Some(top) = stack.pop() {
            if arg_bounds.overlaps(&self.pool[top].bounds) {
                match self.pool[top].node_type {
                    BVHNodeType::Leaf(ref mut val) => {
                        callback(val);
                    },

                    BVHNodeType::Parent(lchild, rchild) => {
                        stack.push(lchild);
                        stack.push(rchild);
                    }
                }
            }
        }
    }

    /// Finds all entries that intersect a ray or segment.
    pub fn raytrace<Arg, F>(&self, arg: &Arg, mut callback: F)
    where
        Arg: Intersects<B>,
        F: FnMut(&V, Intersection)
    {
        if self.empty() {
            return;
        }
        let mut stack = SmallVec::<[usize; 64]>::new();
        stack.push(self.root);
        while let Some(top) = stack.pop() {
            if let Some(inter) = arg.intersection(&self.pool[top].bounds) {
                match self.pool[top].node_type {
                    BVHNodeType::Leaf(ref val) => {
                        callback(val, inter);
                    },

                    BVHNodeType::Parent(lchild, rchild) => {
                        stack.push(lchild);
                        stack.push(rchild);
                    }
                }
            }
        }
    }

    fn balance(&mut self, a: usize) -> usize {
        // This could be really cleaned up by using pointers instead of indices
        // everywhere.
        if self.pool[a].height < 2 {
            return a;
        }
        if let BVHNodeType::Parent(b, c) = self.pool[a].node_type {
            if self.pool[c].height > self.pool[b].height + 1 {
                if let BVHNodeType::Parent(f, g) = self.pool[c].node_type {
                    // Swap A and C
                    self.pool[c].parent = self.pool[a].parent;
                    self.pool[a].parent = c;

                    // A's old parent should point to C
                    if self.root == a {
                        self.root = c;
                    } else if let BVHNodeType::Parent(pchild1, pchild2) =
                        self.pool[self.pool[c].parent].node_type
                    {
                        let parent = self.pool[c].parent;
                        self.pool[parent].node_type = if pchild1 == a {
                            BVHNodeType::Parent(c, pchild2)
                        } else {
                            BVHNodeType::Parent(pchild1, c)
                        };
                    }

                    // Rotate and readjust
                    // This really needs to be extracted into its own function
                    if self.pool[f].height > self.pool[g].height {
                        self.pool[c].node_type = BVHNodeType::Parent(a, f);
                        self.pool[a].node_type = BVHNodeType::Parent(b, g);
                        self.pool[g].parent = a;
                        self.pool[a].bounds =
                            B::combine(&self.pool[b].bounds, &self.pool[g].bounds);
                        self.pool[c].bounds =
                            B::combine(&self.pool[a].bounds, &self.pool[f].bounds);
                        self.pool[a].height =
                            1 + cmp::max(self.pool[b].height, self.pool[g].height);
                        self.pool[c].height =
                            1 + cmp::max(self.pool[a].height, self.pool[f].height);
                    } else {
                        self.pool[c].node_type = BVHNodeType::Parent(a, g);
                        self.pool[a].node_type = BVHNodeType::Parent(b, f);
                        self.pool[f].parent = a;
                        self.pool[a].bounds =
                            B::combine(&self.pool[b].bounds, &self.pool[f].bounds);
                        self.pool[c].bounds =
                            B::combine(&self.pool[a].bounds, &self.pool[g].bounds);
                        self.pool[a].height =
                            1 + cmp::max(self.pool[b].height, self.pool[f].height);
                        self.pool[c].height =
                            1 + cmp::max(self.pool[a].height, self.pool[g].height);
                    }
                }
                return c;
            }
            if self.pool[b].height > self.pool[c].height + 1 {
                if let BVHNodeType::Parent(d, e) = self.pool[b].node_type {
                    // Swap A and B
                    self.pool[b].parent = self.pool[a].parent;
                    self.pool[a].parent = b;

                    // A's old parent should point to B
                    if self.root == a {
                        self.root = b;
                    } else if let BVHNodeType::Parent(pchild1, pchild2) =
                        self.pool[self.pool[b].parent].node_type
                    {
                        let parent = self.pool[b].parent;
                        self.pool[parent].node_type = if pchild1 == a {
                            BVHNodeType::Parent(b, pchild2)
                        } else {
                            BVHNodeType::Parent(pchild1, b)
                        };
                    }

                    // Rotate and readjust
                    // This really needs to be extracted into its own function
                    if self.pool[d].height > self.pool[e].height {
                        self.pool[b].node_type = BVHNodeType::Parent(a, d);
                        self.pool[a].node_type = BVHNodeType::Parent(e, c);
                        self.pool[e].parent = a;
                        self.pool[a].bounds =
                            B::combine(&self.pool[c].bounds, &self.pool[e].bounds);
                        self.pool[b].bounds =
                            B::combine(&self.pool[a].bounds, &self.pool[d].bounds);
                        self.pool[a].height =
                            1 + cmp::max(self.pool[c].height, self.pool[e].height);
                        self.pool[b].height =
                            1 + cmp::max(self.pool[a].height, self.pool[d].height);
                    } else {
                        self.pool[b].node_type = BVHNodeType::Parent(a, e);
                        self.pool[a].node_type = BVHNodeType::Parent(d, c);
                        self.pool[d].parent = a;
                        self.pool[a].bounds =
                            B::combine(&self.pool[c].bounds, &self.pool[d].bounds);
                        self.pool[b].bounds =
                            B::combine(&self.pool[a].bounds, &self.pool[e].bounds);
                        self.pool[a].height =
                            1 + cmp::max(self.pool[c].height, self.pool[d].height);
                        self.pool[b].height =
                            1 + cmp::max(self.pool[a].height, self.pool[e].height);
                    }
                }
                return b;
            }
        }
        a
    }
}

impl<B, V> Index<usize> for BVH<B, V>
where
    B: Bound
{
    type Output = B;

    fn index(&self, i: usize) -> &B {
        &self.pool[i].bounds
    }
}

impl<B, V> BoundedBy<B> for BVH<B, V>
where
    B: Bound
{
    fn bounds(&self) -> B {
        if self.empty() {
            panic!("BVH is empty and thus has no bounds");
        }
        self.pool[self.root].bounds
    }
}

#[cfg(test)]
mod tests {
    mod bvh {
        use cgmath::{Point3};
        use crate::geom::{Sphere, AABB};
        use crate::bvh::BVH;

        #[test]
        fn test_bvh() {
            let sphere_a = Sphere{ c: Point3::new(0.0, 5.0, 0.0), r: 1.0 };
            let sphere_b = Sphere{ c: Point3::new(0.0, 8.0, 0.0), r: 1.0 };
            let sphere_c = Sphere{ c: Point3::new(3.0, 0.0, 0.0), r: 1.0 };

            let mut bvh: BVH<AABB, usize> = BVH::new();
            bvh.insert(&sphere_a, 1);
            bvh.insert(&sphere_b, 2);
            bvh.insert(&sphere_c, 3);

            let mut found: usize = 0;
            bvh.query(&sphere_a, |&id| { found += 1; assert_eq!(id, 1); });
            bvh.query(&sphere_b, |&id| { found += 1; assert_eq!(id, 2); });
            bvh.query(&sphere_c, |&id| { found += 1; assert_eq!(id, 3); });
            assert_eq!(found, 3);
        }
    }
}
