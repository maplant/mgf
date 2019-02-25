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

use std::fmt::Binary;

/// A bit set of fixed, limited capacity.
pub trait FixedSizeBitSet : Default + Binary {
    /// Maximum number of bits allowed in the set.
    const NUM_BITS: usize;

    /// Retrieve a bit value at the index.
    fn get(&self, i: usize) -> bool;

    /// Set the bit at the given index.
    fn insert(&mut self, i: usize);

    /// Unset the bit at the given index.
    fn remove(&mut self, i: usize);
}

macro_rules! impl_bit_set {
    (
        $type:ty, $num_bits:expr
    ) => {
        impl FixedSizeBitSet for $type {
            const NUM_BITS: usize = $num_bits;

            fn get(&self, i: usize) -> bool {
                if i >= Self::NUM_BITS {
                    panic!("index is out of bounds: the len is {} but the index is {}",
                           Self::NUM_BITS, i);
                }
                (*self >> i & 0b_1) == 1
            }

            fn insert(&mut self, i: usize) {
                if i >= Self::NUM_BITS {
                    panic!("index is out of bounds: the len is {} but the index is {}",
                           Self::NUM_BITS, i);
                }
                *self |= 1 << i;
            }

            fn remove(&mut self, i: usize) {
                if i >= Self::NUM_BITS {
                    panic!("index is out of bounds: the len is {} but the index is {}",
                           Self::NUM_BITS, i);
                }
                *self &= !(1 << i);
            }
        }
    };
}

impl_bit_set!(u8, 8);
impl_bit_set!(u16, 16);
impl_bit_set!(u32, 32);
impl_bit_set!(u64, 64);

#[cfg(test)]
mod tests {
    mod bitset {
        #[test]
        fn test_bitset() {
            use crate::bitset::FixedSizeBitSet;

            let mut bitset: u64 = 0;

            bitset.insert(3);
            bitset.insert(10);
            bitset.insert(13);
            bitset.insert(35);

            for i in 0..32 {
                match i {
                    3 | 10 | 13 | 35 => { assert!(bitset.get(i)); },
                    _ => { assert!(!bitset.get(i)); },
                }
            }
        }
    }
}
