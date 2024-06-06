use crate::error::Error;
use crate::input::Input;
use crate::serde;

mod collapsing_lowest;

use crate::index_mapping::IndexMapping;
use crate::output::Output;
use crate::sketch::{Flag, FlagType};
pub use collapsing_lowest::CollapsingLowestDenseStore;

pub trait Store: Send + Sync {
    fn add(&mut self, index: i32, count: f64);
    fn add_bin(&mut self, bin: (i32, f64));
    fn merge_with(&mut self, bins: Vec<(i32, f64)>) {
        for bin in bins {
            self.add_bin(bin)
        }
    }
    fn clear(&mut self);
    fn is_empty(&self) -> bool;
    fn get_total_count(&self) -> f64;
    fn get_offset(&self) -> i32;
    fn get_min_index(&self) -> i32;
    fn get_max_index(&self) -> i32;
    fn get_count(&self, i: i32) -> f64;
    fn encode(&self, output: &mut Output, store_flag_type: FlagType) -> Result<(), Error> {
        if self.is_empty() {
            return Ok(());
        }

        let min_index = self.get_min_index();
        let max_index = self.get_max_index();
        let offset = self.get_offset();

        let mut dense_encoding_size: i64 = 0;
        let num_bins: i64 = max_index as i64 - min_index as i64 + 1;
        dense_encoding_size += serde::unsigned_var_long_encoded_length(num_bins);
        dense_encoding_size += serde::signed_var_long_encoded_length(min_index as i64);
        dense_encoding_size += serde::signed_var_long_encoded_length(1);

        let mut sparse_encoding_size: i64 = 0;
        let mut num_non_empty_bins: i64 = 0;
        let mut previous_index: i64 = 0;

        for i in min_index - offset..max_index - offset + 1 {
            let count = self.get_count(i);
            let count_var_double_encoded_length = serde::var_double_encoded_length(count);
            dense_encoding_size += count_var_double_encoded_length;
            if count != 0.0 {
                num_non_empty_bins += 1;
                let index: i64 = offset as i64 + i as i64;
                sparse_encoding_size +=
                    serde::signed_var_long_encoded_length(index - previous_index);
                sparse_encoding_size += count_var_double_encoded_length;
                previous_index = index;
            }
        }

        if dense_encoding_size <= sparse_encoding_size {
            BinEncodingMode::ContiguousCounts
                .to_flag(store_flag_type)
                .encode(output)?;
            serde::encode_unsigned_var_long(output, num_bins)?;
            serde::encode_signed_var_long(output, min_index as i64)?;
            serde::encode_signed_var_long(output, 1)?;
            for i in min_index - offset..max_index - offset + 1 {
                serde::encode_var_double(output, self.get_count(i))?;
            }
        } else {
            BinEncodingMode::IndexDeltasAndCounts
                .to_flag(store_flag_type)
                .encode(output)?;
            serde::encode_unsigned_var_long(output, num_non_empty_bins)?;
            let mut previous_index = 0;
            for i in min_index - offset..max_index - offset + 1 {
                let count = self.get_count(i);
                if count != 0.0 {
                    let index: i64 = offset as i64 + i as i64;
                    serde::encode_signed_var_long(output, index - previous_index)?;
                    serde::encode_var_double(output, count)?;
                    previous_index = index;
                }
            }
        }
        Ok(())
    }
    fn decode_and_merge_with(
        &mut self,
        input: &mut Input,
        mode: BinEncodingMode,
    ) -> Result<(), Error> {
        match mode {
            BinEncodingMode::IndexDeltasAndCounts => {
                let num_bins = serde::decode_unsigned_var_long(input)?;
                let mut index: i64 = 0;
                let mut i = 0;
                while i < num_bins {
                    let index_delta = serde::decode_signed_var_long(input)?;
                    let count = serde::decode_var_double(input)?;
                    index += index_delta;
                    self.add(serde::i64_to_i32_exact(index)?, count);
                    i += 1;
                }

                Ok(())
            }

            BinEncodingMode::IndexDeltas => {
                let num_bins = serde::decode_unsigned_var_long(input)?;
                let mut index: i64 = 0;
                let mut i = 0;
                while i < num_bins {
                    let index_delta = serde::decode_signed_var_long(input)?;
                    index += index_delta;
                    self.add(serde::i64_to_i32_exact(index)?, 1.0);
                    i += 1;
                }
                Ok(())
            }

            BinEncodingMode::ContiguousCounts => {
                let num_bins = serde::decode_unsigned_var_long(input)?;
                let mut index: i64 = serde::decode_signed_var_long(input)?;
                let index_delta = serde::decode_signed_var_long(input)?;

                let mut i = 0;
                while i < num_bins {
                    let count = serde::decode_var_double(input)?;
                    self.add(serde::i64_to_i32_exact(index)?, count);
                    index += index_delta;
                    i += 1;
                }
                Ok(())
            }
        }
    }
    fn get_descending_stream(&self) -> Vec<(i32, f64)>;
    fn get_descending_iter(&self) -> StoreIter;
    fn get_ascending_iter(&self) -> StoreIter;
    fn get_sum(&self, index_mapping: &IndexMapping) -> f64 {
        let mut sum = 0.0;
        if self.is_empty() {
            return sum;
        }

        for i in self.get_min_index()..self.get_max_index() {
            let value = self.get_count(i - self.get_offset());
            if value != 0.0 {
                sum += index_mapping.value(i) * value;
            }
        }

        let last_count = self.get_count(self.get_max_index() - self.get_offset());
        if last_count != 0.0 {
            sum += index_mapping.value(self.get_max_index()) * last_count;
        }

        sum
    }
}

pub struct StoreIter<'a> {
    min_index: i32,
    max_index: i32,
    offset: i32,
    desc: bool,
    counts: &'a [f64],
}

impl<'a> StoreIter<'a> {
    pub fn new(
        min_index: i32,
        max_index: i32,
        offset: i32,
        desc: bool,
        counts: &'a [f64],
    ) -> StoreIter {
        StoreIter {
            desc,
            min_index,
            max_index,
            offset,
            counts,
        }
    }
}

impl<'a> Iterator for StoreIter<'a> {
    type Item = (i32, f64);
    fn next(&mut self) -> Option<Self::Item> {
        if self.desc {
            if self.max_index < self.min_index {
                return None;
            }

            let index = self.max_index;
            self.max_index -= 1;

            while self.max_index >= self.min_index {
                let count = self.counts[(self.max_index - self.offset) as usize];
                if count != 0.0 {
                    break;
                }
                self.max_index -= 1;
            }

            let count = self.counts[(index - self.offset) as usize];
            Some((index, count))
        } else {
            if self.min_index > self.max_index {
                return None;
            }

            let index = self.min_index;
            self.min_index += 1;

            while self.min_index <= self.max_index {
                let count = self.counts[(self.min_index - self.offset) as usize];
                if count != 0.0 {
                    break;
                }
                self.min_index += 1;
            }

            let count = self.counts[(index - self.offset) as usize];
            Some((index, count))
        }
    }
}

pub enum BinEncodingMode {
    IndexDeltasAndCounts = 1,
    IndexDeltas = 2,
    ContiguousCounts = 3,
}

impl BinEncodingMode {
    pub fn of_flag(marker: u8) -> Result<BinEncodingMode, Error> {
        let index = (marker >> 2) - 1;
        match index {
            0 => Ok(BinEncodingMode::IndexDeltasAndCounts),
            1 => Ok(BinEncodingMode::IndexDeltas),
            2 => Ok(BinEncodingMode::ContiguousCounts),
            _ => Err(Error::InvalidArgument("Unknown BinEncodingMode.")),
        }
    }

    pub fn to_flag(self, store_flag_type: FlagType) -> Flag {
        let sub_flag = self as u8;
        Flag::with_type(store_flag_type, sub_flag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collapsing_lowest_dense_store_add() {
        let mut store = CollapsingLowestDenseStore::with_capacity(10).unwrap();
        let indexes = vec![
            66, 14, 95, 71, 63, 28, 80, 54, 67, 41, 4, 24, 93, 73, 37, 37, 51, 49, 22, 90,
        ];
        for i in indexes {
            store.add(i, 1.0);
        }
        assert_eq!(95, store.get_max_index());
        assert_eq!(86, store.get_min_index());
        assert_eq!(20.0, store.get_total_count());
    }
}
