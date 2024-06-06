use crate::sketch::{Flag, FlagType};
use crate::Error;

use crate::output::Output;

#[derive(PartialEq, Debug, Clone)]
pub enum IndexMapping {
    LogarithmicMapping(f64, f64, f64, f64),
}

const LOGARITHMIC_MAPPING_CORRECTING_FACTOR: f64 = 1.0;
const LOGARITHMIC_MAPPING_BASE: f64 = std::f64::consts::E;

impl IndexMapping {
    pub fn layout(&self) -> IndexMappingLayout {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                _index_offset,
                _multiplier,
                _relative_accuracy,
            ) => IndexMappingLayout::LOG,
        }
    }

    pub fn gamma(&self) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                gamma,
                _index_offset,
                _multiplier,
                _relative_accuracy,
            ) => *gamma,
        }
    }

    pub fn index_offset(&self) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                index_offset,
                _multiplier,
                _relative_accuracy,
            ) => *index_offset,
        }
    }

    pub fn multiplier(&self) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                _index_offset,
                multiplier,
                _relative_accuracy,
            ) => *multiplier,
        }
    }

    pub fn relative_accuracy(&self) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                _index_offset,
                _multiplier,
                relative_accuracy,
            ) => *relative_accuracy,
        }
    }

    fn log(&self, value: f64) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                _index_offset,
                _multiplier,
                _relative_accuracy,
            ) => value.ln(),
        }
    }

    fn log_inverse(&self, index: f64) -> f64 {
        match self {
            IndexMapping::LogarithmicMapping(
                _gamma,
                _index_offset,
                _multiplier,
                _relative_accuracy,
            ) => index.exp(),
        }
    }

    pub fn index(&self, value: f64) -> i32 {
        let index: f64 = self.log(value) * self.multiplier() + self.index_offset();
        if index >= 0.0 {
            index as i32
        } else {
            (index - 1.0) as i32
        }
    }

    pub fn value(&self, index: i32) -> f64 {
        self.lower_bound(index) * (1.0 + self.relative_accuracy())
    }

    fn lower_bound(&self, index: i32) -> f64 {
        self.log_inverse((index as f64 - self.index_offset()) / self.multiplier())
    }

    #[allow(dead_code)]
    fn upper_bound(&self, index: i32) -> f64 {
        self.lower_bound(index + 1)
    }

    pub(crate) fn min_indexable_value(&self) -> f64 {
        f64::max(
            f64::powf(
                2.0,
                (i32::MIN as f64 - self.index_offset()) / self.multiplier() + 1.0,
            ),
            f64::MIN_POSITIVE * (1.0 + self.relative_accuracy()) / (1.0 - self.relative_accuracy()),
        )
    }

    pub(crate) fn max_indexable_value(&self) -> f64 {
        f64::max(
            f64::powf(
                2.0,
                (i32::MAX as f64 - self.index_offset()) / self.multiplier() - 1.0,
            ),
            f64::MAX / (1.0 + self.relative_accuracy()),
        )
    }

    pub fn encode(&self, output: &mut Output) -> Result<(), Error> {
        self.layout().to_flag().encode(output)?;
        output.write_double_le(self.gamma())?;
        output.write_double_le(self.index_offset())?;
        Ok(())
    }

    pub fn with_relative_accuracy(
        index_layout: IndexMappingLayout,
        relative_accuracy: f64,
    ) -> Result<IndexMapping, Error> {
        if relative_accuracy <= 0.0 || relative_accuracy >= 1.0 {
            return Err(Error::InvalidArgument(
                "The relative accuracy must be between 0 and 1.",
            ));
        }

        match index_layout {
            IndexMappingLayout::LOG => {
                if relative_accuracy <= 0.0 || relative_accuracy >= 1.0 {
                    return Err(Error::InvalidArgument(
                        "The relative accuracy must be between 0 and 1.",
                    ));
                }

                let gamma =
                    calculate_gamma(relative_accuracy, LOGARITHMIC_MAPPING_CORRECTING_FACTOR);
                let index_offset: f64 = 0.0;
                let multiplier = LOGARITHMIC_MAPPING_BASE.ln() / (gamma - 1.0).ln_1p();
                let relative_accuracy = calculate_relative_accuracy(gamma, 1.0);
                Ok(IndexMapping::LogarithmicMapping(
                    gamma,
                    index_offset,
                    multiplier,
                    relative_accuracy,
                ))
            }
        }
    }

    pub fn with_gamma_offset(
        index_layout: IndexMappingLayout,
        gamma: f64,
        index_offset: f64,
    ) -> Result<IndexMapping, Error> {
        match index_layout {
            IndexMappingLayout::LOG => {
                let multiplier = LOGARITHMIC_MAPPING_BASE.ln() / gamma.ln();
                let relative_accuracy =
                    calculate_relative_accuracy(gamma, LOGARITHMIC_MAPPING_CORRECTING_FACTOR);
                Ok(IndexMapping::LogarithmicMapping(
                    gamma,
                    index_offset,
                    multiplier,
                    relative_accuracy,
                ))
            }
        }
    }
}

pub enum IndexMappingLayout {
    LOG = 0,
}

impl IndexMappingLayout {
    pub fn of_flag(flag: &Flag) -> Result<IndexMappingLayout, Error> {
        let index = flag.get_marker() >> 2;
        match index {
            0 => Ok(IndexMappingLayout::LOG),
            _ => Err(Error::InvalidArgument("Unknown Index Flag.")),
        }
    }

    pub fn to_flag(self) -> Flag {
        let sub_flag = self as u8;
        Flag::with_type(FlagType::IndexMapping, sub_flag)
    }
}

fn calculate_relative_accuracy(gamma: f64, correcting_factor: f64) -> f64 {
    let exact_log_gamma = gamma.powf(correcting_factor);
    (exact_log_gamma - 1.0) / (exact_log_gamma + 1.0)
}

fn calculate_gamma(relative_accuracy: f64, correcting_factor: f64) -> f64 {
    let exact_log_gamma = (1.0 + relative_accuracy) / (1.0 - relative_accuracy);
    exact_log_gamma.powf(1.0 / correcting_factor)
}

#[cfg(test)]
mod tests {
    use crate::index_mapping::IndexMapping;
    use crate::index_mapping::IndexMappingLayout::LOG;

    const TEST_GAMMAS: [f64; 3] = [1.0 + 1e-6, 1.02, 1.5];
    const TEST_INDEX_OFFSETS: [f64; 4] = [0.0, 1.0, -12.23, 7768.3];
    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_logarithmic_mapping_offset() {
        for gamma in TEST_GAMMAS {
            for index_offset in TEST_INDEX_OFFSETS {
                let index_mapping =
                    IndexMapping::with_gamma_offset(LOG, gamma, index_offset).unwrap();
                let index_of1 = index_mapping.index(1.0) as f64;
                // If 1 is on a bucket boundary, its associated index can be either of the ones of the previous
                // and the next buckets.
                assert!(index_offset.ceil() - 1.0 <= index_of1);
                assert!(index_of1 <= index_offset.floor());
            }
        }
    }

    #[test]
    fn test_cubically_interpolated_mapping_validity_manual_check() {
        let d0: f64 = -0.37469387755102035;
        let d1: f64 = 0.8904489795918369;
        let s1 = d1 * d1 - 4.0_f64 * d0 * d0 * d0;
        let s2 = s1.sqrt();
        let s3 = d1 - s2;
        let s4 = s3 / 2.0;
        let s5 = s4.cbrt();
        eprintln!(
            "test_cubically_interpolated_mapping_validity_manual_check: {} {} {} {} {}",
            s1, s2, s3, s4, s5
        );
    }

    #[test]
    fn test_logarithmic_mapping_validity() {
        let mapping = IndexMapping::with_relative_accuracy(LOG, 1e-2).unwrap();

        println!("LogarithmicMapping: {:?}", mapping);

        let min_index = -50;
        let max_index = 50;

        let mut index = min_index;
        let mut bound = mapping.upper_bound(index - 1);

        while index <= max_index {
            println!(
                "test_logarithmic_mapping_validity {} {} {} {}",
                index,
                mapping.value(index),
                mapping.lower_bound(index),
                mapping.upper_bound(index)
            );

            assert!(f64::abs(mapping.lower_bound(index) - bound) <= 1e10);
            assert!(mapping.value(index) >= mapping.lower_bound(index));

            assert!(mapping.upper_bound(index) >= mapping.value(index));

            assert!(mapping.index(mapping.lower_bound(index) - EPSILON) < index);
            assert!(mapping.index(mapping.lower_bound(index) + EPSILON) >= index);

            assert!(mapping.index(mapping.upper_bound(index) - EPSILON) <= index);
            assert!(mapping.index(mapping.upper_bound(index) + EPSILON) > index);

            bound = mapping.upper_bound(index);
            index += 1;
        }
    }

    #[test]
    fn test_logarithmic_mapping_index() {
        let mapping = IndexMapping::with_relative_accuracy(LOG, 2e-2).unwrap();
        let values: Vec<f64> = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
            17.0, 18.0, 19.0, 228.0, 484.0, 499.0, 559.0, 584.0, 629.0, 722.0, 730.0, 777.0, 805.0,
            846.0, 896.0, 997.0, 999.0, 1065.0, 1178.0, 1189.0, 1218.0, 1255.0, 1308.0, 1343.0,
            1438.0, 1819.0, 2185.0, 2224.0, 2478.0, 2574.0, 2601.0, 2745.0, 2950.0, 3013.0, 3043.0,
            3064.0, 3116.0, 3188.0, 3224.0, 3254.0, 3390.0, 3476.0, 3543.0, 3836.0, 3921.0, 4014.0,
            4074.0, 4332.0, 4344.0, 4456.0, 4736.0, 4984.0, 5219.0, 5244.0, 5259.0, 5341.0, 5467.0,
            5536.0, 5600.0, 6054.0, 6061.0, 6118.0, 6137.0, 6222.0, 6263.0, 6320.0, 6454.0, 6499.0,
            6732.0, 6922.0, 6988.0, 7047.0, 7057.0, 7202.0, 7205.0, 7330.0, 7507.0, 7616.0, 7971.0,
            8056.0, 8381.0, 8416.0, 8684.0, 8784.0, 8790.0, 8823.0, 8841.0, 8945.0, 8967.0, 8982.0,
            9142.0, 9181.0, 9284.0, 9320.0, 9331.0, 9596.0, 9699.0, 9850.0, 9884.0, 9947.0,
        ];
        let indexes = vec![
            0, 17, 27, 34, 40, 44, 48, 51, 54, 57, 59, 62, 64, 65, 67, 69, 70, 72, 73, 135, 154,
            155, 158, 159, 161, 164, 164, 166, 167, 168, 169, 172, 172, 174, 176, 176, 177, 178,
            179, 180, 181, 187, 192, 192, 195, 196, 196, 197, 199, 200, 200, 200, 201, 201, 201,
            202, 203, 203, 204, 206, 206, 207, 207, 209, 209, 210, 211, 212, 213, 214, 214, 214,
            215, 215, 215, 217, 217, 217, 218, 218, 218, 218, 219, 219, 220, 221, 221, 221, 221,
            222, 222, 222, 223, 223, 224, 224, 225, 225, 226, 226, 227, 227, 227, 227, 227, 227,
            227, 228, 228, 228, 228, 229, 229, 229, 229, 230,
        ];
        for i in 0..values.len() {
            assert_eq!(indexes[i], mapping.index(values[i]));
        }
    }
}
