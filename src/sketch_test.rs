use crate::output::Output;
use crate::serde;
use crate::{Bucket, DDSketch};

const RELATIVE_ACCURACY: f64 = 0.02;
const GAMMA: f64 = (1.0 + RELATIVE_ACCURACY) / (1.0 - RELATIVE_ACCURACY);

fn write_mapping(output: &mut Output) {
    output.write_byte(2).unwrap();
    output.write_double_le(GAMMA).unwrap();
    output.write_double_le(0.0).unwrap();
}

fn write_empty_clickhouse_store(output: &mut Output, store_flag: u8) {
    output.write_byte(store_flag).unwrap();
    output.write_byte(12).unwrap();
    serde::encode_unsigned_var_long(output, 0).unwrap();
    serde::encode_signed_var_long(output, 0).unwrap();
    serde::encode_signed_var_long(output, 1).unwrap();
}

fn clickhouse_payload_with_bin_count(bin_count: f64) -> Vec<u8> {
    let mut output = Output::with_capacity(64);
    write_mapping(&mut output);

    output.write_byte(1).unwrap();
    output.write_byte(4).unwrap();
    serde::encode_unsigned_var_long(&mut output, 1).unwrap();
    serde::encode_signed_var_long(&mut output, 0).unwrap();
    output.write_double_le(bin_count).unwrap();

    write_empty_clickhouse_store(&mut output, 3);

    output.write_byte(4).unwrap();
    output.write_double_le(0.0).unwrap();
    output.trim()
}

fn clickhouse_payload_with_zero_count(zero_count: f64) -> Vec<u8> {
    let mut output = Output::with_capacity(64);
    write_mapping(&mut output);
    write_empty_clickhouse_store(&mut output, 1);
    write_empty_clickhouse_store(&mut output, 3);
    output.write_byte(4).unwrap();
    output.write_double_le(zero_count).unwrap();
    output.trim()
}

fn normal_payload_with_bin_count(bin_count: f64) -> Vec<u8> {
    let mut output = Output::with_capacity(64);
    write_mapping(&mut output);
    output.write_byte(5).unwrap();
    serde::encode_unsigned_var_long(&mut output, 1).unwrap();
    serde::encode_signed_var_long(&mut output, 0).unwrap();
    serde::encode_var_double(&mut output, bin_count).unwrap();
    output.trim()
}

fn normal_payload_with_zero_count(zero_count: f64) -> Vec<u8> {
    let mut output = Output::with_capacity(64);
    write_mapping(&mut output);
    output.write_byte(4).unwrap();
    serde::encode_var_double(&mut output, zero_count).unwrap();
    output.trim()
}

#[test]
fn test_sketch_quantile_3() {
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.accept(1.0);
    sketch.accept(2.0);
    sketch.accept(3.0);
    sketch.accept(4.0);
    sketch.accept(5.0);

    assert!((f64::abs(sketch.get_value_at_quantile(0.0).unwrap() - 1.0) / 1.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(0.5).unwrap() - 3.0) / 3.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(1.0).unwrap() - 5.0) / 5.0) < 0.021);

    let encoded = sketch.encode().unwrap();
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.decode_and_merge_with(&encoded).unwrap();
    assert!((f64::abs(sketch.get_value_at_quantile(0.0).unwrap() - 1.0) / 1.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(0.5).unwrap() - 3.0) / 3.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(1.0).unwrap() - 5.0) / 5.0) < 0.021);

    let encoded = sketch.encode_clickhouse().unwrap();
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.decode_clickhouse_and_merge_with(&encoded).unwrap();
    assert!((f64::abs(sketch.get_value_at_quantile(0.0).unwrap() - 1.0) / 1.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(0.5).unwrap() - 3.0) / 3.0) < 0.021);
    assert!((f64::abs(sketch.get_value_at_quantile(1.0).unwrap() - 5.0) / 5.0) < 0.021);
}

#[test]
fn test_sketch_decode_3() {
    let input = vec![
        2, 42, 120, 57, 5, 47, 167, 240, 63, 0, 0, 0, 0, 0, 0, 0, 0, 13, 50, 130, 1, 2, 136, 32, 0,
        3, 0, 0, 0, 3, 0, 2, 0, 0, 3, 3, 2, 2, 3, 3, 2, 0, 0, 0, 0, 2, 0, 2, 2, 2, 4, 4, 132, 64,
        0, 4, 2, 0, 2, 2, 3, 132, 64, 4, 132, 64, 4, 2, 2, 0, 6, 4, 6, 132, 64, 2, 6,
    ];
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(2e-2, 50).unwrap();
    sketch.decode_and_merge_with(&input).unwrap();
    assert_eq!(sketch.get_count(), 100.0);

    let encoded = sketch.encode().unwrap();
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.decode_and_merge_with(&encoded).unwrap();
    assert_eq!(sketch.get_count(), 100.0);

    let encoded = sketch.encode_clickhouse().unwrap();
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.decode_clickhouse_and_merge_with(&encoded).unwrap();
    assert_eq!(sketch.get_count(), 100.0);
}

#[test]
fn test_sketch_buckets() {
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    let values = [-10.0, -1.0, 0.0, 1.0, 10.0];
    for value in values {
        sketch.accept(value);
    }

    let buckets = sketch.get_buckets();
    assert_eq!(buckets.len(), 5);

    let total_count: f64 = buckets.iter().map(|bucket| bucket.count).sum();
    assert_eq!(total_count, 5.0);

    let mut last_begin = f64::NEG_INFINITY;
    for bucket in &buckets {
        assert!(bucket.range_begin <= bucket.range_end);
        assert!(bucket.range_begin >= last_begin);
        last_begin = bucket.range_begin;
    }

    for value in values {
        let bucket = find_bucket(&buckets, value).unwrap();
        assert_eq!(bucket.count, 1.0);
    }
}

#[test]
fn test_accept_rejects_non_finite_values() {
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.accept(f64::NAN);
    sketch.accept(f64::INFINITY);
    sketch.accept(f64::NEG_INFINITY);
    assert!(sketch.is_empty());
}

#[test]
fn test_accept_with_count_uses_count_and_rejects_invalid_counts() {
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.accept_with_count(1.0, 3.0);
    sketch.accept_with_count(2.0, f64::NAN);
    sketch.accept_with_count(3.0, f64::INFINITY);
    sketch.accept_with_count(4.0, -1.0);
    sketch.accept_with_count(5.0, 0.0);
    assert_eq!(sketch.get_count(), 3.0);

    let encoded = sketch.encode().unwrap();
    let mut decoded = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    decoded.decode_and_merge_with(&encoded).unwrap();
    assert_eq!(decoded.get_count(), 3.0);
}

#[test]
fn test_encode_rejects_invalid_total_count() {
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(0.02, 100).unwrap();
    sketch.accept_with_count(1.0, f64::MAX);
    sketch.accept_with_count(2.0, f64::MAX);
    assert!(sketch.encode().is_err());
    assert!(sketch.encode_clickhouse().is_err());
}

#[test]
fn test_clickhouse_decode_rejects_invalid_mapping() {
    let mut encoded = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100)
        .unwrap()
        .encode_clickhouse()
        .unwrap();
    encoded[1..9].copy_from_slice(&f64::NAN.to_le_bytes());

    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
    assert!(sketch.decode_clickhouse_and_merge_with(&encoded).is_err());

    let mut encoded = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100)
        .unwrap()
        .encode_clickhouse()
        .unwrap();
    encoded[9..17].copy_from_slice(&f64::NAN.to_le_bytes());

    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
    assert!(sketch.decode_clickhouse_and_merge_with(&encoded).is_err());
}

#[test]
fn test_clickhouse_decode_rejects_invalid_counts() {
    for count in [f64::NAN, f64::INFINITY, -1.0] {
        let mut sketch =
            DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
        assert!(sketch
            .decode_clickhouse_and_merge_with(&clickhouse_payload_with_bin_count(count))
            .is_err());

        let mut sketch =
            DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
        assert!(sketch
            .decode_clickhouse_and_merge_with(&clickhouse_payload_with_zero_count(count))
            .is_err());
    }
}

#[test]
fn test_clickhouse_decode_rejects_too_many_bins() {
    let mut output = Output::with_capacity(64);
    write_mapping(&mut output);
    output.write_byte(1).unwrap();
    output.write_byte(12).unwrap();
    serde::encode_unsigned_var_long(&mut output, 65_537).unwrap();
    serde::encode_signed_var_long(&mut output, 0).unwrap();
    serde::encode_signed_var_long(&mut output, 1).unwrap();

    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
    assert!(sketch
        .decode_clickhouse_and_merge_with(&output.trim())
        .is_err());
}

#[test]
fn test_clickhouse_encode_empty_stores_round_trips() {
    let encoded = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100)
        .unwrap()
        .encode_clickhouse()
        .unwrap();
    let mut sketch = DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
    sketch.decode_clickhouse_and_merge_with(&encoded).unwrap();
    assert_eq!(sketch.get_count(), 0.0);
}

#[test]
fn test_normal_decode_rejects_invalid_counts() {
    for count in [f64::NAN, f64::INFINITY, -1.0] {
        let mut sketch =
            DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
        assert!(sketch
            .decode_and_merge_with(&normal_payload_with_bin_count(count))
            .is_err());

        let mut sketch =
            DDSketch::logarithmic_collapsing_lowest_dense(RELATIVE_ACCURACY, 100).unwrap();
        assert!(sketch
            .decode_and_merge_with(&normal_payload_with_zero_count(count))
            .is_err());
    }
}

fn find_bucket(buckets: &[Bucket], value: f64) -> Option<&Bucket> {
    buckets
        .iter()
        .find(|bucket| bucket.range_begin <= value && value <= bucket.range_end)
}
