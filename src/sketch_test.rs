use crate::DDSketch;

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
