/// 秦九韶算法求多项式的值
///
/// # Arguments
///
/// - `x`:  自变量的值
/// - `arr`:  多项式的系数数组, 按次数降序
pub fn eval_poly(x: f64, arr: &[f64]) -> f64 {
    arr.iter().rev().fold(1.0, |acc, &a| acc * x + a)
}