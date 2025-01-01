/// 判断两个浮点数是否相等
pub(crate) fn approx_eq(a: f64, b: f64) -> bool {
    if a.is_infinite() || b.is_infinite() {
        false
    } else {
        (a - b).abs() < f64::EPSILON
    }
}

/// 秦九韶算法求多项式的值
///
/// # Arguments
///
/// - `x`:  自变量的值
/// - `arr`:  多项式的系数数组, 按次数降序
/// 
/// # Example 
/// 
/// ```
/// use special_functions::utils::eval_poly;
/// eval_poly(0.5, &[16., 0., 20., 0., 5., 0.]) // 6th first-kind Chebyshev polynomial
/// ```
pub fn eval_poly(x: f64, arr: &[f64]) -> f64 {
    arr.iter().fold(0.0, |acc, &a| acc * x + a)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_approx_eq() {
        assert_ne!(0.1 + 0.2, 0.3);
        assert!(approx_eq(0.1 + 0.2, 0.3));
    }
    #[test]
    fn test_eval_poly() {
        let arr = [
            0.3198453915289723,
            0.9076227501539942,
            0.40138509410337553,
            0.9088787482769067,
            0.7563007138750291,
        ];
        let x = 0.35625260496659283;
        let result = eval_poly(x, &arr);
        assert!(approx_eq(result, 1.1772226211231838));
        assert!(approx_eq(
            eval_poly(2.7172900350129723, &[4., 2., 9., 8.]),
            127.47717934998103
        ));
    }
}
