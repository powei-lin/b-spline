use tiny_solver::na;

pub const fn binomial_coefficient(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else {
        let mut r = 1;
        let mut d = 1;
        let mut n = n;
        while d < k + 1 {
            r *= n;
            n -= 1;
            r /= d;
            d += 1;
        }
        r
    }
}

const fn power(num: f64, i: usize) -> f64 {
    if i == 0 {
        return 1.0;
    }
    let mut n = num;
    let mut c = 1;
    while c < i {
        n *= num;
        c += 1;
    }
    n
}

pub const fn compute_blending_matrix_c<const N: usize>(cumulative: bool) -> [[f64; N]; N] {
    let mut m = [[0.0; N]; N];
    let mut n = 2;
    let mut factorial = 1.0;
    while n < N {
        factorial *= n as f64;
        n += 1;
    }
    let mut i = 0;
    while i < N {
        let mut j = 0;
        while j < N {
            let mut s = j;
            let mut sum = 0.0;
            while s < N {
                let sign = power(-1.0, s - j);
                sum += sign
                    * binomial_coefficient(N, s - j) as f64
                    * power((N - s) as f64 - 1.0, N - 1 - i);
                s += 1;
            }
            m[j][i] = binomial_coefficient(N - 1, N - 1 - i) as f64 * sum / factorial;
            j += 1;
        }
        i += 1;
    }
    let mut i = 0;
    if cumulative {
        while i < N {
            let mut j = i + 1;
            while j < N {
                let mut col = 0;
                while col < N {
                    m[i][col] += m[j][col];
                    col += 1;
                }
                j += 1;
            }
            i += 1;
        }
    }
    m
}

pub fn compute_blending_matrix<T: na::RealField>(n: usize, cumulative: bool) -> na::DMatrix<T> {
    let mut m = na::DMatrix::zeros(n, n);
    let factorial = (1..n).reduce(|acc, e| acc * e).unwrap() as f64;
    for i in 0..n {
        for j in 0..n {
            let sum: f64 = (j..n)
                .map(|s| {
                    let sign = (-1.0_f64).powi((s - j) as i32);
                    sign * binomial_coefficient(n, s - j) as f64
                        * ((n - s) as f64 - 1.0).powi((n - 1 - i) as i32)
                })
                .sum();
            m[(j, i)] = binomial_coefficient(n - 1, n - 1 - i) as f64 * sum / factorial;
        }
    }
    if cumulative {
        for i in 0..n {
            for j in i + 1..n {
                m.set_row(i, &(m.row(i) + m.row(j)));
            }
        }
    }
    m.cast()
}

pub fn compute_base_coefficients(n: usize) -> na::DMatrix<f64> {
    let mut base_coefficients = na::DMatrix::zeros(n, n);
    base_coefficients.row_mut(0).add_scalar_mut(1.0);
    let deg = n - 1;
    let mut order = deg;
    for row in 1..n {
        for i in (deg - order)..n {
            base_coefficients[(row, i)] =
                (order + i - deg) as f64 * base_coefficients[(row - 1, i)];
        }
        order -= 1;
    }
    base_coefficients
}

pub const fn compute_base_coefficients_c<const N: usize>() -> [[f64; N]; N] {
    let mut base_coefficients = [[0.0; N]; N];
    let mut col = 0;
    while col < N {
        base_coefficients[0][col] = 1.0;
        col += 1;
    }
    let deg = N - 1;
    let mut order = deg;

    let mut row = 1;
    while row < N {
        let mut i = deg - order;
        while i < N {
            base_coefficients[row][i] = (order + i - deg) as f64 * base_coefficients[row - 1][i];
            i += 1;
        }
        order -= 1;
        row += 1;
    }
    base_coefficients
}

#[cfg(test)]
mod tests {
    use crate::common::{
        binomial_coefficient, compute_base_coefficients, compute_base_coefficients_c,
        compute_blending_matrix, compute_blending_matrix_c,
    };

    #[test]
    fn test_binomial_coefficient() {
        assert!(binomial_coefficient(1, 1) == 1);
        assert!(binomial_coefficient(3, 2) == 3);
        assert!(binomial_coefficient(4, 2) == 6);
        assert!(binomial_coefficient(5, 2) == 10);
    }
    #[test]
    fn test_compute_blending_matrix() {
        const N: usize = 5;
        let m0 = compute_blending_matrix::<f64>(N, false);
        let m1 = compute_blending_matrix_c::<N>(false);
        for i in 0..N {
            for j in 0..N {
                assert!((m0[(i, j)] - m1[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_compute_base_coefficients() {
        const N: usize = 5;
        let m0 = compute_base_coefficients(N);
        let m1 = compute_base_coefficients_c::<N>();
        for i in 0..N {
            for j in 0..N {
                assert!((m0[(i, j)] - m1[i][j]).abs() < 1e-10);
            }
        }
    }
}
