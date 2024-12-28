use nalgebra as na;

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

pub fn compute_blending_matrix(n: usize, cumulative: bool) -> na::DMatrix<f64> {
    let mut m = na::DMatrix::zeros(n, n);
    let factorial = (1..n).reduce(|acc, e| acc * e).unwrap() as f64;
    // let factorial = 1.0;
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

    m
}

#[cfg(test)]
mod tests {
    use crate::common::{binomial_coefficient, compute_blending_matrix, compute_blending_matrix_c};

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
        let m0 = compute_blending_matrix(N, false);
        let m1 = compute_blending_matrix_c::<N>(false);
        for i in 0..N {
            for j in 0..N {
                assert!((m0[(i, j)] - m1[i][j]).abs() < 1e-10);
            }
        }
    }
}