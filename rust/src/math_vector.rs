///
/// Math functions for N-dimensional vectors
///
/// All functions work with arbitrary-dimension vectors represented as slices.
/// Functions returning arrays use DIMS_MAX for stack allocation, but only
/// operate on the first `dims` elements.

/// Maximum supported dimensions. Arrays are stack-allocated to this size.
/// Configurable via DIMS_MAX environment variable at build time (default: 32).
pub const DIMS_MAX: usize = {
    const fn parse_usize(s: &str) -> usize {
        let mut result = 0usize;
        let bytes = s.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            result = result * 10 + (bytes[i] - b'0') as usize;
            i += 1;
        }
        result
    }
    parse_usize(env!("DIMS_MAX"))
};

const _: () = assert!(DIMS_MAX > 0, "DIMS_MAX must be greater than 0");

const EPS: f64 = 1e-8;

pub fn sq(d: f64) -> f64 { d * d }

pub fn is_finite_vn(v0: &[f64]) -> bool {
    v0.iter().all(|f| f.is_finite())
}

pub fn zero_vn(v0: &mut [f64]) {
    v0.fill(0.0);
}

pub fn negated_vn(v0: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = -v0[j];
    }
    out
}

pub fn copy_vnvn(v0: &mut [f64], v1: &[f64]) {
    debug_assert_eq!(v0.len(), v1.len());
    v0.copy_from_slice(v1);
}

pub fn dot_vnvn(v0: &[f64], v1: &[f64]) -> f64 {
    debug_assert_eq!(v0.len(), v1.len());
    v0.iter().zip(v1).map(|(a, b)| a * b).sum()
}

pub fn add_vnvn(v0: &[f64], v1: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] + v1[j];
    }
    out
}

pub fn sub_vnvn(v0: &[f64], v1: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] - v1[j];
    }
    out
}

#[allow(dead_code)]
pub fn mid_vnvn(v0: &[f64], v1: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = (v0[j] + v1[j]) * 0.5;
    }
    out
}

#[allow(dead_code)]
pub fn interp_vnvn(v0: &[f64], v1: &[f64], t: f64, dims: usize) -> [f64; DIMS_MAX] {
    let s = 1.0 - t;
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] * s + v1[j] * t;
    }
    out
}

pub fn madd_vnvn_fl(v0: &[f64], v1: &[f64], f: f64, dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] + v1[j] * f;
    }
    out
}

pub fn msub_vnvn_fl(v0: &[f64], v1: &[f64], f: f64, dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] - v1[j] * f;
    }
    out
}

pub fn mul_vn_fl(v0: &[f64], f: f64, dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    for j in 0..dims {
        out[j] = v0[j] * f;
    }
    out
}

fn imul_vn_fl(v0: &mut [f64], f: f64) {
    v0.iter_mut().for_each(|v| *v *= f);
}

pub fn len_squared_vnvn(v0: &[f64], v1: &[f64]) -> f64 {
    debug_assert_eq!(v0.len(), v1.len());
    v0.iter().zip(v1).map(|(a, b)| sq(a - b)).sum()
}

pub fn len_squared_vn(v0: &[f64]) -> f64 {
    v0.iter().map(|v| sq(*v)).sum()
}

pub fn len_vnvn(v0: &[f64], v1: &[f64]) -> f64 {
    len_squared_vnvn(v0, v1).sqrt()
}

pub fn len_squared_negated_vnvn(v0: &[f64], v1: &[f64]) -> f64 {
    debug_assert_eq!(v0.len(), v1.len());
    v0.iter().zip(v1).map(|(a, b)| sq(a + b)).sum()
}

// special case, save us negating a copy, then getting the length
pub fn len_negated_vnvn(v0: &[f64], v1: &[f64]) -> f64 {
    len_squared_negated_vnvn(v0, v1).sqrt()
}

pub fn normalize_vn(v0: &mut [f64]) -> f64 {
    let mut d = len_squared_vn(v0);
    if (d != 0.0) && ({d = d.sqrt(); d} != 0.0) {
        imul_vn_fl(v0, 1.0 / d);
    }
    d
}

pub fn normalized_vn(v0: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut out = [0.0; DIMS_MAX];
    out[..dims].copy_from_slice(&v0[..dims]);
    normalize_vn(&mut out[..dims]);
    out
}

// v_out = (v0 - v1).normalized()
pub fn normalized_vnvn(v0: &[f64], v1: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mut v = sub_vnvn(v0, v1, dims);
    normalize_vn(&mut v[..dims]);
    v
}

pub fn normalized_vnvn_with_len(v0: &[f64], v1: &[f64], dims: usize) -> ([f64; DIMS_MAX], f64) {
    let mut v = sub_vnvn(v0, v1, dims);
    let d = normalize_vn(&mut v[..dims]);
    (v, d)
}

pub fn is_almost_zero_ex(val: f64, eps: f64) -> bool {
    (-eps < val) && (val < eps)
}

pub fn is_almost_zero(val: f64) -> bool {
    is_almost_zero_ex(val, EPS)
}

pub fn project_vnvn_normalized(p: &[f64], v_proj: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let mul = dot_vnvn(&p[..dims], &v_proj[..dims]);
    mul_vn_fl(v_proj, mul, dims)
}

pub fn project_plane_vnvn_normalized(v: &[f64], v_plane: &[f64], dims: usize) -> [f64; DIMS_MAX] {
    let proj = project_vnvn_normalized(v, v_plane, dims);
    sub_vnvn(v, &proj[..dims], dims)
}
