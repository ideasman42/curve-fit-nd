//! Curve Fitting Module
//!
//! Provides algorithms for fitting cubic bezier curves to point data.
//! This includes single-segment fitting, polygon fitting, and iterative
//! refinement with corner detection.

/// Polygon to bezier curve conversion with parallel processing support.
mod curve_fit_from_polys;
/// Single-segment cubic bezier fitting with multiple solver strategies.
mod curve_fit_cubic;
/// Iterative curve refinement: knot removal, repositioning, and corner detection.
mod curve_fit_cubic_refit;

/// Re-export math_vector for external access to vector operations.
#[allow(unused_imports)]
pub use crate::math_vector;

/// Public API for polygon-to-curve fitting.
#[allow(unused_imports)]
pub use self::curve_fit_from_polys::{
    TraceMode,
    fit_poly_single,
    fit_poly_list,
};

/// Re-export refinement types for external use.
#[allow(unused_imports)]
pub use self::curve_fit_cubic_refit::{
    Knot,
    PointData,
    refine_remove,
    refine_refit,
    refine_corner,
};

