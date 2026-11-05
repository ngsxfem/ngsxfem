#pragma once

/*
  Plain C++ interface to the quadrature algorithms of the bundled Algoim library
  (external_dependencies/algoim, R. I. Saye, https://algoim.github.io/).

  Algoim computes high-order quadrature rules for domains in the unit
  hypercube [0,1]^dim that are implicitly defined by a multivariate polynomial
  in Bernstein form (R. I. Saye, J. Comput. Phys. 448, 110720, 2022).

  Algoim headers must only be included in one translation unit
  (cutint/algoim_wrapper.cpp), hence this interface does not expose any Algoim
  (or NGSolve) types. The NGSolve-facing part is in cutint/algoimcutrule.hpp.
*/

#include <vector>

namespace xintegration
{
  namespace algoim_wrapper
  {
    /// classification of the cell w.r.t. the zero level of phi
    enum CUT_STATUS { UNCUT_NEGATIVE = -1, CUT = 0, UNCUT_POSITIVE = 1 };

    /// classification of the cell w.r.t. a part defined by several level sets:
    /// empty, a proper part of the cell (cut), or the whole cell
    enum MULTI_STATUS { EMPTY = -2, PARTIAL = 0, FULL = 2 };

    /// Interpolation points (in [0,1]^dim) for a tensor-product polynomial of
    /// degree 'degree' in each coordinate direction, (degree+1)^dim points.
    /// On exit, points[i*dim+d] is coordinate d of point i. The ordering is
    /// the one expected for 'nodal_values' in CutQuadrature.
    void InterpolationPoints (int dim, int degree, std::vector<double> & points);

    /// Quadrature rule on the part of the cell C given by
    ///   domain == 0 : { phi < 0 },
    ///   domain == 1 : { phi > 0 },
    ///   domain == 2 : { phi = 0 } (surface quadrature; for dim == 1 the
    ///                 roots of phi with weight 1),
    /// where phi is the tensor-product polynomial (of degree 'degree' in each
    /// direction) that interpolates 'nodal_values' at the points of
    /// InterpolationPoints. (The domain values coincide with ngsxfem's NEG,
    /// POS, IF.) The cell C is the unit hypercube [0,1]^dim or, if simplex is
    /// true, the unit simplex {x in [0,1]^dim : x_0 + ... + x_{dim-1} < 1}
    /// (treated as a second polynomial constraint on [0,1]^dim; the level set
    /// values are needed in the whole hypercube in this case).
    ///
    /// q is the number of Gauss-Legendre points on every one-dimensional
    /// sub-interval of Algoim's dimension-reduction approach.
    ///
    /// Algoim is applied on sub-boxes of [0,1]^dim: [0,1]^dim is uniformly
    /// subdivided subdivlvl times (once more for simplices) and every sub-box,
    /// in which Algoim detects vertical tangents of the interface (w.r.t. any
    /// coordinate direction), is subdivided once more. (The subdivisions are
    /// affine.) On meshes that are coarse compared to the curvature of the
    /// interface, (nearly) vertical tangents can still limit the convergence
    /// w.r.t. q in single cells; subdivlvl > 0 helps in this case.
    ///
    /// If phi has (provably) no zero in C, nothing is computed and
    /// UNCUT_NEGATIVE / UNCUT_POSITIVE is returned. Otherwise CUT is returned
    /// and points (n*dim), weights (n) and - for domain == 2 - the gradient of
    /// phi at the points (n*dim, not normalized) are filled. Surface weights
    /// are with respect to the surface measure in [0,1]^dim.
    /// Quadrature rule for the part of the cell C where n level set functions
    /// phi_l (interpolating nodal_values[l] with degree degrees[l], see
    /// CutQuadrature) satisfy the sign conditions domains[l] (0: phi_l < 0,
    /// 1: phi_l > 0, 2: phi_l = 0; at most one interface condition). Returns
    /// EMPTY, FULL (nothing is computed) or PARTIAL (with the rule as in
    /// CutQuadrature; for an interface, the gradients are those of the level set
    /// with the interface condition).
    ///
    /// The cell C is the unit hypercube (simplex_dim == 0) or the product of
    /// the unit simplex in the first simplex_dim coordinates with the unit
    /// hypercube in the remaining ones (e.g. a triangle or tetrahedron times a
    /// time interval); dim <= 4. The subdivisions are as described for
    /// CutQuadrature, the additional subdivision of simplices only applies to
    /// the simplex directions.
    ///
    /// last_axis_time: the last coordinate is the time of a space-time cell.
    /// Then, Algoim's choice of the height directions (in all levels of the
    /// dimension reduction) prefers the spatial directions: the time is used
    /// as height direction only if it is clearly the steepest direction (by a
    /// factor LastAxisScoreFactor^{-1}) or if the spatial directions have
    /// vertical tangents. As for the standard space-time rules, time is then
    /// mostly the variable of the outermost integral. The height functions in
    /// time direction are often close to singular, e.g. if an interface
    /// changes its direction of motion in the cell, whereas the spatial height
    /// functions are smooth (and exact for level sets that are linear in
    /// space).
    MULTI_STATUS MultiCutQuadrature (int dim, int n, const int * degrees, const double * const * nodal_values,
                                     const int * domains, int simplex_dim, int q, int subdivlvl,
                                     std::vector<double> & points,
                                     std::vector<double> & weights,
                                     std::vector<double> & gradients,
                                     bool last_axis_time = false);

    /// factor for the score of the time axis, see MultiCutQuadrature
    constexpr double LastAxisScoreFactor = 0.25;

    CUT_STATUS CutQuadrature (int dim, int degree, const double * nodal_values, bool simplex,
                              int domain, int q, int subdivlvl,
                              std::vector<double> & points,
                              std::vector<double> & weights,
                              std::vector<double> & gradients);
  }
}
