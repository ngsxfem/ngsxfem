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
    /// classification of the unit hypercube w.r.t. the zero level of phi
    enum CUT_STATUS { UNCUT_NEGATIVE = -1, CUT = 0, UNCUT_POSITIVE = 1 };

    /// Interpolation points (in [0,1]^dim) for a tensor-product polynomial of
    /// degree 'degree' in each coordinate direction, (degree+1)^dim points.
    /// On exit, points[i*dim+d] is coordinate d of point i. The ordering is
    /// the one expected for 'nodal_values' in CutQuadrature.
    void InterpolationPoints (int dim, int degree, std::vector<double> & points);

    /// Quadrature rule on the part of the cell C given by
    ///   domain == 0 : { phi < 0 },
    ///   domain == 1 : { phi > 0 },
    ///   domain == 2 : { phi = 0 } (surface quadrature, dim > 1),
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
    /// in which Algoim detects vertical tangents (w.r.t. the height directions
    /// it chooses), is subdivided once more. (The subdivisions are affine.)
    ///
    /// If phi has (provably) no zero in C, nothing is computed and
    /// UNCUT_NEGATIVE / UNCUT_POSITIVE is returned. Otherwise CUT is returned
    /// and points (n*dim), weights (n) and - for domain == 2 - the gradient of
    /// phi at the points (n*dim, not normalized) are filled. Surface weights
    /// are with respect to the surface measure in [0,1]^dim.
    CUT_STATUS CutQuadrature (int dim, int degree, const double * nodal_values, bool simplex,
                              int domain, int q, int subdivlvl,
                              std::vector<double> & points,
                              std::vector<double> & weights,
                              std::vector<double> & gradients);
  }
}
