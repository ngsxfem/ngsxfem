#pragma once

/*
  Minimal LAPACKE replacement for the bundled Algoim library
  (external_dependencies/algoim).

  Algoim includes <lapacke.h> and uses exactly two LAPACKE routines:
    * LAPACKE_dgesvd  (SVD of the Bernstein-Vandermonde matrix used for
                       polynomial interpolation),
    * LAPACKE_dggevx  (generalized eigenvalues, used as a robust fallback
                       for the real roots of univariate Bernstein polynomials).

  ngsxfem does not link against LAPACK(E) itself (on some platforms NGSolve
  provides LAPACK only internally). This header provides self-contained
  implementations of the two routines for exactly the arguments Algoim uses
  (square matrices, row-major layout, no eigenvectors). The directory of this
  header is put first on the include path of cutint/algoim_wrapper.cpp, the
  only translation unit that includes Algoim.

  The functions live in a C++ namespace (and are made visible by
  using-declarations) so that they can never be confused with the C symbols of
  a "real" LAPACKE library that may be loaded in the same process.

  The implementations are in cutint/algoim_lapacke.cpp.
*/

/// marker that allows cutint/algoim_wrapper.cpp to check that this header (and
/// not a lapacke.h of the system) has been picked up
#define NGSXFEM_ALGOIM_LAPACKE 1

#ifndef LAPACK_ROW_MAJOR
#define LAPACK_ROW_MAJOR 101
#endif
#ifndef LAPACK_COL_MAJOR
#define LAPACK_COL_MAJOR 102
#endif

namespace xintegration
{
  namespace algoim_lapacke
  {
    typedef int lapack_int;

    /// Singular value decomposition A = U diag(s) Vt of a square n x n matrix
    /// (only m == n, jobu == jobvt == 'A' is supported). Singular values are
    /// returned in descending order. On exit, A is overwritten.
    /// Returns 0 on success, -1 for unsupported arguments and > 0 if the
    /// iteration did not converge.
    lapack_int LAPACKE_dgesvd (int matrix_layout, char jobu, char jobvt,
                               lapack_int m, lapack_int n, double * a, lapack_int lda,
                               double * s, double * u, lapack_int ldu,
                               double * vt, lapack_int ldvt, double * superb);

    /// Generalized eigenvalues lambda = (alphar + i alphai) / beta of the
    /// matrix pencil (A,B) (only jobvl == jobvr == 'N', sense == 'N' is
    /// supported, balancing is not performed). Infinite eigenvalues are
    /// returned with beta == 0. ilo, ihi, lscale, rscale, abnrm, bbnrm are set
    /// to trivial values. On exit, A and B are overwritten.
    /// Returns 0 on success, -1 for unsupported arguments and > 0 if the
    /// iteration did not converge.
    lapack_int LAPACKE_dggevx (int matrix_layout, char balanc, char jobvl, char jobvr, char sense,
                               lapack_int n, double * a, lapack_int lda, double * b, lapack_int ldb,
                               double * alphar, double * alphai, double * beta,
                               double * vl, lapack_int ldvl, double * vr, lapack_int ldvr,
                               lapack_int * ilo, lapack_int * ihi,
                               double * lscale, double * rscale, double * abnrm, double * bbnrm,
                               double * rconde, double * rcondv);

    /// Eigenvalues wr + i wi of a general real n x n matrix (row-major,
    /// leading dimension lda); the matrix is overwritten. Hessenberg reduction
    /// followed by the Francis double-shift QR algorithm.
    /// Returns 0 on success and > 0 if the iteration did not converge.
    int RealEigenvalues (int n, double * a, int lda, double * wr, double * wi);
  }
}

using xintegration::algoim_lapacke::lapack_int;
using xintegration::algoim_lapacke::LAPACKE_dgesvd;
using xintegration::algoim_lapacke::LAPACKE_dggevx;
