/*
  Self-contained replacements of the two LAPACKE routines used by the bundled
  Algoim library, see algoim_lapacke/lapacke.h.

  * LAPACKE_dgesvd: wrapper of NGSolve's (LAPACK-free) SVD ngbla::CalcSVD.
    Algoim decomposes Bernstein-Vandermonde matrices of size degree+1, where
    the degree of the resultants that Algoim computes in 3D can be several
    hundred, so that an O(n^3) algorithm with a small constant is mandatory.
    (The decompositions are cached by Algoim.)
  * LAPACKE_dggevx: generalized eigenvalues of a pencil (A,B) by a
    shift-and-invert transformation to the standard eigenvalue problem of
    M = (A - sigma B)^{-1} B (eigenvalues mu = 1/(lambda - sigma), mu = 0 for
    infinite lambda), followed by a reduction of M to Hessenberg form and the
    Francis double-shift QR algorithm (cf. G. H. Golub, C. F. Van Loan,
    Matrix Computations, Sec. 7.4/7.5). Algoim uses the eigenvalues only to
    find real roots in [0,1] of univariate Bernstein polynomials, so the shifts
    are chosen close to (but not exactly at) this interval.

  This file must not include Algoim headers.
*/

#include <bla.hpp>   // NGSolve's dense linear algebra (CalcSVD)
#include "algoim_lapacke/lapacke.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace xintegration
{
  namespace algoim_lapacke
  {
    namespace
    {
      const double eps = std::numeric_limits<double>::epsilon();

      // Householder vector v (length m <= 3) and factor beta with
      // (I - beta v v^T) x = alpha e_0; returns false if x == 0
      bool Householder (int m, const double * x, double * v, double & beta)
      {
        double nrm = 0;
        for (int i = 0; i < m; i++)
          nrm += x[i]*x[i];
        nrm = std::sqrt(nrm);
        if (nrm == 0)
          return false;
        double alpha = (x[0] > 0) ? -nrm : nrm;
        for (int i = 0; i < m; i++)
          v[i] = x[i];
        v[0] -= alpha;
        double vtv = 0;
        for (int i = 0; i < m; i++)
          vtv += v[i]*v[i];
        if (vtv == 0)
          return false;
        beta = 2.0 / vtv;
        return true;
      }

      // eigenvalues of the 2x2 matrix [[a,b],[c,d]]
      void Eigenvalues2x2 (double a, double b, double c, double d,
                           double & wr1, double & wi1, double & wr2, double & wi2)
      {
        double p = 0.5 * (a - d);
        double disc = p*p + b*c;
        if (disc >= 0)
        {
          double z = p + std::copysign(std::sqrt(disc), p);
          wr1 = d + z;
          wr2 = (z != 0) ? d - b*c/z : d;
          wi1 = wi2 = 0.0;
        }
        else
        {
          wr1 = wr2 = d + p;
          wi1 = std::sqrt(-disc);
          wi2 = -wi1;
        }
      }
    }

    int RealEigenvalues (int n, double * a, int lda, double * wr, double * wi)
    {
      auto H = [a,lda] (int i, int j) -> double & { return a[i*lda+j]; };

      // 1. reduction to upper Hessenberg form by Householder reflections
      std::vector<double> x(n), v(n);
      for (int k = 0; k < n-2; k++)
      {
        int m = n-k-1;
        double nrm = 0;
        for (int i = 0; i < m; i++)
        {
          x[i] = H(k+1+i,k);
          nrm += x[i]*x[i];
        }
        nrm = std::sqrt(nrm);
        if (nrm == 0)
          continue;
        double alpha = (x[0] > 0) ? -nrm : nrm;
        for (int i = 0; i < m; i++)
          v[i] = x[i];
        v[0] -= alpha;
        double vtv = 0;
        for (int i = 0; i < m; i++)
          vtv += v[i]*v[i];
        if (vtv == 0)
          continue;
        double beta = 2.0 / vtv;
        // from the left: rows k+1..n-1, columns k..n-1
        for (int j = k; j < n; j++)
        {
          double s = 0;
          for (int i = 0; i < m; i++)
            s += v[i]*H(k+1+i,j);
          s *= beta;
          for (int i = 0; i < m; i++)
            H(k+1+i,j) -= s*v[i];
        }
        // from the right: rows 0..n-1, columns k+1..n-1
        for (int i = 0; i < n; i++)
        {
          double s = 0;
          for (int j = 0; j < m; j++)
            s += H(i,k+1+j)*v[j];
          s *= beta;
          for (int j = 0; j < m; j++)
            H(i,k+1+j) -= s*v[j];
        }
        for (int i = k+2; i < n; i++)
          H(i,k) = 0.0;
      }

      double hnorm = 0;
      for (int i = 0; i < n; i++)
        for (int j = std::max(0,i-1); j < n; j++)
          hnorm += std::fabs(H(i,j));

      // 2. Francis double-shift QR iteration on the active block [lo,hi]
      int hi = n-1;
      int iter = 0;
      const int maxiter = 30 * std::max(n,1);
      int totaliter = 0;
      while (hi >= 0)
      {
        // look for a negligible subdiagonal entry
        int lo = hi;
        while (lo > 0)
        {
          // relative criterion (neighbouring diagonal entries) and absolute
          // criterion w.r.t. the norm of H (both backward stable); the latter
          // is needed if the neighbouring diagonal entries are (nearly) zero
          double s = std::fabs(H(lo-1,lo-1)) + std::fabs(H(lo,lo));
          if (s == 0)
            s = hnorm;
          double sub = std::fabs(H(lo,lo-1));
          if (sub <= eps * s || sub <= 0.5 * eps * hnorm)
          {
            H(lo,lo-1) = 0.0;
            break;
          }
          lo--;
        }

        if (lo == hi) // 1x1 block converged
        {
          wr[hi] = H(hi,hi);
          wi[hi] = 0.0;
          hi--;
          iter = 0;
          continue;
        }
        if (lo == hi-1) // 2x2 block converged
        {
          Eigenvalues2x2(H(hi-1,hi-1), H(hi-1,hi), H(hi,hi-1), H(hi,hi),
                         wr[hi-1], wi[hi-1], wr[hi], wi[hi]);
          hi -= 2;
          iter = 0;
          continue;
        }

        if (++totaliter > maxiter)
          return hi+1;
        iter++;

        // shifts: eigenvalues of the trailing 2x2 block (exceptional shifts
        // in case of slow convergence)
        double tr, det;
        if (iter % 10 == 0)
        {
          double s = std::fabs(H(hi,hi-1)) + std::fabs(H(hi-1,hi-2));
          double h11 = 0.75*s + H(hi,hi);
          double h12 = -0.4375*s;
          tr = 2*h11;
          det = h11*h11 - h12*s;
        }
        else
        {
          tr = H(hi-1,hi-1) + H(hi,hi);
          det = H(hi-1,hi-1)*H(hi,hi) - H(hi-1,hi)*H(hi,hi-1);
        }

        // first column of (H - s1 I)(H - s2 I)
        double xyz[3];
        xyz[0] = H(lo,lo)*H(lo,lo) + H(lo,lo+1)*H(lo+1,lo) - tr*H(lo,lo) + det;
        xyz[1] = H(lo+1,lo) * (H(lo,lo) + H(lo+1,lo+1) - tr);
        xyz[2] = H(lo+1,lo) * H(lo+2,lo+1);

        // chase the bulge
        for (int k = lo; k <= hi-2; k++)
        {
          double hv[3], beta;
          if (Householder(3, xyz, hv, beta))
          {
            int jstart = std::max(lo, k-1);
            for (int j = jstart; j <= hi; j++)
            {
              double s = beta * (hv[0]*H(k,j) + hv[1]*H(k+1,j) + hv[2]*H(k+2,j));
              H(k,j) -= s*hv[0];
              H(k+1,j) -= s*hv[1];
              H(k+2,j) -= s*hv[2];
            }
            int iend = std::min(k+3, hi);
            for (int i = lo; i <= iend; i++)
            {
              double s = beta * (H(i,k)*hv[0] + H(i,k+1)*hv[1] + H(i,k+2)*hv[2]);
              H(i,k) -= s*hv[0];
              H(i,k+1) -= s*hv[1];
              H(i,k+2) -= s*hv[2];
            }
          }
          xyz[0] = H(k+1,k);
          xyz[1] = H(k+2,k);
          xyz[2] = (k < hi-2) ? H(k+3,k) : 0.0;
        }
        // final 2x2 reflection
        {
          double hv[2], beta;
          int k = hi-1;
          if (Householder(2, xyz, hv, beta))
          {
            for (int j = std::max(lo, k-1); j <= hi; j++)
            {
              double s = beta * (hv[0]*H(k,j) + hv[1]*H(k+1,j));
              H(k,j) -= s*hv[0];
              H(k+1,j) -= s*hv[1];
            }
            for (int i = lo; i <= hi; i++)
            {
              double s = beta * (H(i,k)*hv[0] + H(i,k+1)*hv[1]);
              H(i,k) -= s*hv[0];
              H(i,k+1) -= s*hv[1];
            }
          }
        }
        // clean up the bulge remainders (exact zeros in exact arithmetic)
        for (int i = lo+2; i <= hi; i++)
        {
          H(i,i-2) = 0.0;
          if (i > lo+2)
            H(i,i-3) = 0.0;
        }
      }
      return 0;
    }


    lapack_int LAPACKE_dgesvd (int matrix_layout, char jobu, char jobvt,
                               lapack_int m, lapack_int n, double * a, lapack_int lda,
                               double * s, double * u, lapack_int ldu,
                               double * vt, lapack_int ldvt, double * superb)
    {
      if (matrix_layout != LAPACK_ROW_MAJOR || jobu != 'A' || jobvt != 'A' || m != n || n < 0)
        return -1;
      if (n == 0)
        return 0;

      // NGSolve's SVD (Householder bidiagonalization + QR iteration):
      // A = U diag(A) V^T, A is overwritten by the singular values (in descending order)
      ngbla::Matrix<double> A(n, n);
      ngbla::Matrix<double, ngbla::ColMajor> U(n, n), V(n, n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          A(i,j) = a[i*lda+j];
      ngbla::CalcSVD(A, U, V);

      // LAPACK conventions: non-negative singular values in descending order
      std::vector<int> perm(n);
      for (int j = 0; j < n; j++)
      {
        perm[j] = j;
        if (A(j,j) < 0)
        {
          A(j,j) = -A(j,j);
          for (int i = 0; i < n; i++)
            U(i,j) = -U(i,j);
        }
      }
      std::stable_sort(perm.begin(), perm.end(), [&A] (int i, int j) { return A(i,i) > A(j,j); });
      for (int jj = 0; jj < n; jj++)
      {
        int j = perm[jj];
        s[jj] = A(j,j);
        for (int i = 0; i < n; i++)
        {
          u[i*ldu+jj] = U(i,j);
          vt[jj*ldvt+i] = V(i,j);
        }
      }
      for (int i = 0; i < n-1; i++)
        superb[i] = 0.0;
      return 0;
    }


    lapack_int LAPACKE_dggevx (int matrix_layout, char balanc, char jobvl, char jobvr, char sense,
                               lapack_int n, double * a, lapack_int lda, double * b, lapack_int ldb,
                               double * alphar, double * alphai, double * beta,
                               double * vl, lapack_int ldvl, double * vr, lapack_int ldvr,
                               lapack_int * ilo, lapack_int * ihi,
                               double * lscale, double * rscale, double * abnrm, double * bbnrm,
                               double * rconde, double * rcondv)
    {
      if (matrix_layout != LAPACK_ROW_MAJOR || jobvl != 'N' || jobvr != 'N' || sense != 'N' || n < 0)
        return -1;
      if (n == 0)
        return 0;

      if (ilo) *ilo = 1;
      if (ihi) *ihi = n;
      for (int i = 0; i < n; i++)
      {
        if (lscale) lscale[i] = 1.0;
        if (rscale) rscale[i] = 1.0;
      }

      double anrm = 0, bnrm = 0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
          anrm = std::max(anrm, std::fabs(a[i*lda+j]));
          bnrm = std::max(bnrm, std::fabs(b[i*ldb+j]));
        }
      if (abnrm) *abnrm = anrm;
      if (bbnrm) *bbnrm = bnrm;

      if (bnrm == 0) // all eigenvalues are infinite (or the pencil is singular)
      {
        for (int i = 0; i < n; i++)
        {
          alphar[i] = 1.0; alphai[i] = 0.0; beta[i] = 0.0;
        }
        return 0;
      }

      // shift-and-invert: M = (A - sigma B)^{-1} B; the shifts are close to [0,1]
      // (where Algoim looks for roots) and are chosen such that A - sigma B is
      // well-conditioned (in the sense of the LU pivots)
      const double shifts[] = { 0.31830988618379067, 0.70710678118654752, -0.27182818284590452,
                                1.41421356237309505, 0.57721566490153286, -1.61803398874989485,
                                2.23606797749978970 };
      std::vector<double> LU(n*n), M(n*n);
      std::vector<int> piv(n);
      double sigma = 0;
      bool found = false;
      double best_quality = -1;
      for (double shift : shifts)
      {
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            LU[i*n+j] = a[i*lda+j] - shift * b[i*ldb+j];
        double cnrm = 0;
        for (int k = 0; k < n*n; k++)
          cnrm = std::max(cnrm, std::fabs(LU[k]));
        // LU decomposition with partial pivoting
        double minpivot = std::numeric_limits<double>::max();
        for (int k = 0; k < n; k++)
        {
          int p = k;
          for (int i = k+1; i < n; i++)
            if (std::fabs(LU[i*n+k]) > std::fabs(LU[p*n+k]))
              p = i;
          piv[k] = p;
          if (p != k)
            for (int j = 0; j < n; j++)
              std::swap(LU[k*n+j], LU[p*n+j]);
          double pivot = LU[k*n+k];
          minpivot = std::min(minpivot, std::fabs(pivot));
          if (pivot == 0)
            break;
          for (int i = k+1; i < n; i++)
          {
            double f = LU[i*n+k] / pivot;
            LU[i*n+k] = f;
            for (int j = k+1; j < n; j++)
              LU[i*n+j] -= f * LU[k*n+j];
          }
        }
        double quality = (cnrm > 0) ? minpivot / cnrm : 0.0;
        if (quality > best_quality)
        {
          best_quality = quality;
          sigma = shift;
        }
        if (quality > 1e-8)
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        // recompute the factorization for the best shift found
        if (best_quality <= 0)
        {
          // singular pencil (det(A - lambda B) == 0 for all tested lambda): report
          // alpha = beta = 0 as LAPACK does (Algoim then finds no root)
          for (int i = 0; i < n; i++)
          {
            alphar[i] = 0.0; alphai[i] = 0.0; beta[i] = 0.0;
          }
          return 0;
        }
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            LU[i*n+j] = a[i*lda+j] - sigma * b[i*ldb+j];
        for (int k = 0; k < n; k++)
        {
          int p = k;
          for (int i = k+1; i < n; i++)
            if (std::fabs(LU[i*n+k]) > std::fabs(LU[p*n+k]))
              p = i;
          piv[k] = p;
          if (p != k)
            for (int j = 0; j < n; j++)
              std::swap(LU[k*n+j], LU[p*n+j]);
          double pivot = LU[k*n+k];
          for (int i = k+1; i < n; i++)
          {
            double f = LU[i*n+k] / pivot;
            LU[i*n+k] = f;
            for (int j = k+1; j < n; j++)
              LU[i*n+j] -= f * LU[k*n+j];
          }
        }
      }

      // M = LU^{-1} B (column by column)
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          M[i*n+j] = b[i*ldb+j];
      for (int k = 0; k < n; k++)
        if (piv[k] != k)
          for (int j = 0; j < n; j++)
            std::swap(M[k*n+j], M[piv[k]*n+j]);
      for (int j = 0; j < n; j++)
      {
        for (int i = 1; i < n; i++)
          for (int k = 0; k < i; k++)
            M[i*n+j] -= LU[i*n+k] * M[k*n+j];
        for (int i = n-1; i >= 0; i--)
        {
          for (int k = i+1; k < n; k++)
            M[i*n+j] -= LU[i*n+k] * M[k*n+j];
          M[i*n+j] /= LU[i*n+i];
        }
      }

      std::vector<double> mur(n), mui(n);
      int info = RealEigenvalues(n, M.data(), n, mur.data(), mui.data());
      if (info != 0)
        return info;

      double mnrm = 0;
      for (int i = 0; i < n; i++)
        mnrm = std::max(mnrm, std::hypot(mur[i], mui[i]));
      for (int i = 0; i < n; i++)
      {
        double absmu2 = mur[i]*mur[i] + mui[i]*mui[i];
        if (absmu2 == 0 || std::sqrt(absmu2) <= 10 * n * eps * mnrm)
        {
          // mu = 0 <=> lambda = infinity
          alphar[i] = 1.0; alphai[i] = 0.0; beta[i] = 0.0;
        }
        else
        {
          // lambda = sigma + 1/mu
          alphar[i] = sigma + mur[i] / absmu2;
          alphai[i] = (mui[i] == 0) ? 0.0 : -mui[i] / absmu2;
          beta[i] = 1.0;
        }
      }
      return 0;
    }
  }
}
