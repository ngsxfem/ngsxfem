#ifndef FILE_ERROR_HPP
#define FILE_ERROR_HPP

#include <ngstd.hpp> // for Array
#include <fem.hpp>   // for ScalarFiniteElement
#include <solve.hpp>

// #include "xintegration.hpp"
#include "../spacetime/spacetimefespace.hpp"
#include "../xfem/xfemIntegrators.hpp"
#include "../xfem/stxfemIntegrators.hpp"
#include "../xfem/setvaluesx.hpp"
#include "../utils/error.hpp"

using namespace ngsolve;

/// from ngsolve

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"

namespace ngfem
{
  template<int D>
  void CalcDxShapeOfCoeff(const CoefficientFunction * coef, const MappedIntegrationPoint<D,D>& mip, 
                          double time,
                          Vec<D>& der, LocalHeap& lh);

  template <int D>
  class SolutionCoefficients
  {
  protected:
    shared_ptr<CoefficientFunction> conv_n = NULL;
    shared_ptr<CoefficientFunction> conv_p = NULL;
    shared_ptr<CoefficientFunction> coef_n = NULL;
    shared_ptr<CoefficientFunction> coef_p = NULL;
    shared_ptr<CoefficientFunction> coef_d_n = NULL;
    shared_ptr<CoefficientFunction> coef_d_p = NULL;
    shared_ptr<CoefficientFunction> coef_jumprhs = NULL;
    shared_ptr<CoefficientFunction> lset = NULL;

  public:
    // SolutionCoefficients(const CoefficientFunction * a_coef_n,
    //                      const CoefficientFunction * a_coef_p,
    //                      const CoefficientFunction * a_coef_d_n,
    //                      const CoefficientFunction * a_coef_d_p,
    //                      const CoefficientFunction * a_lset);
    SolutionCoefficients(shared_ptr<PDE> pde, const Flags & flags);
    ~SolutionCoefficients();

    bool HasConvectionNeg() { return ! (conv_n == NULL); }
    bool HasConvectionPos() { return ! (conv_p == NULL); }
    bool HasSolutionNeg() { return ! (coef_n == NULL); }
    bool HasSolutionPos() { return ! (coef_p == NULL); }
    bool HasSolutionDNeg() { return ! (coef_d_n == NULL); }
    bool HasSolutionDPos() { return ! (coef_d_p == NULL); }
    bool HasLevelSet() { return ! (lset == NULL); }
    bool HasJumpRhs() { return ! (coef_jumprhs == NULL); }

    const CoefficientFunction & GetConvectionNeg() { return *conv_n; }
    const CoefficientFunction & GetConvectionPos() { return *conv_p; }
    const CoefficientFunction & GetSolutionNeg() { return *coef_n; }
    const CoefficientFunction & GetSolutionPos() { return *coef_p; }
    const CoefficientFunction & GetSolutionDNeg() { return *coef_d_n; }
    const CoefficientFunction & GetSolutionDPos() { return *coef_d_p; }
    const CoefficientFunction & GetLevelSet() { return *lset; }
    const CoefficientFunction & GetJumpRhs() { return *coef_jumprhs; }

  };

  class ErrorTable
  {
  public:
    Array<double> l2err_p;
    Array<double> l2err_n;
    Array<double> l2err;
    Array<double> h1err_p;
    Array<double> h1err_n;
    Array<double> h1err;
    Array<double> iferr;  
    Array<double> ifdudnerr;  
    Array<double> ifsigmanerr;  

    ErrorTable();

    void Reset()
    {
      l2err_p.SetSize(0);
      l2err_n.SetSize(0);
      l2err.SetSize(0);  
      h1err_p.SetSize(0);
      h1err_n.SetSize(0);
      h1err.SetSize(0);  
      iferr.SetSize(0);  
      ifdudnerr.SetSize(0);  
      ifsigmanerr.SetSize(0);  
    }
  };

  template<int D>
  void CalcXError (shared_ptr<GridFunction> gfu, 
                   shared_ptr<GridFunction> gfu2, 
                   SolutionCoefficients<D> & solcoef, 
                   int intorder, 
                   double a_neg, 
                   double a_pos, 
                   double b_neg, 
                   double b_pos, 
                   double time, 
                   ErrorTable & errtab, 
                   LocalHeap & lh,
                   bool output,
                   const Flags & flags);
  
  template<int D>
  void CalcTraceDiff (shared_ptr<GridFunction> gfu,
                      shared_ptr<CoefficientFunction> coef,
                      int intorder, Array<double>& errors,
                      LocalHeap & lh);

}
#endif
