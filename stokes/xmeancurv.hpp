#ifndef FILE_XMEANCURV_HPP
#define FILE_XMEANCURV_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
// #include "../spacetime/spacetimeintegrators.hpp"
#include "../utils/stcoeff.hpp"

namespace ngfem
{
  template <int D, bool improved>
  class XLBMeanCurvIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
    shared_ptr<CoefficientFunction> coef_lset;
  public:
    XLBMeanCurvIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0])
    {
      if (improved) coef_lset=coeffs[1]; else coef_lset = NULL;
    }
    virtual ~XLBMeanCurvIntegrator(){ ; };

    virtual string Name () const { return "XLBMeanCurvIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual VorB VB () const { return VOL; }

    // Calculates the element matrix
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> elvec,
                       LocalHeap & lh) const;
  };

}

#endif

