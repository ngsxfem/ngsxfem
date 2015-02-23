#ifndef FILE_XMEANCURV_HPP
#define FILE_XMEANCURV_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"
#include "../utils/stcoeff.hpp"

namespace ngfem
{
  template <int D>
  class XLBMeanCurvIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    XLBMeanCurvIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0])
    {
      // if (D==3)
      //   throw Exception("Implementation only 2D for now");
    }
    virtual ~XLBMeanCurvIntegrator(){ ; };

    virtual string Name () const { return "XLBMeanCurvIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> elvec,
                       LocalHeap & lh) const;
  };

}

#endif

