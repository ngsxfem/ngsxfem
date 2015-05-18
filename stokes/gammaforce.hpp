#ifndef FILE_GAMMAFORCE_HPP
#define FILE_GAMMAFORCE_HPP

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
  class GammaForceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
    shared_ptr<CoefficientFunction> coef_lset;
  public:
    GammaForceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0])
    {
      ;
    }
    virtual ~GammaForceIntegrator(){ ; };

    virtual string Name () const { return "GammaForceIntegrator"; }

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

