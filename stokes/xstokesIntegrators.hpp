#ifndef FILE_XSTOKESINTEGRATORS_HPP
#define FILE_XSTOKESINTEGRATORS_HPP

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
  class XStokesIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    XStokesIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    {
      // if (D==3)
      //   throw Exception("Implementation only 2D for now");
    }
    virtual ~XStokesIntegrator(){ ; };

    virtual string Name () const { return "XStokesIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }
    virtual bool IsSymmetric () const { return true; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;
  };

}

#endif

