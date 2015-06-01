#ifndef FILE_XSTOKESNITSCHE_HPP
#define FILE_XSTOKESNITSCHE_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"
#include "xfemIntegrators.hpp"
#include "stxfemIntegrators.hpp"

namespace ngfem
{

  template <int D>
  class XStokesNitscheIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> alpha_neg;
    shared_ptr<CoefficientFunction> alpha_pos;
    shared_ptr<CoefficientFunction> lambda;   
    shared_ptr<XLaplaceIntegrator<D>> laplace;
  public:
    XStokesNitscheIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        lambda(coeffs.Size() > 2 ? coeffs[2]: NULL)
    { 
      laplace = make_shared<XLaplaceIntegrator<D>>(coeffs);
      ;
    }
    
    virtual ~XStokesNitscheIntegrator()
    { 
      ;
    }

    virtual string Name () const { return "XStokesNitscheIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;
 
  };
}
#endif

