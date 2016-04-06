#ifndef FILE_HDGTRACEINTEGRATORS_HPP
#define FILE_HDGTRACEINTEGRATORS_HPP

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
  class HDGTraceLaplaceBeltramiIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_alpha;
    double param_IP_edge;
    shared_ptr<CoefficientFunction> coef_param_normaldiff;
    double param_IP_facet;
  public:
    HDGTraceLaplaceBeltramiIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_alpha(coeffs[0]),
        param_IP_edge(coeffs[1] -> EvaluateConst()),
        coef_param_normaldiff(coeffs[2]),
        param_IP_facet(coeffs[3] -> EvaluateConst())
    {;}
    virtual ~HDGTraceLaplaceBeltramiIntegrator(){ ; };
    virtual string Name () const { return "HDGTraceLaplaceBeltramiIntegrator"; }
    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    virtual bool BoundaryForm () const { return false; }
    virtual bool IsSymmetric () const { return true; }
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                         const ElementTransformation & eltrans,
                         FlatMatrix<double> elmat,
                         LocalHeap & lh) const;
  };

}

#endif
