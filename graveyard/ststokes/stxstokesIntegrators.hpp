#ifndef FILE_SPACETIMESTOKESINTEGRATORS_HPP
#define FILE_SPACETIMESTOKESINTEGRATORS_HPP

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
  class SpaceTimeXStokesIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
    double t1;
    double t0;
    double tau;
  public:
    SpaceTimeXStokesIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;

      if (D==3)
        throw Exception("Implementation only 2D for now");
    }
    virtual ~SpaceTimeXStokesIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXStokesIntegrator"; }

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

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

  };

}

#endif

