#pragma once

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
#include "calcpointshift.hpp"

namespace ngfem
{

  template <int D>
  class ShiftIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_lset_p1;
    shared_ptr<CoefficientFunction> coef_lset_ho;
    shared_ptr<CoefficientFunction> coef_blending = nullptr;
    double max_deform = -1;
    double lower_lset_bound = 0.0;
    double upper_lset_bound = 0.0;
    shared_ptr<CoefficientFunction> qn;
  public:
    ShiftIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs);
    virtual ~ShiftIntegrator(){; };
    virtual string Name () const { return "ShiftIntegrator"; }
    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // virtual bool BoundaryForm () const { return false; }
    virtual VorB VB () const { return VOL; }
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans,
                                    FlatVector<double> elvec,
                                    LocalHeap & lh,
                                    shared_ptr<LsetEvaluator<D>> lseteval) const;
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans,
                                    FlatVector<double> elvec,
                                    LocalHeap & lh) const
    {
      CalcElementVector(fel,eltrans,elvec,lh,nullptr);
    }
  };

}
