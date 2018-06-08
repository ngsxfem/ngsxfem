#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
using namespace xintegration;

namespace ngfem
{


  class SymbolicCutLinearFormIntegrator : public SymbolicLinearFormIntegrator
  {
    shared_ptr<GridFunction> gf_lset = nullptr;
    shared_ptr<CoefficientFunction> cf_lset;
    DOMAIN_TYPE dt = NEG;
    int force_intorder = -1;
    int subdivlvl = 0;
    int time_order = -1;
    SWAP_DIMENSIONS_POLICY pol;

  public:

    SymbolicCutLinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                     shared_ptr<CoefficientFunction> acf,
                                     DOMAIN_TYPE adt,
                                     int aforce_intorder = -1,
                                     int asubdivlvl = 0,
                                     SWAP_DIMENSIONS_POLICY apol = FIND_OPTIMAL,
                                     VorB vb = VOL);

    void SetTimeIntegrationOrder(int tiorder) { time_order = tiorder; }
    virtual VorB VB () const { return VOL; }
    virtual string Name () const { return string ("Symbolic Cut LFI"); }


    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;
      
    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const;

    template <typename SCAL> 
    void T_CalcElementVector (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatVector<SCAL> elvec,
                              LocalHeap & lh) const;
  };
  
}
