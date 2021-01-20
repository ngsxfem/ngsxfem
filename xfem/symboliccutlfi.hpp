#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
using namespace xintegration;

namespace ngfem
{


  class SymbolicCutLinearFormIntegrator : public SymbolicLinearFormIntegrator
  {
    LevelsetIntegrationDomain lsetintdom;    
  public:
    SymbolicCutLinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                     shared_ptr<CoefficientFunction> acf,
                                     VorB vb = VOL);
    
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
