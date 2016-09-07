#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xdecompose.hpp"
using namespace xintegration;

// #include "xfiniteelement.hpp"
// #include "../spacetime/spacetimefe.hpp"
// #include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{

  class SymbolicCutBilinearFormIntegrator : public SymbolicBilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> cf_lset;
    DOMAIN_TYPE dt = NEG;
    int force_intorder = -1;
    int subdivlvl = 0;
  public:
    
    SymbolicCutBilinearFormIntegrator (shared_ptr<CoefficientFunction> & cf_lset,
                                       shared_ptr<CoefficientFunction> & acf,
                                       DOMAIN_TYPE adt,
                                       int aforce_intorder = -1,
                                       int asubdivlvl = 0);

    virtual bool BoundaryForm() const { return false; }
    virtual bool IsSymmetric() const { return true; }  // correct would be: don't know
    virtual string Name () const { return string ("Symbolic Cut BFI"); }

    template <typename SCAL, typename SCAL_SHAPES = double>
    void T_CalcElementMatrix (const FiniteElement & fel,
                              const ElementTransformation & trafo, 
                              FlatMatrix<SCAL> elmat,
                              LocalHeap & lh) const;

    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_CalcElementMatrixEB (const FiniteElement & fel,
                                const ElementTransformation & trafo, 
                                FlatMatrix<SCAL> elmat,
                                LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::T_CalcElementMatrixEB not yet implemented");
    }

    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::CalcLinearizedElementMatrix not yet implemented");
    }

    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                          const ElementTransformation & trafo, 
                                          FlatVector<double> elveclin,
                                          FlatMatrix<double> elmat,
                                          LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::T_CalcLinearizedElementMatrixEB not yet implemented");
    }
    
    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & trafo, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::ApplyElementMatrix not yet implemented");
    }

      
    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_ApplyElementMatrixEB (const FiniteElement & fel, 
                                 const ElementTransformation & trafo, 
                                 const FlatVector<double> elx, 
                                 FlatVector<double> ely,
                                 void * precomputed,
                                 LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::T_ApplyElementMatrixEB not yet implemented");
    }

  };

}
