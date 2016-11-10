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
    
    SymbolicCutBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                       shared_ptr<CoefficientFunction> acf,
                                       DOMAIN_TYPE adt,
                                       int aforce_intorder = -1,
                                       int asubdivlvl = 0);

    virtual bool BoundaryForm() const { return false; }
    virtual bool IsSymmetric() const { return true; }  // correct would be: don't know
    virtual string Name () const { return string ("Symbolic Cut BFI"); }

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const;

    virtual void 
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & trafo, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const;    

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



  class SymbolicCutFacetBilinearFormIntegrator : public SymbolicFacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> cf_lset;
    DOMAIN_TYPE dt = NEG;
    int force_intorder = -1;
    int subdivlvl = 0;
  public:
    SymbolicCutFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                            shared_ptr<CoefficientFunction> acf,
                                            DOMAIN_TYPE adt,
                                            int aforce_intorder,
                                            int asubdivlvl);

    virtual bool BoundaryForm() const { return vb == BND; }
    virtual bool IsSymmetric() const { return true; }  // correct would be: don't know
    
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(neighbor_testfunction,
                                                                          element_boundary); }
    
    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<double> & elmat,
                     LocalHeap & lh) const;

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans,  
                     FlatMatrix<double> & elmat,
                     LocalHeap & lh) const
    {
      throw Exception("SymbolicCutFacetBilinearFormIntegrator::CalcFacetMatrix on boundary not yet implemented");
    }

    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      throw Exception("SymbolicCutFacetBilinearFormIntegrator::ApplyFacetMatrix not yet implemented");
    }


    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      throw Exception("SymbolicCutFacetBilinearFormIntegrator::ApplyFacetMatrix not yet implemented");
    }

  };


}
