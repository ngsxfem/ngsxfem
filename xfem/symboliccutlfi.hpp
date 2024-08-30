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


  class SymbolicCutFacetLinearFormIntegrator : public SymbolicFacetLinearFormIntegrator
  {
  protected:
    int time_order = -1;
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;    
  public:
    SymbolicCutFacetLinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                     shared_ptr<CoefficientFunction> acf,
                                     VorB vb = BND);
    
    virtual string Name () const { return string ("Symbolic Cut Facet LFI"); }

    virtual void
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatVector<double> elvec,
			 LocalHeap & lh) const
    {
	 throw Exception("SymbolicCutFacetLinearFormIntegrator::CalcFacetVector not yet implemented");   
    }

    virtual void 
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatVector<Complex> elvec,
			 LocalHeap & lh) const

    {
	 throw Exception("SymbolicCutFacetLinearFormIntegrator::CalcFacetVector not yet implemented");   
    }

    virtual void
    CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                    const ElementTransformation & seltrans,
                    FlatVector<double> elvec,
                    LocalHeap & lh) const
    {
      T_CalcFacetVector (volumefel, LocalFacetNr, eltrans, ElVertices, seltrans, elvec, lh);
    }

    virtual void
    CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                    const ElementTransformation & seltrans,
                    FlatVector<Complex> elvec,
                    LocalHeap & lh) const
    {
      T_CalcFacetVector (volumefel, LocalFacetNr, eltrans, ElVertices, seltrans, elvec, lh);
    }

    template<typename TSCAL>
    void T_CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr,
                            const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                            const ElementTransformation & strafo,
                            FlatVector<TSCAL> elvec,
                            LocalHeap & lh) const;

  };


}
