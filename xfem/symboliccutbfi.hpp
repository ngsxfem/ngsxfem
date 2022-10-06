#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
using namespace xintegration;

// #include "xfiniteelement.hpp"
// #include "../spacetime/spacetimefe.hpp"
// #include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{

  class SymbolicCutBilinearFormIntegrator : public SymbolicBilinearFormIntegrator
  {
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;    
  public:

    SymbolicCutBilinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                       shared_ptr<CoefficientFunction> acf,
                                       VorB vb = VOL,
                                       VorB element_vb = VOL);
    
    virtual xbool IsSymmetric() const { return maybe; }
    virtual string Name () const { return string ("Symbolic Cut BFI"); }


    
    virtual void 
    CalcElementMatrix     (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<double> elmat,
                          LocalHeap & lh) const;

    virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<double> elmat,
                          bool & symmetric_so_far,                          
                          LocalHeap & lh) const;

    virtual void 
    CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<Complex> elmat,
                          bool & symmetric_so_far,                          
                          LocalHeap & lh) const;

    
    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixAdd (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
                                 FlatMatrix<SCAL_RES> elmat,
                                 LocalHeap & lh) const;

    template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
    void T_CalcElementMatrixEBAdd (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatMatrix<SCAL_RES> elmat,
                                   LocalHeap & lh) const;

    virtual void 
    CalcLinearizedElementMatrix (const FiniteElement & fel,
                                 const ElementTransformation & trafo, 
				 FlatVector<double> elveclin,
                                 FlatMatrix<double> elmat,
                                 LocalHeap & lh) const;

    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                          const ElementTransformation & trafo, 
                                          FlatVector<SCAL> elveclin,
                                          FlatMatrix<SCAL> elmat,
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
			LocalHeap & lh) const;

      
    template <int D, typename SCAL, typename SCAL_SHAPES>
    void T_ApplyElementMatrixEB (const FiniteElement & fel, 
                                 const ElementTransformation & trafo, 
                                 const FlatVector<SCAL> elx, 
                                 FlatVector<SCAL> ely,
                                 void * precomputed,
                                 LocalHeap & lh) const
    {
      throw Exception("SymbolicCutBilinearFormIntegrator::T_ApplyElementMatrixEB not yet implemented");
    }

  };
  
  class SymbolicCutFacetBilinearFormIntegrator : public SymbolicFacetBilinearFormIntegrator
  {
  protected:
    int time_order = -1;
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;    
    
  public:
    SymbolicCutFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                            shared_ptr<CoefficientFunction> acf,
                                            DOMAIN_TYPE adt,
                                            int asubdivlvl);
    SymbolicCutFacetBilinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                            shared_ptr<CoefficientFunction> acf);

    virtual VorB VB () const { return vb; }
    virtual xbool IsSymmetric() const { return maybe; }
    void SetTimeIntegrationOrder(int tiorder) { time_order = tiorder; }
    
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(neighbor_testfunction,
                                                                          element_boundary); }

    template <typename SCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix(const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<SCAL> elmat,
                     LocalHeap & lh) const;

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }


    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }

    
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
  class SymbolicFacetBilinearFormIntegrator2 : public SymbolicFacetBilinearFormIntegrator
  {
  // uncut integrator for facets with space-time capabilities 
  // (may become redundant if SymbolicCutFacetBLFI is fully implemented)
  // use case: "traditional" ghost penalty (higher order derivatives) 
  protected:
    int time_order = -1;
  public:
    SymbolicFacetBilinearFormIntegrator2 (shared_ptr<CoefficientFunction> acf);
    void SetTimeIntegrationOrder(int tiorder) { time_order = tiorder; }

    virtual VorB VB () const { return vb; }
    virtual xbool IsSymmetric() const { return maybe; }
    
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(neighbor_testfunction,
                                                                          element_boundary); }

    template <typename SCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix(const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<SCAL> elmat,
                     LocalHeap & lh) const;

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }


    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans,  
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
    {
      throw Exception("SymbolicFacetBilinearFormIntegrator2::CalcFacetMatrix on boundary not yet implemented");
    }

    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      throw Exception("SymbolicFacetBilinearFormIntegrator2::ApplyFacetMatrix not yet implemented");
    }


    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      throw Exception("SymbolicFacetBilinearFormIntegrator2::ApplyFacetMatrix not yet implemented");
    }

  };


  class SymbolicFacetPatchBilinearFormIntegrator : public SymbolicFacetBilinearFormIntegrator
  {
  protected:
    int time_order = -1;
    bool has_tref = false;
    double tref = 0.0;
  public:
    SymbolicFacetPatchBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf);
    void SetTimeIntegrationOrder(int tiorder) { time_order = tiorder; }
    void SetReferenceTime(double _tref) { has_tref = true; tref = _tref; }

    virtual VorB VB () const { return vb; }
    virtual xbool IsSymmetric() const { return maybe; }  // correct would be: don't know
    
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(neighbor_testfunction,
                                                                          element_boundary); }
    
    template <typename SCAL, typename SCAL_SHAPES = double>
    void T_CalcFacetMatrix(const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<SCAL> elmat,
                     LocalHeap & lh) const;

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
                     {
                         T_CalcFacetMatrix(volumefel1,LocalFacetNr1,eltrans1,ElVertices1,volumefel2,LocalFacetNr2,eltrans2,ElVertices2,elmat,lh);
                     }

    virtual void
    CalcLinearizedFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                     const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                     const FiniteElement & volumefel2, int LocalFacetNr2,
                     const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                     FlatVector<double> elvec,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
    {
      static bool warned = false;
      if (!warned)
      {
        cout << IM(3) << "WARNING: SymbolicFacetPatchBilinearFormIntegrator::CalcLinearizedFacetMatrix called. The form is treated as bilinear and hence CalcFacetMatri is called.";
        warned = true;
      }
      CalcFacetMatrix(volumefel1, LocalFacetNr1, eltrans1, ElVertices1, volumefel2, LocalFacetNr2, eltrans2, ElVertices2, elmat, lh);
    }
    

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                     const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                     const ElementTransformation & seltrans, FlatArray<int> & SElVertices, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
    {
      throw Exception("SymbolicFacetPatchBilinearFormIntegrator::CalcFacetMatrix on boundary not yet implemented");
    }

    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      static bool warned = false;
      if (!warned)
      {
        cout << IM(3) << "WARNING: SymbolicFacetPatchBilinearFormIntegrator::ApplyFacetMatrix called. The application is done through the computation of the element matrices (i.e. slower than possible).";
        warned = true;
      }
      FlatMatrix<double> mat(ely.Size(), elx.Size(), lh);
      CalcFacetMatrix (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                       volumefel2, LocalFacetNr2, eltrans2, ElVertices2, mat, lh);
      ely = mat * elx;
    }

    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                      const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                      const FiniteElement & volumefel2, int LocalFacetNr2,
                      const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                      FlatVector<Complex> elx, FlatVector<Complex> ely,
                      LocalHeap & lh) const
    {
      static bool warned = false;
      if (!warned)
      {
        cout << IM(3) << "WARNING: SymbolicFacetPatchBilinearFormIntegrator::ApplyFacetMatrix called. The application is done through the computation of the element matrices (i.e. slower than possible).";
        warned = true;
      }
      FlatMatrix<Complex> mat(ely.Size(), elx.Size(), lh);
      CalcFacetMatrix (volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                       volumefel2, LocalFacetNr2, eltrans2, ElVertices2, mat, lh);
      ely = mat * elx;
    }


    virtual void
    ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                      const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                      const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                      FlatVector<double> elx, FlatVector<double> ely,
                      LocalHeap & lh) const
    {
      static bool warned = false;
      if (!warned)
      {
        cout << IM(3) << "WARNING: SymbolicFacetPatchBilinearFormIntegrator::ApplyFacetMatrix called. The application is done through the computation of the element matrices (i.e. slower than possible).";
        warned = true;
      }
      FlatMatrix<double> mat(ely.Size(), elx.Size(), lh);
      CalcFacetMatrix (volumefel, LocalFacetNr,
                       eltrans, ElVertices,
                       seltrans, SElVertices, mat, lh);
      ely = mat * elx;
    }

  };


}
