#pragma once
/*********************************************************************/
/* File:   wrapperfespace.hpp                                        */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   12. July 2018                                             */
/*********************************************************************/
#include <comp.hpp>
#include "compressedfespace.hpp"

namespace ngcomp
{

  class RestrictedFESpace : public CompressedFESpace
  {
    // shared_ptr<FESpace> space;
    // Array<DofId> comp2all;
    // Array<DofId> all2comp;
    // shared_ptr<BitArray> active_dofs = nullptr;
  protected:
    shared_ptr<BitArray> active_els = nullptr;

  public:
    RestrictedFESpace (shared_ptr<FESpace> bfes, shared_ptr<BitArray> _active_els = nullptr);
    virtual ~RestrictedFESpace () {};
    void Update() override;
    shared_ptr<FESpace> GetBaseSpace() const { return space; }

    void WrapDofs(Array<DofId> & dnums) const
    {
      CompressedFESpace::WrapDofs(dnums);
      // for (DofId & d : dnums)
      //   if (IsRegularDof (d))
      //     d = all2comp[d];
    }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    /// get dof-nrs of the element of certain coupling type
    virtual void GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const override;

    virtual void SetActiveElements(shared_ptr<BitArray> actels) 
    {
      active_els = actels;
    }

    shared_ptr<BitArray> GetActiveDofs() const { return active_dofs; }
    shared_ptr<BitArray> GetActiveElements() const { return active_els; }

    // a name for our new fe-space
    virtual string GetClassName () const override
    {
      return "RestrictedFESpace(" + space->GetClassName() + ")";
    }

    /// update element coloring
    void FinalizeUpdate() override
    {
      CompressedFESpace::FinalizeUpdate();
    }

    virtual ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override;

  };
}

namespace ngfem
{
  class RestrictedDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
  public:
    RestrictedDifferentialOperator (shared_ptr<DifferentialOperator> adiffop)
      : DifferentialOperator(adiffop->Dim(), adiffop->BlockDim(),
                             adiffop->VB(), adiffop->DiffOrder()), diffop(adiffop) { 
      SetDimensions(adiffop->Dimensions());
    }

    virtual ~RestrictedDifferentialOperator (){ ; };
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return diffop->UsedDofs(fel); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())
        return make_shared<RestrictedDifferentialOperator> (diffoptrace);
      else
        return nullptr;
    }
    
    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		BareSliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    
    
    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
    BareSliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const override;    

    virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override;
    
    virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<Complex>> mat) const override;
    
    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override;
    
    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;
    
    virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<SIMD<Complex>> flux) const override;
    
    virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override;

    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const override;

    virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const override;

    virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;

    virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const override;

    shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                               shared_ptr<CoefficientFunction> dir,
                                               bool Eulerian) const override;
  };


}
