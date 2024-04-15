#include <comp.hpp>
#include "restrictedfespace.hpp"

namespace ngcomp
{

  RestrictedFESpace::RestrictedFESpace (shared_ptr<FESpace> bfes, shared_ptr<BitArray> _active_els)
      : CompressedFESpace (bfes), active_els(_active_els)
  {
    type = "restricted-" + space->type;

    for (auto vb : { VOL,BND, BBND, BBBND })
    {
      if (space->GetEvaluator(vb))
        evaluator[vb] = make_shared<RestrictedDifferentialOperator>(space->GetEvaluator(vb));
      if (space->GetFluxEvaluator(vb))
        flux_evaluator[vb] = make_shared<RestrictedDifferentialOperator>(space->GetFluxEvaluator(vb));
      integrator[vb] = FESpace::GetIntegrator(vb);
    }

  }

  void RestrictedFESpace::Update()
  { 
    //space->Update(lh); // removed as it may override changed doftypes
    int ne = space->GetMeshAccess()->GetNE();
    int ndof = space->GetNDof();
    active_dofs = make_shared<BitArray> (ndof);
    if (active_els)
    {
      active_dofs->Clear();

	    ParallelForRange (ne, [&](IntRange r)
      {      
        Array<DofId> dnums;
        for (auto elnr : r)
        {
          ElementId elid(VOL,elnr);
          if (active_els->Test(elnr))
          {
            space->GetDofNrs(elid,dnums);
            for (auto dof : dnums)
              active_dofs->SetBitAtomic(dof);
          }
        }
      });
    }
    else
      active_dofs->Set();
    CompressedFESpace::Update();
  }

  FiniteElement & RestrictedFESpace::GetFE (ElementId ei, Allocator & lh) const
  {
    if ((ei.VB() != VOL) || (!active_els || active_els->Test(ei.Nr())))
      return space->GetFE(ei,lh);
    else
      return SwitchET(ma->GetElType(ei), [&lh] (auto type) -> FiniteElement&
                      { return * new (lh) DummyFE<type.ElementType()>(); });
  }

  void RestrictedFESpace::GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    if ((ei.VB() != VOL) || (!active_els || active_els->Test(ei.Nr())))
    {
      space->GetDofNrs(ei,dnums);
      WrapDofs(dnums);
    }
    else
    {
      dnums.SetSize(0);
    }
  }

  void RestrictedFESpace::GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const
  {
    if ((ei.VB() != VOL) || (!active_els || active_els->Test(ei.Nr())))
    {
      space->GetElementDofsOfType(ei,dnums,ctype);
      WrapDofs(dnums);
    }
    else  
      dnums.SetSize(0);
  }


  ProxyNode RestrictedFESpace::MakeProxyFunction (bool testfunction,
                                                  const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const
  {
    auto proxybase = GetBaseSpace()->GetProxyFunction(testfunction);
    auto compound = dynamic_pointer_cast<CompoundDifferentialOperator>(proxybase);
    shared_ptr<DifferentialOperator> ev1 = proxybase->Evaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->Evaluator()) : nullptr;
    shared_ptr<DifferentialOperator> ev2 = proxybase->DerivEvaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->DerivEvaluator()) : nullptr;
    shared_ptr<DifferentialOperator> ev3 = proxybase->TraceEvaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->TraceEvaluator()) : nullptr;
    shared_ptr<DifferentialOperator> ev4 = proxybase->TraceDerivEvaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->TraceDerivEvaluator()) : nullptr;
    shared_ptr<DifferentialOperator> ev5 = proxybase->TTraceEvaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->TTraceEvaluator()) : nullptr;
    shared_ptr<DifferentialOperator> ev6 = proxybase->TTraceDerivEvaluator() ? make_shared<RestrictedDifferentialOperator>(proxybase->TTraceDerivEvaluator()) : nullptr;

    shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<RestrictedFESpace*>(this)->shared_from_this());
    auto proxy = make_shared<ProxyFunction>(fes, testfunction, fes->IsComplex(), ev1, ev2, ev3, ev4, ev5, ev6); 
    ProxyNode proxyn(proxy);
    proxyn.SetFESpace(fes);

    auto add_diffops = proxybase->GetAdditionalEvaluators();
    for (int i = 0; i < add_diffops.Size(); i++)
      proxyn->SetAdditionalEvaluator (add_diffops.GetName(i), make_shared<RestrictedDifferentialOperator>(add_diffops[i]));
    proxyn = addblock(proxyn);

    return proxyn;
  }


}

#include <fem.hpp>
#include "diffop.hpp"
#include "diffop_impl.hpp"

namespace ngfem
{
  void RestrictedDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              BareSliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    if (fel.GetNDof() == 0)
      ; // mat = 0;
    else
      diffop->CalcMatrix (fel, mip, mat, lh);
  }
  
  void RestrictedDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              BareSliceMatrix<Complex,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    if (fel.GetNDof() == 0)
      ; // mat = 0;
    else
      diffop->CalcMatrix (fel, mip, mat, lh);
  }
  
  void RestrictedDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    if (fel.GetNDof() == 0)
      mat.AddSize(fel.GetNDof(),fel.GetNDof()) = 0.0;
    else
      diffop->CalcMatrix(fel, mir, mat);
  }
  
  void RestrictedDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              BareSliceMatrix<SIMD<Complex>> mat) const
  {
    if (fel.GetNDof() == 0)
      mat.AddSize(fel.GetNDof(),fel.GetNDof()) = 0.0;
    else
      diffop->CalcMatrix(fel, mir, mat);
  }

  void RestrictedDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->Apply(fel, mip, x, flux, lh);
  }
  
  void RestrictedDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->Apply(fel, mip, x, flux, lh);
  }

  void RestrictedDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<SIMD<double>> flux) const
  // LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux.AddSize(fel.GetNDof(),diffop->Dim()) = 0.0;
    else
      diffop->Apply(fel, mir, x, flux);
  }
 
  void RestrictedDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<Complex> x, 
         BareSliceMatrix<SIMD<Complex>> flux) const
  // LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux.AddSize(fel.GetNDof(),diffop->Dim()) = 0.0;
    else
      diffop->Apply(fel, mir, x, flux);
  }
 
  void RestrictedDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> x, 
              LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->ApplyTrans(fel, mip, flux, x, lh);
  }

  void RestrictedDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              BareSliceVector<Complex> x, 
              LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->ApplyTrans(fel, mip, flux, x, lh);
  }

  void RestrictedDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      FlatMatrix<double> flux,
	      BareSliceVector<double> x, 
	      LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->ApplyTrans(fel, mir, flux, x, lh);
  }

  void RestrictedDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      FlatMatrix<Complex> flux,
	      BareSliceVector<Complex> x, 
	      LocalHeap & lh) const
  {
    if (fel.GetNDof() == 0)
      flux = 0.0;
    else
      diffop->ApplyTrans(fel, mir, flux, x, lh);
  }

  void RestrictedDifferentialOperator ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    if (fel.GetNDof() != 0)
      diffop->AddTrans(fel, mir, flux, x);
  }

  void RestrictedDifferentialOperator ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<Complex>> flux,
            BareSliceVector<Complex> x) const
  {
    if (fel.GetNDof() != 0)
      diffop->AddTrans(fel, mir, flux, x);
  }


  shared_ptr<CoefficientFunction> RestrictedDifferentialOperator ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir,
             bool Eulerian) const 
  {
    return diffop->DiffShape(proxy, dir, Eulerian);
  }


}
