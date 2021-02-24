#include <comp.hpp>
#include "restrictedfespace.hpp"

namespace ngcomp
{

  RestrictedFESpace::RestrictedFESpace (shared_ptr<FESpace> bfes, shared_ptr<BitArray> _active_els)
      : CompressedFESpace (bfes), active_els(_active_els)
  {
    type = "restricted-" + space->type;
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
    if (!active_els || active_els->Test(ei.Nr()))
      return space->GetFE(ei,lh);
    else
      return SwitchET(ma->GetElType(ei), [&lh] (auto type) -> FiniteElement&
                      { return * new (lh) DummyFE<type.ElementType()>(); });
  }

  void RestrictedFESpace::GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    if (!active_els || active_els->Test(ei.Nr()))
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
    if (!active_els || active_els->Test(ei.Nr()))
      space->GetElementDofsOfType(ei,dnums,ctype);
    else  
      dnums.SetSize(0);
  }



}
