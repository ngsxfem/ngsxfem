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

    ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override
    {
      return GetBaseSpace()->MakeProxyFunction (testfunction, addblock);
    }

  };

}
