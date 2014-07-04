/**
   Finite Element Space for Tensor product Space Time Finite Elements
*/
/* spacetimefespace.*pp
 * An fespace for functions living on D+1-prisms
 *   - the space has nxm dofs with n the dofs of the space fespace and n the dofs of the time fe
 *   - the ordering of dofs is space-wise, i.e. the first m dofs are correspondig to one time fe ansatz
 *   - no low-order space
 *   - no prolongation so far (will probably never be needed)
 */

#ifndef SPACETIME_FESPACE_HPP
#define SPACETIME_FESPACE_HPP

#include <comp.hpp>
#include <fem.hpp> 

namespace ngcomp
{

  class SpaceTimeFESpace : public FESpace 
  {
  protected:  
    /// space fespace
    FESpace * fes_space;
    /// space fespace
    DGFiniteElement<1> * fel_time;
    ///
    int spacedim;  
    /// 
    int order_time;
    /// 
    int order_space;
    ///
    int nel;
    ///
    int nfa;
    ///
    int ndof;
    ///
    int ndof_space;
    ///
    int ndof_time;
    ///
    bool gaussradau;
    ///
    Array<int> ndlevel;
  public:
    ///
    SpaceTimeFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags=false);
    ///
    virtual ~SpaceTimeFESpace ();
    ///
    virtual string GetClassName () const
    {
      return "SpaceTimeFESpace";
    }
  
    ///
    virtual int OrderTime() const { return order_time; }
    ///
    virtual int OrderSpace() const { return order_space; }
    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual void UpdateCouplingDofArray();    
    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int level) const;
    ///
    virtual const FiniteElement & GetSpaceFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetTimeFE (LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const; 
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
    ///
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & precflags) const;

    virtual void GetVertexDofNrs ( int nr, Array<int> & dnums ) const;

    virtual void GetEdgeDofNrs ( int nr, Array<int> & dnums ) const;

    virtual void GetFaceDofNrs (int nr, Array<int> & dnums) const;
  
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;
  
  };

}

#endif
