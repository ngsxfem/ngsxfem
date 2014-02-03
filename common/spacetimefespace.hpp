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

namespace ngcomp
{

  class SpaceTimeFESpace : public FESpace 
  {
  protected:  
    /// space fespace
    FESpace * fes_space;
    /// space fespace
    FiniteElement * fel_time;
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
  public:
    ///
    SpaceTimeFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~SpaceTimeFESpace ();
    ///
    virtual string GetClassName () const
    {
      return "SpaceTimeFESpace";
    }
  
    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual void UpdateCouplingDofArray();    
    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int level) const;
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
