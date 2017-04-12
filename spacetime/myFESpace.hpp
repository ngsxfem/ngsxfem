#ifndef FILE_MYFESPACE_HPP
#define FILE_MYFESPACE_HPP

/*********************************************************************/
/* File:   myFESpace.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

  My own FESpace for linear and quadratic triangular elements An
  fe-space is the connection between the local reference element, and
  the global mesh.

*/


namespace ngcomp
{

  class MyFESpace : public FESpace
  {
    bool use_quads, tp_quads;
    int ndof, nvert,order;
    Array<int> first_edge_dof;
    Array<int> first_cell_dof;

  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // destructor
    virtual ~MyFESpace ();

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "MyFESpace";
    }


    virtual void Update(LocalHeap & lh);
    virtual size_t GetNDof () const { return ndof; }
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;
  };

}    

#endif
