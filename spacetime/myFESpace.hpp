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

  class SpaceTimeFESpace : public FESpace
  {
    int ndof;
    FESpace* Vh;
    ScalarFiniteElement<1>* tfe;
    double time;
    bool override_time = false;

  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    //MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);
    SpaceTimeFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> aVh, shared_ptr<ScalarFiniteElement<1>> atfe, const Flags & flags);

    //virtual FESpace get_V() const { return *Vh; }

    // destructor
    virtual ~SpaceTimeFESpace ();


    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "SpaceTimeFESpace";
    }


    virtual void Update(LocalHeap & lh);
    virtual size_t GetNDof () const { return ndof; }
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;
    FiniteElement* GetTimeFE() { return tfe; }

    // For debugging
    void SetTime(double a) {time = a; override_time = true;}
    void SetOverrideTime(bool a) {override_time = a;}
    // Provide Info for Python
    int order_time()
    {    NodalTimeFE* time_FE = dynamic_cast< NodalTimeFE*>(tfe);
         return time_FE->order_time();
    }
    void TimeFE_nodes(Vector<>& intp_pts)
    {    NodalTimeFE* time_FE = dynamic_cast< NodalTimeFE*>(tfe);
         time_FE->GetIntpPts (intp_pts);
    }
  };

}    

#endif
