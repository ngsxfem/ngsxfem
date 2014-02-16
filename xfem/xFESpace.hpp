#ifndef FILE_XFESPACE_HPP
#define FILE_XFESPACE_HPP

/// from ngsolve
#include <comp.hpp>    // provides FESpace, ...
#include <fem.hpp>

/// from ngxfem
// #include "../cuttriang/geom.hpp"
#include "../cuttriang/cuttriang.hpp"
#include "../fem/xfemVisInts.hpp"
#include "../xfiniteelement.hpp"

using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{

  template <int D, int SD>
  class XFESpace : public FESpace
  {
  protected:  
    int ndof;

    Table<int> * el2dofs;
    Table<int> * sel2dofs;
    Array<int> domofdof;
    Array<int> basedof2xdof;
    // Table<int> sel2dofs;
    const FESpace * basefes;
    BitArray activeelem;
    BitArray activeselem;

    const GridFunction * gf_lset;

  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    XFESpace (const MeshAccess & ama, const Flags & flags);
    
    // destructor
    virtual ~XFESpace ();

    // function which can create fe-space (needed by pde-parser)
    static FESpace * Create (const MeshAccess & ma, const Flags & flags)
      {
        // Creator should be outside xFESpace (and translate flags to according template parameters)
        return new XFESpace(ma,flags);
      }

    // a name for our new fe-space
    virtual string GetClassName () const
      {
        return "XFESpace";
      }


    virtual void Update(LocalHeap & lh);
    virtual void UpdateCouplingDofArray();

    virtual int GetNDof () const { return ndof; }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    void GetDomainNrs (int elnr, Array<int> & domnums) const;
    void GetSurfaceDomainNrs (int selnr, Array<int> & domnums) const;

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;

    void SetGlobalCutInfo(AdLinCutTriang* gci_){ gci = gci_;};

    void SetBaseFESpace(const FESpace* basefes_){basefes = basefes_;};
    void SetLevelSet(const GridFunction* lset_){ lset = lset_;};
    void SetBaseFESpace(const FESpace& basefes_){basefes = &basefes_;};
    void SetLevelSet(const GridFunction& lset_){ lset = &lset_;};
  };


  // class XH1FESpace : public CompoundFESpace
  // {
  // public:
  //   XH1FESpace (const MeshAccess & ama, 
  //               const Array<FESpace*> & aspaces,
  //               const Flags & flags);
  //   virtual ~XH1FESpace () { ; }
  //   static FESpace * Create (const MeshAccess & ma, const Flags & flags)
  //     {
  //       Array<FESpace*> spaces(2);
  //       spaces[0] = new H1HighOrderFESpace (ma, flags);    
  //       spaces[1] = new XFESpace (ma, flags);        
  //       XH1FESpace * fes = new XH1FESpace (ma, spaces, flags);
  //       return fes;
  //     }
  //   virtual string GetClassName () const { return "XH1FESpace"; }
    
  // };
  
}    

#endif
