#ifndef FILE_XFESPACE_HPP
#define FILE_XFESPACE_HPP

/// from ngsolve
#include <solve.hpp>    // provides FESpace, ...
#include <comp.hpp>    // provides FESpace, ...
#include <fem.hpp>

/// from ngxfem
#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimefespace.hpp"

using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{

  template <int D, int SD>
  class XFESpace : public FESpace
  {
  protected:  
    int ndof;

    int order_space = 1;
      
    bool spacetime = false;
    TimeInterval ti;
      
    Table<int> * el2dofs = NULL;
    Table<int> * sel2dofs = NULL;

    Array<DOMAIN_TYPE> domofdof;
    Array<DOMAIN_TYPE> domofel;
    Array<DOMAIN_TYPE> domofsel;
    Array<DOMAIN_TYPE> domofface;
    Array<DOMAIN_TYPE> domofedge;
    Array<DOMAIN_TYPE> domofvertex;

    Array<int> basedof2xdof;
    // Table<int> sel2dofs;
    const FESpace * basefes = NULL;
    BitArray activeelem;
    BitArray activeselem;
      
    const GridFunction * gf_lset = NULL;
    const CoefficientFunction * coef_lset = NULL;
    EvalFunction * eval_lset = NULL;

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
      bool spacetime = flags.GetDefineFlag("spacetime");
      int mD = ma.GetDimension();
      int mSD = spacetime ? mD+1 : mD;
      if (mD == 2)
        if (mSD == 2)
          return new XFESpace<2,2>(ma,flags);
        else
          return new XFESpace<2,3>(ma,flags);
      else
        if (mSD == 2)
          return new XFESpace<3,2>(ma,flags);
        else
          return new XFESpace<3,3>(ma,flags);
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

    void GetDomainNrs (int elnr, Array<DOMAIN_TYPE> & domnums) const;
    void GetSurfaceDomainNrs (int selnr, Array<DOMAIN_TYPE> & domnums) const;

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;

    // void SetGlobalCutInfo(AdLinCutTriang* gci_){ gci = gci_;};

    void SetBaseFESpace(const FESpace* basefes_){basefes = basefes_;};
    void SetLevelSet(const GridFunction* lset_){ gf_lset = lset_;};
    void SetBaseFESpace(const FESpace& basefes_){basefes = &basefes_;};
    void SetLevelSet(const GridFunction& lset_){ gf_lset = &lset_;};
    void SetTimeInterval( const TimeInterval & a_ti){ ti = a_ti;};
  };

  
  class NumProcInformXFESpace : public NumProc
  {
  public:
    NumProcInformXFESpace (PDE & apde, const Flags & flags);
    ~NumProcInformXFESpace();
    virtual string GetClassName () const;
    virtual void Do (LocalHeap & lh);
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
