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
    int order_time = 1;

    int ref_lvl_space = 0;
    int ref_lvl_time = 0;
      
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
    Array<int> xdof2basedof;

    // Table<int> sel2dofs;
    const FESpace * basefes = NULL;
    BitArray activeelem;
    BitArray activeselem;
      
    const GridFunction * gf_lset = NULL;
    const CoefficientFunction * coef_lset = NULL;
    EvalFunction * eval_lset = NULL;

    double vmax = 1e99;

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
    
    void CleanUp();

    virtual void Update(LocalHeap & lh);
    virtual void UpdateCouplingDofArray();

    virtual int GetNDof () const { return ndof; }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    void GetDomainNrs (int elnr, Array<DOMAIN_TYPE> & domnums) const;
    void GetSurfaceDomainNrs (int selnr, Array<DOMAIN_TYPE> & domnums) const;

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      Array<int> ldnums;
      basefes->GetVertexDofNrs(vnr,ldnums);
      for (int i = 0; i < ldnums.Size(); ++i)
      {
        int dof = basedof2xdof[ldnums[i]];
        if (dof!=-1)
          dnums.Append(dof);
      }
    }

    virtual void GetEdgeDofNrs (int vnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      Array<int> ldnums;
      basefes->GetEdgeDofNrs(vnr,ldnums);
      for (int i = 0; i < ldnums.Size(); ++i)
      {
        int dof = basedof2xdof[ldnums[i]];
        if (dof!=-1)
          dnums.Append(dof);
      }
    }

    virtual void GetFaceDofNrs (int vnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      Array<int> ldnums;
      basefes->GetFaceDofNrs(vnr,ldnums);
      for (int i = 0; i < ldnums.Size(); ++i)
      {
        int dof = basedof2xdof[ldnums[i]];
        if (dof!=-1)
          dnums.Append(dof);
      }
    }

    virtual void GetInnerDofNrs (int vnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      Array<int> ldnums;
      basefes->GetInnerDofNrs(vnr,ldnums);
      for (int i = 0; i < ldnums.Size(); ++i)
      {
        int dof = basedof2xdof[ldnums[i]];
        if (dof!=-1)
          dnums.Append(dof);
      }
    }


    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;

    // void SetGlobalCutInfo(AdLinCutTriang* gci_){ gci = gci_;};

    void SetBaseFESpace(const FESpace* basefes_){basefes = basefes_;};
    void SetLevelSet(const GridFunction* lset_){ gf_lset = lset_;};
    void SetBaseFESpace(const FESpace& basefes_){basefes = &basefes_;};
    void SetLevelSet(const GridFunction& lset_){ gf_lset = &lset_;};
    void SetLevelSetCoefficient(const CoefficientFunction* _coef_lset){ coef_lset = _coef_lset;};
    void SetTimeInterval( const TimeInterval & a_ti){ ti = a_ti;};
    
    void XToNegPos(const GridFunction & gf, GridFunction & gf_neg_pos) const;
  };



  class LevelsetContainerFESpace : public FESpace
  {
    const CoefficientFunction * coef_lset = NULL;
    double told;
    double tnew;
  public:
    LevelsetContainerFESpace (const MeshAccess & ama, const Flags & flags);
    virtual ~LevelsetContainerFESpace () { ; }
    static FESpace * Create (const MeshAccess & ma, const Flags & flags)
    {
      return new LevelsetContainerFESpace(ma,flags);
    }
    virtual void Update(LocalHeap & lh) { ; }
    virtual void UpdateCouplingDofArray() { ; }

    virtual int GetNDof () const { return 0; }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const { dnums.SetSize(0); }
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const { dnums.SetSize(0); }

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const 
    { return *new (lh) LevelsetContainerFE(coef_lset,told,tnew); }
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const 
    { return *new (lh) LevelsetContainerFE(coef_lset,told,tnew); }

    void SetLevelSetCoefficient(const CoefficientFunction* _coef_lset)
    { coef_lset = _coef_lset;}
    void SetTime(double ta, double tb) { told=ta; tnew=tb; }

    virtual string GetClassName () const { return "LevelsetContainerFESpace"; }
  };  


  class NumProcInformXFESpace : public NumProc
  {
    const CoefficientFunction * coef = NULL;
  public:
    NumProcInformXFESpace (PDE & apde, const Flags & flags);
    ~NumProcInformXFESpace();
    virtual string GetClassName () const;
    virtual void Do (LocalHeap & lh);
  };


  class XH1FESpace : public CompoundFESpace
  {
  public:
    XH1FESpace (const MeshAccess & ama, 
                const Array<FESpace*> & aspaces,
                const Flags & flags);
    virtual ~XH1FESpace () { ; }
    static FESpace * Create (const MeshAccess & ma, const Flags & flags)
    {
        Array<FESpace*> spaces(2);
        spaces[0] = new H1HighOrderFESpace (ma, flags);    
        if (ma.GetDimension() == 2)
          spaces[1] = new XFESpace<2,2> (ma, flags);        
        else
          spaces[1] = new XFESpace<3,3> (ma, flags);        
        XH1FESpace * fes = new XH1FESpace (ma, spaces, flags);
        return fes;
    }

    Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    Array<int> * CreateDirectSolverClusters (const Flags & flags) const;

    virtual string GetClassName () const { return "XH1FESpace"; }
    
  };
  
}    

#endif
