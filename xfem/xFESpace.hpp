#ifndef FILE_XFESPACE_HPP
#define FILE_XFESPACE_HPP

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngxfem
#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimefespace.hpp"

using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{

  // Base class for extended finite elements with data for 
  // mappings between degrees of freedoms and level set information
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
      
    shared_ptr<Table<int>> el2dofs = NULL;
    shared_ptr<Table<int>> sel2dofs = NULL;

    Array<DOMAIN_TYPE> domofdof;
    Array<DOMAIN_TYPE> domofel;
    Array<DOMAIN_TYPE> domofsel;
    Array<DOMAIN_TYPE> domofface;
    Array<DOMAIN_TYPE> domofedge;
    Array<DOMAIN_TYPE> domofvertex;
    Array<DOMAIN_TYPE> domofinner;

    Array<int> basedof2xdof;
    Array<int> xdof2basedof;

    // Table<int> sel2dofs;
    shared_ptr<FESpace> basefes = NULL;
    BitArray activeelem;
    BitArray activeselem;
      
    shared_ptr<CoefficientFunction> coef_lset = NULL;

    double vmax = 1e99;

    bool trace = false;
    bool empty = false;
  public:
    XFESpace (shared_ptr<MeshAccess> ama, const Flags & flags): FESpace(ama, flags){;}
    
    void CleanUp();
    virtual ~XFESpace(){CleanUp();};

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "XFESpace";
    }

    // function which can create fe-space (needed by pde-parser)
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

    virtual int GetNDof () const { return ndof; }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    void GetDomainNrs (int elnr, Array<DOMAIN_TYPE> & domnums) const;
    void GetSurfaceDomainNrs (int selnr, Array<DOMAIN_TYPE> & domnums) const;

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      if (activeelem.Size() == 0) return;
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
      if (activeelem.Size() == 0) return;
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
      if (activeelem.Size() == 0) return;
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
      if (activeelem.Size() == 0) return;
      Array<int> ldnums;
      basefes->GetInnerDofNrs(vnr,ldnums);
      for (int i = 0; i < ldnums.Size(); ++i)
      {
        int dof = basedof2xdof[ldnums[i]];
        if (dof!=-1)
          dnums.Append(dof);
      }
    }

    virtual void UpdateCouplingDofArray();

    virtual bool DefinedOn (ElementId id) const
    {
      if (activeelem.Size() == 0) return false;

      if (id.IsBoundary())
        return activeselem.Test(int(id));
      else
        return activeelem.Test(int(id));
    }

    void SetBaseFESpace(shared_ptr<FESpace> basefes_){basefes = basefes_;};
    void SetLevelSet(shared_ptr<GridFunction> lset_){ 
      coef_lset = make_shared<GridFunctionCoefficientFunction>(lset_);
    };
    void SetLevelSet(shared_ptr<CoefficientFunction> _coef_lset){ coef_lset = _coef_lset;};
    void SetTimeInterval( const TimeInterval & a_ti){ ti = a_ti;};
    
    bool IsElementCut(int elnr) const { return activeelem.Test(elnr); }
    bool IsNeighborElementCut(int elnr) const 
    { 
      Array<int> faces(0);
      ma->GetElFacets(elnr, faces);
      for (int fa : faces)
      {
        Array<int> els(0);
        ma->GetFacetElements(fa,els);
        for (int el : els)
          if (el == elnr)
            continue;
          else
            if (IsElementCut(el)) return true;
      }
      return false;
    }
    // void SetGlobalCutInfo(AdLinCutTriang* gci_){ gci = gci_;};

    DOMAIN_TYPE GetDomOfDof(int n) const { return domofdof[n];}
    int GetBaseDofOfXDof(int n) const { return xdof2basedof[n];}
    int GetXDofOfBaseDof(int n) const { return basedof2xdof[n];}
    
    static void XToNegPos(shared_ptr<GridFunction> gf, shared_ptr<GridFunction> gf_neg_pos);
  };


  // XFESpace with implementations for current (space-time) dimensions D, SD
  template <int D, int SD>
  class T_XFESpace : public XFESpace
  {
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    T_XFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);
    
    // destructor
    virtual ~T_XFESpace ();

    virtual void Update(LocalHeap & lh);

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;

  };



  class LevelsetContainerFESpace : public FESpace
  {
    shared_ptr<CoefficientFunction> coef_lset = NULL;
    double told;
    double tnew;
  public:
    LevelsetContainerFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);
    virtual ~LevelsetContainerFESpace () { ; }
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<LevelsetContainerFESpace>(ma,flags);
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

    void SetLevelSet(shared_ptr<CoefficientFunction> _coef_lset)
    { coef_lset = _coef_lset;}
    void SetTime(double ta, double tb) { told=ta; tnew=tb; }
    virtual string GetClassName () const { return "LevelsetContainerFESpace"; }
  };  


  class NumProcInformXFESpace : public NumProc
  {
    const CoefficientFunction * coef = NULL;
  public:
    NumProcInformXFESpace (shared_ptr<PDE> apde, const Flags & flags);
    ~NumProcInformXFESpace();
    virtual string GetClassName () const;
    virtual void Do (LocalHeap & lh);
  };

  class NumProcXToNegPos : public NumProc
  {
    shared_ptr<GridFunction> gfxstd = nullptr;
    shared_ptr<GridFunction> gfnegpos = nullptr;
  public:
    NumProcXToNegPos (shared_ptr<PDE> apde, const Flags & flags);
    ~NumProcXToNegPos(){;};
    virtual string GetClassName () const{ return "NumProcXToNegPos";};
    virtual void Do (LocalHeap & lh);
  };


  class XStdFESpace : public CompoundFESpace
  {
    bool spacetime = false;
  public:
    XStdFESpace (shared_ptr<MeshAccess> ama, 
                const Array<shared_ptr<FESpace> > & aspaces,
                const Flags & flags);
    virtual ~XStdFESpace () { ; }
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      bool spacetime = flags.GetDefineFlag("spacetime");
      Array<shared_ptr<FESpace> > spaces(2);
      if (spacetime)
      {
        shared_ptr<FESpaceClasses::FESpaceInfo> info;
        string fet_space = flags.GetStringFlag("type_std","h1ho");
        Flags fespaceflags(flags);
        fespaceflags.SetFlag("type_space",fet_space);
        spaces[0] = make_shared<SpaceTimeFESpace> (ma, fespaceflags);    
      }
      else
      {
        shared_ptr<FESpaceClasses::FESpaceInfo> info;
        string fet_space = flags.GetStringFlag("type_std","h1ho");
        info = GetFESpaceClasses().GetFESpace(fet_space);
        if (!info) throw Exception("XStdFESpace ::  XStdFESpace : fespace not given ");
        Flags fespaceflags(flags);
        spaces[0] = info->creator(ma, fespaceflags);
      }

      spaces[1] = XFESpace::Create(ma,flags);
      shared_ptr<XStdFESpace> fes = make_shared<XStdFESpace> (ma, spaces, flags);
      return fes;
    }

    Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    Array<int> * CreateDirectSolverClusters (const Flags & flags) const;
    bool IsSpaceTime() const { return spacetime;}
    virtual string GetClassName () const { return "XStdFESpace"; }
    void SetLevelSet(shared_ptr<GridFunction> lset_){ 
      dynamic_pointer_cast<XFESpace>(spaces[1])->SetLevelSet(lset_);
    };
    void SetLevelSet(shared_ptr<CoefficientFunction> lset_){ 
      dynamic_pointer_cast<XFESpace>(spaces[1])->SetLevelSet(lset_);
    };
    shared_ptr<XFESpace> GetXFESpace() const { return dynamic_pointer_cast<XFESpace>(spaces[1]);}
    shared_ptr<FESpace> GetStdFESpace() const { return spaces[0];}
    
  };
  
}    

#endif
