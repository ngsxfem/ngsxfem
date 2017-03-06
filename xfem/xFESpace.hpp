#ifndef FILE_XFESPACE_HPP
#define FILE_XFESPACE_HPP

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngxfem
#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../xfem/cutinfo.hpp"

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
    int nvertdofs; // <- number of x vertex dofs
    int order_space = 1;

    // to be removed:
    int ref_lvl_space = 0;

    shared_ptr<Table<int>> el2dofs = NULL;
    shared_ptr<Table<int>> sel2dofs = NULL;

    Array<DOMAIN_TYPE> domofdof;

    // to be removed:
    Array<DOMAIN_TYPE> domofel;
    // to be removed:
    Array<DOMAIN_TYPE> domofsel;
    // to be removed:
    Array<DOMAIN_TYPE> domofface;
    // to be removed:
    Array<DOMAIN_TYPE> domofedge;
    // to be removed:
    Array<DOMAIN_TYPE> domofvertex;
    // to be removed:
    Array<DOMAIN_TYPE> domofinner;

    Array<int> basedof2xdof;
    Array<int> xdof2basedof;

    // Table<int> sel2dofs;
    shared_ptr<FESpace> basefes = NULL;

    // to be removed:
    BitArray activeelem;
    // to be removed:
    BitArray activeselem;

    shared_ptr<CoefficientFunction> coef_lset = NULL;
    shared_ptr<CutInformation> cutinfo = NULL;
    bool private_cutinfo = true;  // <-- am I responsible for the cutinformation (or is it an external one)

    double vmax = 1e99;

    bool trace = false;   // xfespace is a trace fe space (special case for further optimization (CouplingDofTypes...))
  public:
    void SetBaseFESpace(shared_ptr<FESpace> basefes_){basefes = basefes_;};
    shared_ptr<FESpace> GetBaseFESpace() const { return basefes;};

    void SetLevelSet(shared_ptr<GridFunction> lset_){
      coef_lset = make_shared<GridFunctionCoefficientFunction>(lset_);
    };
    void SetLevelSet(shared_ptr<CoefficientFunction> _coef_lset){ coef_lset = _coef_lset;};
    shared_ptr<CoefficientFunction> GetLevelSet() const { return coef_lset;};

    XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
              shared_ptr<CoefficientFunction> lset, const Flags & flags) : FESpace(ama, flags){
      SetBaseFESpace(basefes);
      SetLevelSet(lset);
    }

    XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
              shared_ptr<CutInformation> acutinfo, const Flags & flags)
      : FESpace(ama, flags), cutinfo(acutinfo)
    {
      SetBaseFESpace(basefes);
      // SetLevelSet(lset);
    }

    shared_ptr<CutInformation> GetCutInfo()
    {
      return cutinfo;
    }

    void CleanUp();
    virtual ~XFESpace(){CleanUp();};

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "XFESpace";
    }

    /// update element coloring
    virtual void FinalizeUpdate(LocalHeap & lh)
    {
      if ( basefes == NULL )
      {
        cout << " no basefes, FinalizeUpdate postponed " << endl;
        return;
      }
      if ( coef_lset == NULL )
      {
        cout << " no lset, FinalizeUpdate postponed " << endl;
        return;
      }
      FESpace::FinalizeUpdate (lh);
    }

    virtual size_t GetNDof () const { return ndof; }
    // virtual int GetNVertexDof () const { return nvertdofs; }

    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const;

    DOMAIN_TYPE GetDomainOfDof (int dof) const { return domofdof[dof]; }
    void GetDomainNrs (ElementId ei, Array<DOMAIN_TYPE> & domnums) const;

    virtual int GetRelOrder() const
    {
      // cout << "virtual GetRelOrder called for FiniteElementSpace, not available ! " << endl;
      return 0;
    }

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
        return activeselem.Test(id.Nr());
      else
        return activeelem.Test(id.Nr());
    }
    DOMAIN_TYPE GetDomainOfElement(int elnr) const {return domofel[elnr];}

    // bool IsElementCut(int elnr) const { return activeelem.Test(elnr); }
    // const BitArray & CutElements() const { return activeelem; }
    // const BitArray & CutSurfaceElements() const { return activeselem; }

    // bool IsNeighborElementCut(int elnr) const
    // {
    //   Array<int> faces(0);
    //   ma->GetElFacets(elnr, faces);
    //   for (int fa : faces)
    //   {
    //     Array<int> els(0);
    //     ma->GetFacetElements(fa,els);
    //     for (int el : els)
    //       if (el == elnr)
    //         continue;
    //       else
    //       if (IsElementCut(el)) return true;
    //   }
    //   return false;
    // }
    // void SetGlobalCutInfo(AdLinCutTriang* gci_){ gci = gci_;};

    DOMAIN_TYPE GetDomOfDof(int n) const { return domofdof[n];}
    int GetBaseDofOfXDof(int n) const { return xdof2basedof[n];}
    int GetXDofOfBaseDof(int n) const { return basedof2xdof[n];}

    static void XToNegPos(shared_ptr<GridFunction> gf, shared_ptr<GridFunction> gf_neg_pos);
  };


  // XFESpace with implementations for dimension D
  template <int D>
  class T_XFESpace : public XFESpace
  {
  public:
    T_XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
                shared_ptr<CoefficientFunction> lset, const Flags & flags);

    T_XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
                shared_ptr<CutInformation> cutinfo, const Flags & flags);

    // destructor
    virtual ~T_XFESpace ();

    virtual void Update(LocalHeap & lh);

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;

    SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const;

  };

}

#endif
