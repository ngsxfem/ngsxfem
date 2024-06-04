#ifndef FILE_XFESPACE_HPP
#define FILE_XFESPACE_HPP

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngsxfem
#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../xfem/cutinfo.hpp"

using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{

  // Base class for extended finite elements with data for
  // mappings between degrees of freedoms and cut information
  class XFESpace : public FESpace
  {
  protected:
    int ndof;

    shared_ptr<Table<int>> el2dofs = NULL;
    shared_ptr<Table<int>> sel2dofs = NULL;

    Array<DOMAIN_TYPE> domofdof;

    Array<int> basedof2xdof;
    Array<int> xdof2basedof;

    shared_ptr<FESpace> basefes = NULL;

    shared_ptr<CoefficientFunction> coef_lset = NULL; // <- only used if cutinfo is private
    shared_ptr<CutInformation> cutinfo = NULL;
    bool private_cutinfo = true;  // <-- am I responsible for the cutinformation (or is it an external one)
    bool trace = false;   // xfespace is a trace fe space (special case for further optimization (CouplingDofTypes...))
  public:
    shared_ptr<FESpace> GetBaseFESpace() const { return basefes;};

    XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> abasefes,
              shared_ptr<CoefficientFunction> lset, const Flags & flags)
      : FESpace(ama, flags), basefes(abasefes)
    {
      type = "xfes(" + abasefes->type + ")";
    }

    XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> abasefes,
              shared_ptr<CutInformation> acutinfo, const Flags & flags)
      : FESpace(ama, flags), basefes(abasefes), cutinfo(acutinfo)
    {
      type = "xfes(" + abasefes->type + ")";
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
    virtual void FinalizeUpdate();

    virtual size_t GetNDof () const { return ndof; }

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
      if (cutinfo->GetElementsOfDomainType(IF,VOL)->Size() == 0) return;
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
      if (cutinfo->GetElementsOfDomainType(IF,VOL)->Size() == 0) return;
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
      if (cutinfo->GetElementsOfDomainType(IF,VOL)->Size() == 0) return;
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
      if (cutinfo->GetElementsOfDomainType(IF,VOL)->Size() == 0) return;
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
      if (cutinfo->GetElementsOfDomainType(IF,VOL)->Size() == 0)
        return false;
      return cutinfo->GetElementsOfDomainType(IF,id.VB())->Test(id.Nr());
    }

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

    virtual void Update();

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;

    SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const;

  };

}

#endif
