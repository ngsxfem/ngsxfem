#pragma once

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngsxfem
#include "../cutint/xintegration.hpp"
// #include "../xfem/xfiniteelement.hpp"

using namespace ngsolve;
using namespace xintegration;
// using namespace cutinfo;

namespace ngcomp
{

  template <int D> class T_XFESpace; // forward declaration

  class CutInformation
  {
    template <int D> friend class T_XFESpace;
  protected:
    shared_ptr<MeshAccess> ma;
    shared_ptr<VVector<double>> cut_ratio_of_element [2] = {nullptr, nullptr};
    shared_ptr<BitArray> elems_of_domain_type [N_COMBINED_DOMAIN_TYPES] = {nullptr, nullptr, nullptr};
    shared_ptr<BitArray> selems_of_domain_type [N_COMBINED_DOMAIN_TYPES] = {nullptr, nullptr, nullptr};
    shared_ptr<BitArray> facets_of_domain_type [N_COMBINED_DOMAIN_TYPES] = {nullptr, nullptr, nullptr};
    shared_ptr<BitArray> cut_neighboring_node [6] = {nullptr, nullptr, nullptr,
                                                     nullptr, nullptr, nullptr};
    shared_ptr<Array<DOMAIN_TYPE>> dom_of_node [6] = {nullptr, nullptr, nullptr,
                                                      nullptr, nullptr, nullptr};
  public:
    CutInformation (shared_ptr<MeshAccess> ama);
    void Update(shared_ptr<CoefficientFunction> lset, int subdivlvl, int time_order, LocalHeap & lh);
    
    shared_ptr<MeshAccess> GetMesh () const { return ma; }

    shared_ptr<BaseVector> GetCutRatios (VorB vb) const
    {
      return cut_ratio_of_element[vb];
    }

    INLINE DOMAIN_TYPE DomainTypeOfElement(ElementId elid)
    {
      int elnr = elid.Nr();
      VorB vb = elid.VB();
      shared_ptr<BitArray> * ba = nullptr;
      if (vb == VOL)
        ba = elems_of_domain_type;
      else
        ba = selems_of_domain_type;
      if (ba[CDOM_IF]->Test(elnr))
        return IF;
      else
      {
        if (ba[CDOM_NEG]->Test(elnr))
          return NEG;
        else
          return POS;
      }
    }

    // template <NODE_TYPE NT>
    // DOMAIN_TYPE GetDomainOfNode (int nr)
    // {
    //   if (GetCutRatioOfNode<NT>(nr) > 0.0 && GetCutRatioOfNode<NT>(nr) < 1.0)
    //     return IF;
    //   else
    //     return GetCutRatioOfNode<NT>(nr) == 0.0 ? NEG : POS;
    // }

    shared_ptr<BitArray> GetElementsOfDomainType(COMBINED_DOMAIN_TYPE dt, VorB vb) const
    {
      if (vb == VOL)
        return elems_of_domain_type[dt];
      else
        return selems_of_domain_type[dt];
    }

    shared_ptr<BitArray> GetElementsOfDomainType(DOMAIN_TYPE dt, VorB vb) const
    {
      return GetElementsOfDomainType(TO_CDT(dt),vb);
    }
    
    shared_ptr<BitArray> GetFacetsOfDomainType(COMBINED_DOMAIN_TYPE dt) const { return facets_of_domain_type[dt]; }
    shared_ptr<BitArray> GetFacetsOfDomainType(DOMAIN_TYPE dt) const { return facets_of_domain_type[TO_CDT(dt)]; }

    shared_ptr<BitArray> GetElementsWithThresholdContribution(DOMAIN_TYPE dt, double threshold, VorB vb);
  };



  class MultiLevelsetCutInformation
  {
  protected:
    shared_ptr<MeshAccess> ma;
    Array<shared_ptr<GridFunction>> lsets;
    vector<tuple<shared_ptr<BitArray>, Array<Array<DOMAIN_TYPE>>, VorB >> collect_elements_with_contribution;
    vector<tuple<shared_ptr<BitArray>, Array<Array<DOMAIN_TYPE>>, VorB >> collect_elements_of_domain_type;
  public:
    MultiLevelsetCutInformation (shared_ptr<MeshAccess> ama, 
                                 const Array<shared_ptr<GridFunction>> & lsets_in);
    
    void Update(const Array<shared_ptr<GridFunction>> & lsets_in, LocalHeap & lh);

    shared_ptr<MeshAccess> GetMesh () const { return ma; }
    
    int GetLen() const {return lsets.Size();}

    void UpdateElementsOfDomainType(const shared_ptr<BitArray> & elems_of_domain_type, 
                                    const Array<Array<DOMAIN_TYPE>> & cdt, 
                                    VorB vb, 
                                    LocalHeap & lh) const;

    shared_ptr<BitArray> GetElementsOfDomainType(const Array<Array<DOMAIN_TYPE>> & cdt,
                                                 VorB vb, 
                                                 LocalHeap & lh);

    shared_ptr<BitArray> GetElementsOfDomainType(const Array<DOMAIN_TYPE> & cdt_,
                                                 VorB vb, 
                                                 LocalHeap & lh)
    {
      Array<Array<DOMAIN_TYPE>> cdt(1);
      cdt[0] = cdt_;
      return GetElementsOfDomainType(cdt,vb,lh);
    }

    void UpdateElementsWithContribution(const shared_ptr<BitArray> & elems_of_domain_type, 
                                        const Array<Array<DOMAIN_TYPE>> & cdt, 
                                        VorB vb, 
                                        LocalHeap & lh) const;

    shared_ptr<BitArray> GetElementsWithContribution(const Array<Array<DOMAIN_TYPE>> & cdt,
                                                     VorB vb, 
                                                     LocalHeap & lh);

    shared_ptr<BitArray> GetElementsWithContribution(const Array<DOMAIN_TYPE> & cdt_,
                                                     VorB vb, 
                                                     LocalHeap & lh)
    {
      Array<Array<DOMAIN_TYPE>> cdt(1);
      cdt[0] = cdt_;
      return GetElementsWithContribution(cdt,vb,lh);
    }

    bool CombinedDomainTypesEqual(const Array<Array<DOMAIN_TYPE>> & cdta, 
                                  const Array<Array<DOMAIN_TYPE>> & cdtb) const;
  };
  
  shared_ptr<BitArray> GetFacetsWithNeighborTypes(shared_ptr<MeshAccess> ma,
                                                  shared_ptr<BitArray> a,
                                                  shared_ptr<BitArray> b,
                                                  bool bound_val_a,
                                                  bool bound_val_b,
                                                  bool ask_and,
                                                  LocalHeap & lh);

  shared_ptr<BitArray> GetElementsWithNeighborFacets(shared_ptr<MeshAccess> ma,
                                                     shared_ptr<BitArray> a,
                                                     LocalHeap & lh);

  shared_ptr<BitArray> GetElementsWithSharedVertex(shared_ptr<MeshAccess> ma, shared_ptr<BitArray> a, LocalHeap & lh);

  shared_ptr<BitArray> GetDofsOfElements(shared_ptr<FESpace> fes,
                                         shared_ptr<BitArray> a,
                                         LocalHeap & lh);

  shared_ptr<BitArray> GetDofsOfFacets(shared_ptr<FESpace> fes,
                                       shared_ptr<BitArray> a,
                                       LocalHeap & lh);


}
