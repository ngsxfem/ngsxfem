/// from ngxfem
#include "../xfem/cutinfo.hpp"
#include "../cutint/xintegration.hpp"
using namespace ngsolve;
using namespace xintegration;
using namespace ngfem;

namespace ngcomp
{

  CutInformation::CutInformation (shared_ptr<MeshAccess> ama)
    : ma(ama)
  {

    for (auto cdt : all_cdts)
    {
      elems_of_domain_type[cdt] = make_shared<BitArray>(ma->GetNE(VOL));
      selems_of_domain_type[cdt] = make_shared<BitArray>(ma->GetNE(BND));
      facets_of_domain_type[cdt] = make_shared<BitArray>(ma->GetNFacets());
    }
    facets_of_domain_type[NEG]->Set();
    facets_of_domain_type[POS]->Clear();
    facets_of_domain_type[IF]->Clear();

    for (NODE_TYPE nt : {NT_VERTEX,NT_EDGE,NT_FACE,NT_CELL})
    {
      cut_neighboring_node[nt] = make_shared<BitArray>(ma->GetNNodes (nt));
      cut_neighboring_node[nt]->Clear();
      dom_of_node[nt] = make_shared<Array<DOMAIN_TYPE>>(ma->GetNNodes (nt));
      (*dom_of_node[nt]) = NEG;
    }

    if (ma->GetDimension() == 3)
    {
      cut_neighboring_node[NT_ELEMENT] = cut_neighboring_node[NT_CELL];
      cut_neighboring_node[NT_FACET] = cut_neighboring_node[NT_FACE];
      dom_of_node[NT_ELEMENT] = dom_of_node[NT_CELL];
      dom_of_node[NT_FACET] = dom_of_node[NT_FACE];
    }
    else
    {
      cut_neighboring_node[NT_ELEMENT] = cut_neighboring_node[NT_FACE];
      cut_neighboring_node[NT_FACET] = cut_neighboring_node[NT_EDGE];
      dom_of_node[NT_ELEMENT] = dom_of_node[NT_FACE];
      dom_of_node[NT_FACET] = dom_of_node[NT_EDGE];
    }
    for (VorB vb : {VOL,BND})
    {
      int ne = ma->GetNE(vb);
      cut_ratio_of_element[vb] = make_shared<VVector<double>>(ne);
    }
  }

  void CutInformation::Update(shared_ptr<CoefficientFunction> cf_lset,int time_order, LocalHeap & lh)
  {
    shared_ptr<GridFunction> gf_lset;
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(cf_lset,subdivlvl);

    for (auto cdt : all_cdts)
    {
      elems_of_domain_type[cdt]->Clear();
      selems_of_domain_type[cdt]->Clear();
    }
    elems_of_domain_type[CDOM_ANY]->Set();
    selems_of_domain_type[CDOM_ANY]->Set();

    for (VorB vb : {VOL,BND})
    {
      int ne = ma->GetNE(vb);
      IterateRange
        (ne, lh,
        [&] (int elnr, LocalHeap & lh)
      {
        ElementId ei = ElementId(vb,elnr);
        Ngs_Element ngel = ma->GetElement(ei);
        ELEMENT_TYPE eltype = ngel.GetType();
        ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

        double part_vol [] = {0.0, 0.0};
        for (DOMAIN_TYPE np : {POS, NEG})
        {
          const IntegrationRule * ir_np = CreateCutIntegrationRule(cf_lset, gf_lset, eltrans, np, 0,time_order, lh, subdivlvl);
          // If(time_order > -1 && vb == BND) should have part_vol[NEG] == 0, which will lead to
          // the BND element being marked as POS.
          if (ir_np)
            for (auto ip : *ir_np)
              part_vol[np] += ip.Weight();
        }
        (*cut_ratio_of_element[vb])(elnr) = part_vol[NEG]/(part_vol[NEG]+part_vol[POS]);                 
        if (vb == VOL)
        {
          if (part_vol[NEG] > 0.0)
            if (part_vol[POS] > 0.0)
              (*elems_of_domain_type[CDOM_IF]).Set(elnr);
            else
              (*elems_of_domain_type[CDOM_NEG]).Set(elnr);
          else
            (*elems_of_domain_type[CDOM_POS]).Set(elnr);
        }
        else
        {
          if (part_vol[NEG] > 0.0)
            if (part_vol[POS] > 0.0)
              (*selems_of_domain_type[CDOM_IF]).Set(elnr);
            else
              (*selems_of_domain_type[CDOM_NEG]).Set(elnr);
          else
            (*selems_of_domain_type[CDOM_POS]).Set(elnr);
        }

      });
      *elems_of_domain_type[CDOM_UNCUT] = *elems_of_domain_type[CDOM_NEG] | *elems_of_domain_type[CDOM_POS];
      *elems_of_domain_type[CDOM_HASNEG] = *elems_of_domain_type[CDOM_NEG] | *elems_of_domain_type[CDOM_IF];
      *elems_of_domain_type[CDOM_HASPOS] = *elems_of_domain_type[CDOM_POS] | *elems_of_domain_type[CDOM_IF];
      *selems_of_domain_type[CDOM_UNCUT] = *selems_of_domain_type[CDOM_NEG] | *selems_of_domain_type[CDOM_POS];
      *selems_of_domain_type[CDOM_HASNEG] = *selems_of_domain_type[CDOM_NEG] | *selems_of_domain_type[CDOM_IF];
      *selems_of_domain_type[CDOM_HASPOS] = *selems_of_domain_type[CDOM_POS] | *selems_of_domain_type[CDOM_IF];
    }

    int ne = ma -> GetNE();
    IterateRange
      (ne, lh,
      [&] (int elnr, LocalHeap & lh)
    {
      if ((*elems_of_domain_type[CDOM_IF]).Test(elnr))
      {
        ElementId elid(VOL,elnr);

        Array<int> nodenums(0,lh);

        nodenums = ma->GetElVertices(elid);
        for (int node : nodenums)
          cut_neighboring_node[NT_VERTEX]->Set(node);

        nodenums = ma->GetElEdges(elid);
        for (int node : nodenums)
          cut_neighboring_node[NT_EDGE]->Set(node);

        if (ma->GetDimension() == 3)
        {
          nodenums = ma->GetElFaces(elid.Nr());
          for (int node : nodenums)
            cut_neighboring_node[NT_FACE]->Set(node);
        }
        cut_neighboring_node[NT_ELEMENT]->Set(elnr);
      }
    });

    for (NODE_TYPE nt : {NT_VERTEX,NT_EDGE,NT_FACE,NT_CELL})
      *(dom_of_node[nt]) = IF;

    IterateRange
      (ne, lh,
      [&] (int elnr, LocalHeap & lh)
    {
      if (elems_of_domain_type[CDOM_UNCUT]->Test(elnr))
      {
        ElementId elid(VOL,elnr);
        Array<int> nodenums(0,lh);

        DOMAIN_TYPE dt = elems_of_domain_type[CDOM_NEG]->Test(elnr) ? NEG : POS;

        nodenums = ma->GetElVertices(elid);
        for (int node : nodenums)
          (*dom_of_node[NT_VERTEX])[node] = dt;

        nodenums = ma->GetElEdges(elid);
        for (int node : nodenums)
          (*dom_of_node[NT_EDGE])[node] = dt;

        if (ma->GetDimension() == 3)
        {
          nodenums = ma->GetElFaces(elid.Nr());
          for (int node : nodenums)
            (*dom_of_node[NT_FACE])[node] = dt;
          (*dom_of_node[NT_CELL])[elnr] = (*cut_ratio_of_element[VOL])(elnr) > 0.5 ? NEG : POS;
        }
        else
        {
          (*dom_of_node[NT_FACE])[elnr] = (*cut_ratio_of_element[VOL])(elnr) > 0.5 ? NEG : POS;
        }
      }

    });

  }


  shared_ptr<BitArray> GetFacetsWithNeighborTypes(shared_ptr<MeshAccess> ma,
                                                  shared_ptr<BitArray> a,
                                                  shared_ptr<BitArray> b,
                                                  bool bound_val_a,
                                                  bool bound_val_b,
                                                  bool ask_and,
                                                  LocalHeap & lh)
  {
    int nf = ma->GetNFacets();
    shared_ptr<BitArray> ret = make_shared<BitArray> (nf);
    ret->Clear();

    BitArray fine_facet(nf);
    fine_facet.Clear();
    IterateRange
      (ma->GetNE(VOL), lh,
      [&] (int elnr, LocalHeap & lh)
    {
      Array<int> fanums(0,lh);
      fanums = ma->GetElFacets (ElementId(VOL,elnr));
      for (int j=0; j<fanums.Size(); j++)
        fine_facet.Set(fanums[j]);
    });

    IterateRange
      (nf, lh,
      [&] (int facnr, LocalHeap & lh)
    {
      if (fine_facet.Test(facnr))
      {
        Array<int> elnums(0,lh);
        ma->GetFacetElements (facnr, elnums);

        if(elnums.Size() < 2)
        {
          int facet2 = ma->GetPeriodicFacet(facnr);
          if(facet2 > facnr)
          {
            Array<int> elnums_per(1,lh);
            ma->GetFacetElements (facet2, elnums_per);
            elnums.Append(elnums_per[0]);
          }
          else
            return;
        }
        
        bool a_left = a->Test(elnums[0]);
        bool a_right = elnums.Size() > 1 ? a->Test(elnums[1]) : bound_val_a;
        bool b_left = b->Test(elnums[0]);
        bool b_right = elnums.Size() > 1 ? b->Test(elnums[1]) : bound_val_b;

        if (ask_and)
        {
          if ((a_left && b_right) || (a_right && b_left))
            ret->Set(facnr);
        }
        else
        {
          if ((a_left || b_right) || (a_right || b_left))
            ret->Set(facnr);
        }
      }
    });
    return ret;
  }

  shared_ptr<BitArray> GetElementsWithNeighborFacets(shared_ptr<MeshAccess> ma,
                                                     shared_ptr<BitArray> a,
                                                     LocalHeap & lh)
  {
    int nf = ma->GetNFacets();
    int ne = ma->GetNE();
    shared_ptr<BitArray> ret = make_shared<BitArray> (ne);
    ret->Clear();

    IterateRange
      (nf, lh,
      [&] (int facnr, LocalHeap & lh)
    {
      if (a->Test(facnr))
      {
        Array<int> elnums(0,lh);
        ma->GetFacetElements (facnr, elnums);
        for (auto elnr : elnums)
          ret->Set(elnr);
      }
    });
    return ret;
  }

  shared_ptr<BitArray> GetDofsOfElements(shared_ptr<FESpace> fes,
                                         shared_ptr<BitArray> a,
                                         LocalHeap & lh)
  {
    int ne = fes->GetMeshAccess()->GetNE();
    int ndof = fes->GetNDof();
    shared_ptr<BitArray> ret = make_shared<BitArray> (ndof);
    ret->Clear();

    IterateRange
      (ne, lh,
      [&] (int elnr, LocalHeap & lh)
    {
      ElementId elid(VOL,elnr);
      if (a->Test(elnr))
      {
        Array<int> dnums(0,lh);
        fes->GetDofNrs(elid,dnums);
        for (auto dof : dnums)
          ret->Set(dof);
      }
    });
    return ret;
  }


}
