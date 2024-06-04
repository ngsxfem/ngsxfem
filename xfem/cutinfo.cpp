/// from ngsxfem
#include "../xfem/cutinfo.hpp"
#include "../cutint/xintegration.hpp"
#include "../cutint/straightcutrule.hpp"
#include <unordered_map>

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

    if( ma->GetDimension() == 1) {
      cut_neighboring_node[NT_ELEMENT] = cut_neighboring_node[NT_EDGE];
      cut_neighboring_node[NT_FACET] = cut_neighboring_node[NT_VERTEX];
      dom_of_node[NT_ELEMENT] = dom_of_node[NT_EDGE];
      dom_of_node[NT_FACET] = dom_of_node[NT_VERTEX];
    }
    else if (ma ->GetDimension() == 2) {
      cut_neighboring_node[NT_ELEMENT] = cut_neighboring_node[NT_FACE];
      cut_neighboring_node[NT_FACET] = cut_neighboring_node[NT_EDGE];
      dom_of_node[NT_ELEMENT] = dom_of_node[NT_FACE];
      dom_of_node[NT_FACET] = dom_of_node[NT_EDGE];
    }
    else if (ma->GetDimension() == 3) {
      cut_neighboring_node[NT_ELEMENT] = cut_neighboring_node[NT_CELL];
      cut_neighboring_node[NT_FACET] = cut_neighboring_node[NT_FACE];
      dom_of_node[NT_ELEMENT] = dom_of_node[NT_CELL];
      dom_of_node[NT_FACET] = dom_of_node[NT_FACE];
    }
    for (VorB vb : {VOL,BND})
    {
      int ne = ma->GetNE(vb);
      cut_ratio_of_element[vb] = make_shared<VVector<double>>(ne);
    }
  }

  void CutInformation::Update(shared_ptr<CoefficientFunction> cf_lset, int subdivlvl, int time_order, LocalHeap & lh)
  {
    static Timer timer("CutInformation::Update");
    RegionTimer reg (timer);

    if (time_order > 0)
      cout << IM(3) << "Warning: you used CutInfo.Update with time_order > 0.\ntime_order = 0 suffices to detect the cut topology\nand should be prefered.\n";

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
        // Ngs_Element ngel = ma->GetElement(ei);
        // ELEMENT_TYPE eltype = ngel.GetType();
        ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

        double part_vol [] = {0.0, 0.0};
        for (DOMAIN_TYPE np : {POS, NEG})
        {
            const IntegrationRule * ir_np;
            Array<double> wei_arr;
            tie (ir_np, wei_arr) = CreateCutIntegrationRule(cf_lset, gf_lset, eltrans, np, 0,time_order, lh, subdivlvl);
          // If(time_order > -1 && vb == BND) should have part_vol[NEG] == 0, which will lead to
          // the BND element being marked as POS.
          if (ir_np)
            for (auto w : wei_arr) //for (auto ip : *ir_np)
              part_vol[np] += w; // ... += ip.Weight();
        }
        (*cut_ratio_of_element[vb])(elnr) = part_vol[NEG]/(part_vol[NEG]+part_vol[POS]);
        if (vb == VOL)
        {
          if (part_vol[NEG] > 0.0)
            if (part_vol[POS] > 0.0)
              (*elems_of_domain_type[CDOM_IF]).SetBitAtomic(elnr);
            else
              (*elems_of_domain_type[CDOM_NEG]).SetBitAtomic(elnr);
          else
            (*elems_of_domain_type[CDOM_POS]).SetBitAtomic(elnr);
        }
        else
        {
          if (part_vol[NEG] > 0.0)
            if (part_vol[POS] > 0.0)
              (*selems_of_domain_type[CDOM_IF]).SetBitAtomic(elnr);
            else
              (*selems_of_domain_type[CDOM_NEG]).SetBitAtomic(elnr);
          else
            (*selems_of_domain_type[CDOM_POS]).SetBitAtomic(elnr);
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
          cut_neighboring_node[NT_VERTEX]->SetBitAtomic(node);

        if(ma -> GetDimension() >=2){
          nodenums = ma->GetElEdges(elid);
          for (int node : nodenums)
            cut_neighboring_node[NT_EDGE]->SetBitAtomic(node);
        }
        if (ma->GetDimension() == 3)
        {
          nodenums = ma->GetElFaces(elid.Nr());
          for (int node : nodenums)
            cut_neighboring_node[NT_FACE]->SetBitAtomic(node);
        }
        cut_neighboring_node[NT_ELEMENT]->SetBitAtomic(elnr);
      }
    });

    for (NODE_TYPE nt : {NT_VERTEX,NT_EDGE,NT_FACE,NT_CELL})
      *(dom_of_node[nt]) = IF;

    Array<bool> neg_vertex;
    Array<bool> neg_edge;
    Array<bool> neg_face;
    neg_vertex.SetSize (ma->GetNV());
    if (ma->GetDimension() >= 2)
        neg_edge.SetSize (ma->GetNEdges());
    if (ma->GetDimension() == 3)
      neg_face.SetSize (ma->GetNFaces());
    neg_vertex = false;
    neg_edge = false;
    neg_face = false;
      
    ma->IterateElements(VOL, lh,[&](auto el, LocalHeap& mlh) {
      if (elems_of_domain_type[CDOM_NEG]->Test(el.Nr()))
      {
        for (auto v : el.Vertices())
          neg_vertex[v] = true;
        if (ma->GetDimension() >= 2)
          for (auto e : el.Edges())
            neg_edge[e] = true;
        if (ma->GetDimension() == 3)
          for (auto f : el.Faces())
            neg_face[f] = true;
      }
    });
      
    if (ma->GetCommunicator().Size() > 1)
    {
      ma->AllReduceNodalData (NT_VERTEX, neg_vertex, NG_MPI_LOR);
      ma->AllReduceNodalData (NT_EDGE, neg_edge, NG_MPI_LOR);
      ma->AllReduceNodalData (NT_FACE, neg_face, NG_MPI_LOR);
    }
    
    if (ma->GetDimension() == 3) {
      for (int i = 0; i < ma->GetNV(); i++)
        (*dom_of_node[NT_VERTEX])[i] = neg_vertex[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNEdges(); i++)
        (*dom_of_node[NT_EDGE])[i] = neg_edge[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNFaces(); i++)
        (*dom_of_node[NT_FACE])[i] = neg_face[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNE(); i++)
        (*dom_of_node[NT_CELL])[i] = (*cut_ratio_of_element[VOL])(i) > 0.5 ? NEG : POS;
    }
    else if(ma -> GetDimension() == 2){
      for (int i = 0; i < ma->GetNV(); i++)
        (*dom_of_node[NT_VERTEX])[i] = neg_vertex[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNEdges(); i++)
        (*dom_of_node[NT_EDGE])[i] = neg_edge[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNE(); i++)
        (*dom_of_node[NT_FACE])[i] = (*cut_ratio_of_element[VOL])(i) > 0.5 ? NEG : POS;
    }
    else if(ma -> GetDimension() == 1){
      for (int i = 0; i < ma->GetNV(); i++)
        (*dom_of_node[NT_VERTEX])[i] = neg_vertex[i] ? NEG : POS;
      for (int i = 0; i < ma->GetNEdges(); i++)
        (*dom_of_node[NT_EDGE])[i] = (*cut_ratio_of_element[VOL])(i) > 0.5 ? NEG : POS;
    }
    /* old        
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
    */
        
  }

  shared_ptr<BitArray> CutInformation::GetElementsWithThresholdContribution(DOMAIN_TYPE dt,
                                                                              double threshold,
                                                                              VorB vb){
    int ne = ma->GetNE(vb);
    shared_ptr<BitArray> elems_with_threshold = make_shared<BitArray>(ne);
    elems_with_threshold->Clear();
    if (dt == POS)
        threshold = 1 - threshold;

    LocalHeap dummy_lh (1000, "GetElementsWithThresholdContribution-heap", true);
    IterateRange
      (ne, dummy_lh,
      [&] (int elnr, LocalHeap & lh)
    {
      if (dt == NEG && (*cut_ratio_of_element[vb])(elnr) >= threshold)
        (*elems_with_threshold).SetBitAtomic(elnr);
      else if (dt == POS && (*cut_ratio_of_element[vb])(elnr) <= threshold)
        (*elems_with_threshold).SetBitAtomic(elnr);
    });

    return elems_with_threshold;
  }

  MultiLevelsetCutInformation::MultiLevelsetCutInformation (shared_ptr<MeshAccess> ama,
                                                            const Array<shared_ptr<GridFunction>> & lsets_in)
    : ma(ama), lsets(lsets_in)
  {
    ;
  }

  void MultiLevelsetCutInformation::Update(const Array<shared_ptr<GridFunction>> & lsets_in, LocalHeap & lh)
  {
    for(int i=0; i < lsets_in.Size(); i++)
      this->lsets[i]->GetVectorPtr()->Set(1.0, lsets_in[i]->GetVector());

    for (auto tup : collect_elements_with_contribution)
      UpdateElementsWithContribution(get<0>(tup), get<1>(tup), get<2>(tup), lh);

    for (auto tup : collect_elements_of_domain_type)
      UpdateElementsOfDomainType(get<0>(tup), get<1>(tup), get<2>(tup), lh);

  }

  void MultiLevelsetCutInformation::UpdateElementsWithContribution(const shared_ptr<BitArray> & elems_of_domain_type, 
                                                                   const Array<Array<DOMAIN_TYPE>> & cdt, 
                                                                   VorB vb, 
                                                                   LocalHeap & lh) const
  {

    LevelsetIntegrationDomain lsetintdom(lsets, cdt);
    elems_of_domain_type->Clear();
    
    int ne = ma->GetNE(vb);
    IterateRange
      (ne, lh,
       [&] (int elnr, LocalHeap & lh)
       {
         ElementId ei = ElementId(vb,elnr);
         ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

         const IntegrationRule * ir_np;
         Array<double> wei_arr;
         tie (ir_np, wei_arr) = CreateCutIntegrationRule(lsetintdom, eltrans, lh);
         double part_vol = 0;
         if (ir_np)
           for (auto w : wei_arr) 
             part_vol += w; 
         if (part_vol > 0)
           elems_of_domain_type->SetBitAtomic(elnr);
         
       });
  }

  shared_ptr<BitArray> MultiLevelsetCutInformation::GetElementsWithContribution(const Array<Array<DOMAIN_TYPE>> & cdt,
                                                                                VorB vb, 
                                                                                LocalHeap & lh)
  {
    // Check whether BitArray with requested cdt has already been computed
    // and return that BitArray if found
    for (auto tup : collect_elements_with_contribution)
      if (CombinedDomainTypesEqual(get<1>(tup), cdt))
        return get<0>(tup);

    // Compute the new BitArray
    shared_ptr<BitArray> elems_of_domain_type = make_shared<BitArray>(ma->GetNE(vb));
    UpdateElementsWithContribution(elems_of_domain_type, cdt, vb, lh);

    // Add BitArray to instance collection
    collect_elements_with_contribution.push_back(make_tuple(elems_of_domain_type, cdt, vb));

    return elems_of_domain_type;
  }
  
  void MultiLevelsetCutInformation::UpdateElementsOfDomainType(const shared_ptr<BitArray> & elems_of_domain_type, 
                                                               const Array<Array<DOMAIN_TYPE>> & cdt, 
                                                               VorB vb, 
                                                               LocalHeap & lh) const
  {
    LevelsetIntegrationDomain lsetintdom(lsets, cdt);
    elems_of_domain_type->Clear();
    
    int ne = ma->GetNE(vb);
    IterateRange
      (ne, lh,
       [&] (int elnr, LocalHeap & lh)
       {
         ElementId ei = ElementId(vb,elnr);
         // Ngs_Element ngel = ma->GetElement(ei);
         // ELEMENT_TYPE eltype = ngel.GetType();
         ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

         Array<DofId> dnums(0,lh);
         int M =lsetintdom.GetLevelsetGFs().Size();
         Array<DOMAIN_TYPE> dts_el(M,lh);
         for (int j = 0; j < M; j++)
         {
           auto gflset = lsetintdom.GetLevelsetGFs()[j];
           gflset->GetFESpace()->GetDofNrs(eltrans.GetElementId(), dnums);
           FlatVector<> elvec(dnums.Size(), lh);
           gflset->GetVector().GetIndirect(dnums, elvec);

           dts_el[j] = CheckIfStraightCut(elvec);
         }

         for (int i = 0; i < lsetintdom.GetDomainTypes().Size(); i++)
         {
           bool matching = true;
           for (int j = 0; j < M; j++)
             if (lsetintdom.GetDomainTypes()[i][j] != dts_el[j])
               matching = false;

           if (matching)
             elems_of_domain_type->SetBitAtomic(elnr);
         }
       });
  }

  shared_ptr<BitArray> MultiLevelsetCutInformation::GetElementsOfDomainType(const Array<Array<DOMAIN_TYPE>> & cdt, 
                                                                            VorB vb, 
                                                                            LocalHeap & lh)
  {
    // Check whether BitArray with requested cdt has already been computed
    // and return that BitArray if found
    for (auto tup : collect_elements_of_domain_type)
      if (CombinedDomainTypesEqual(get<1>(tup), cdt))
        return get<0>(tup);

    // Compute the new BitArray
    shared_ptr<BitArray> elems_of_domain_type = make_shared<BitArray>(ma->GetNE(vb));
    UpdateElementsOfDomainType(elems_of_domain_type, cdt, vb, lh);
    
    // Add BitArray to instance collection
    collect_elements_of_domain_type.push_back(make_tuple(elems_of_domain_type, cdt, vb));

    return elems_of_domain_type;
  }
  

  bool MultiLevelsetCutInformation::CombinedDomainTypesEqual(const Array<Array<DOMAIN_TYPE>> & cdta,
                                                             const Array<Array<DOMAIN_TYPE>> & cdtb) const
  {
    int len_a = cdta.Size();
    int len_b = cdtb.Size();
    if (len_a != len_b)
        return false;

    unordered_map<string, int> count;
    string dtt_encoded;
    
    for (const Array<DOMAIN_TYPE> & dtt : cdta)
    {
      dtt_encoded.clear();
      for(const DOMAIN_TYPE & dt : dtt)
      {
        if (dt == NEG)
          dtt_encoded.push_back('n');
        else if (dt == IF)
          dtt_encoded.push_back('i');
        else if (dt == POS)
          dtt_encoded.push_back('p');
        else
          throw Exception("Unable to encode domain domain type");
      }
      count[dtt_encoded]++;
    }

    for (const Array<DOMAIN_TYPE> & dtt : cdtb)
    {
      dtt_encoded.clear();
      for(const DOMAIN_TYPE & dt : dtt)
      {
        if (dt == NEG)
          dtt_encoded.push_back('n');
        else if (dt == IF)
          dtt_encoded.push_back('i');
        else if (dt == POS)
          dtt_encoded.push_back('p');
        else
          throw Exception("Unable to encode domain domain type");
      }
      count[dtt_encoded]--;
    }

    for (auto i : count)
      if (i.second != 0)
        return false;

    return true;
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
        fine_facet.SetBitAtomic(fanums[j]);
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
        }
        
        bool a_left = a->Test(elnums[0]);
        bool a_right = elnums.Size() > 1 ? a->Test(elnums[1]) : bound_val_a;
        bool b_left = b->Test(elnums[0]);
        bool b_right = elnums.Size() > 1 ? b->Test(elnums[1]) : bound_val_b;

        if (ask_and)
        {
          if ((a_left && b_right) || (a_right && b_left))
            ret->SetBitAtomic(facnr);
        }
        else
        {
          if ((a_left || b_right) || (a_right || b_left))
            ret->SetBitAtomic(facnr);
        }
      }
    });
    return ret;
  }

  shared_ptr<BitArray> GetElementsWithNeighborFacets(shared_ptr<MeshAccess> ma,
                                                     shared_ptr<BitArray> a,
                                                     LocalHeap & lh)
  {
    if (ma->GetCommunicator().Size() > 1)
      throw Exception("GetElementsWithNeighborFacets:: No Ghost-Markers for MPI yet");

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
          ret->SetBitAtomic(elnr);
      }
    });
    return ret;
  }

  shared_ptr<BitArray> GetElementsWithSharedVertex(shared_ptr<MeshAccess> ma, shared_ptr<BitArray> a, LocalHeap & lh)
  {
    if (ma->GetCommunicator().Size() > 1)
      throw Exception("GetElementsWithNeighborFacets:: No GetElementsWithSharedVertex for MPI yet");
    int ne = ma->GetNE();
    shared_ptr<BitArray> ret = make_shared<BitArray> (ne);
    ret->Clear();
    IterateRange
      (ne, lh,
      [&] (int elnr, LocalHeap & lh)
    {
        if(a->Test(elnr)){
            ElementId elid(VOL,elnr);
            Array<int> nodenums(0,lh);
            nodenums = ma->GetElVertices(elid);
            for (int node : nodenums) {
              /*
                Array<int> elnums(0,lh);
                ma->GetVertexElements(node, elnums);
                for (auto elnr : elnums)
                    ret->SetBitAtomic(elnr);
              */
              for (auto elnr : ma->GetVertexElements(node))
                ret->SetBitAtomic(elnr);
              
            }
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
          ret->SetBitAtomic(dof);
      }
    });
    return ret;
  }

  shared_ptr<BitArray> GetDofsOfFacets(shared_ptr<FESpace> fes,
                                       shared_ptr<BitArray> a,
                                       LocalHeap & lh)
  {
    int nf = fes->GetMeshAccess()->GetNFacets();
    int ndof = fes->GetNDof();
    shared_ptr<BitArray> ret = make_shared<BitArray> (ndof);
    ret->Clear();

    IterateRange
      (nf, lh,
      [&] (int fanr, LocalHeap & lh)
    {
      NodeId nodeid(NT_FACET,fanr);
      if (a->Test(fanr))
      {
        Array<int> dnums(0,lh);
        fes->GetDofNrs(nodeid,dnums);
        for (auto dof : dnums)
          ret->SetBitAtomic(dof);
      }
    });
    return ret;
  }


}
