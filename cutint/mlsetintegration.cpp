#include "mlsetintegration.hpp"
#include "straightcutrule.hpp"

namespace xintegration
{
  tuple<const IntegrationRule *, Array<double> > CreateMultiLevelsetCutIntegrationRule(const LevelsetIntegrationDomain & lsetintdom,
                                                                                       const ElementTransformation & trafo,
                                                                                       LocalHeap & lh)
  {
    bool debug_out = false;

    ArrayMem<const IntegrationRule *, 20> ir_parts; //static/dynamic memory allocation
      
    const Array<shared_ptr<GridFunction>> & gflsets (lsetintdom.GetLevelsetGFs());
    int intorder = lsetintdom.GetIntegrationOrder();
    int time_intorder = lsetintdom.GetTimeIntegrationOrder();
    SWAP_DIMENSIONS_POLICY quad_dir_policy = lsetintdom.GetSwapDimensionPolicy();
    
    int M = gflsets.Size();
    if (time_intorder >= 0)
      throw Exception("no space-time integration with multiple level sets, yet.");
    
    ArrayMem<DofId,10> dnums;
    gflsets[0]->GetFESpace()->GetDofNrs(trafo.GetElementId(), dnums);
    FlatVector<> elvec(dnums.Size(), lh);
    FlatMatrix<> elvecs(dnums.Size(), M, lh);

    ArrayMem<DOMAIN_TYPE,10> lset_dts(M); //<- domain types corresponding to levelset

    // Check if element is relevant for the given integration domain
    for (int i = 0; i < M; i++)
    {
      gflsets[i]->GetVector().GetIndirect(dnums, elvec);
      elvecs.Col(i) = elvec;
      lset_dts[i] = CheckIfStraightCut(elvec);
    }

    if (debug_out)
    {
      cout << IM(3) << "#########################################################################" << endl;
      cout << IM(3) << "#########################################################################" << endl << endl;
      cout << IM(3) << "intorder : " << intorder << endl;
      cout << IM(3) << "currently on element with ID: " << trafo.GetElementId() << endl << endl;
      cout << IM(3) << "level set vertex values (one column per lset):\n" << elvecs << endl;
      cout << IM(3) << "domain types current of element (corresp. to mlset): \n" << lset_dts << endl;
    }
    
    for (const Array<DOMAIN_TYPE> & dts : lsetintdom.GetDomainTypes())
    {
      bool compatible = true;
      bool cut_element = false;

      const IntegrationRule* ir = nullptr;
    
      for (int i = 0; i < M; i++)
      {
        if ((lset_dts[i] != IF) && (lset_dts[i] != dts[i])) compatible = false; //loop could break here now
        else if (lset_dts[i] == IF) cut_element = true;
      }

      if (debug_out)
      {
        cout << IM(3) << "###################################### " << endl;
        cout << IM(3) << "domain types for integration: \n" << dts << endl;
      }

      ArrayMem<shared_ptr<GridFunction>,10> condense_gflsets;
      ArrayMem<DOMAIN_TYPE,10> condense_dts;// only for debugging
      ArrayMem<DOMAIN_TYPE,10> condense_target_dts;
      ArrayMem<int,10> condense_id_to_full_id;

      if (!compatible)
      {
        if (debug_out)
        {
          cout << IM(3) << "levelset configuration not compatible with integration domain: skipping" << endl;
          cout << IM(3) << "----------------------------------------------------------------------" << endl << endl;
        }
        continue;
      }
 
    
      if (!cut_element)
      {
        if (debug_out)
        {
          cout << IM(3) << "relevant, uncut element: standard integration rule" << endl;
          cout << IM(3) << "--------------------------------------------------" << endl << endl;
        }      
        // Uncut elements simply return the standard integration rule
        ir = & (SelectIntegrationRule (trafo.GetElementType(), intorder));
        ir_parts.Append(ir);
      }

      else
      {
        if (debug_out)
          cout << IM(3) << "This element is cut and relevant ";

        // We reduce the domain type array to the necessary domain types and corresponding level sets 
        for (int i = 0; i < M; i++)
        {
          if ((lset_dts[i] == IF) || (lset_dts[i] != dts[i]))
          {
            condense_dts.Append(lset_dts[i]); // only for debugging
            condense_target_dts.Append(dts[i]);
            condense_gflsets.Append(gflsets[i]);
            condense_id_to_full_id.Append(i);
          }
        }
        if (debug_out)
          cout << IM(3) << "condensed ids:\n" << condense_id_to_full_id << endl;// only for debugging
        
        // If the element is cut by only one levelset, then we can use standard cut integration rules
        if (condense_id_to_full_id.Size() == 1)
        {      
          if (debug_out)
          {
            cout << IM(3) << "relevant element, cut by only one levelset: standard cut integration rule" << endl;
            cout << IM(3) << "-------------------------------------------------------------------------" << endl << endl;
          }
          // condense_gflsets[0]->GetVector().GetIndirect(dnums, elvec);
          int j = condense_id_to_full_id[0];
          elvec = elvecs.Col(j);
          if (debug_out)
          {
            cout << IM(3) << elvec << endl;
          }
          const IntegrationRule * ir = StraightCutIntegrationRule(elvec,
                                                                  trafo, dts[j], intorder, quad_dir_policy, lh);
          if(ir != nullptr) 
            ir_parts.Append(ir);
        }
        else 
        { 
          // An element cut by multiple level sets
          // At least for debugging, this else case should also be able to handle the "if case"
          FlatMatrix<> elvecs_cond(dnums.Size(), condense_id_to_full_id.Size(), lh);
          for (int i = 0; i < condense_id_to_full_id.Size(); i++)
            elvecs_cond.Col(i) = elvecs.Col(condense_id_to_full_id[i]);

          const IntegrationRule * ir = StraightCutsIntegrationRule(elvecs_cond, trafo, condense_target_dts, intorder, quad_dir_policy, lh);
          if (debug_out)
          {
            cout << IM(3) << "relevant element, custom integration rule" << endl;
            cout << IM(3) << "-------------------------------------------------------------------------" << endl << endl;
          }
          if(ir != nullptr)
            ir_parts.Append(ir);
        }
      }
    }
    IntegrationRule * sum_ir = nullptr;

    int size = 0;
    for ( auto part : ir_parts )
      size += part->Size();

    // ArrayMem<double,100> sum_wei_arr(size);
    Array<double> sum_wei_arr(size);
    
    if (ir_parts.Size() > 0) //collect parts:
    {
      sum_ir = new (lh) IntegrationRule(size , lh);
      int cnt = 0;
      for ( auto part : ir_parts )
        for ( auto ip : *part )
        {
          (*sum_ir)[cnt] = ip;
          sum_wei_arr[cnt] = ip.Weight();
          cnt++;
        }
    }

    if (debug_out)
    {
      double sum = 0;
      for (auto w : sum_wei_arr)
      {
        sum += w;
      }
      cout << IM(3) << "Sum of weights on element: " << sum << endl << endl;
    }
    
    return make_tuple(sum_ir, sum_wei_arr);
  }
  
  /// Mlset-integrationRule Factory
  tuple<const IntegrationRule *, Array<double> > CreateCutIntegrationRule(Array<shared_ptr<GridFunction>> & gflsets,
                                                                          const ElementTransformation & trafo,
                                                                          Array<DOMAIN_TYPE> & dts_in,
                                                                          int intorder,
                                                                          int time_intorder,
                                                                          LocalHeap & lh,
                                                                          SWAP_DIMENSIONS_POLICY quad_dir_policy)
  {
    Array<Array<DOMAIN_TYPE>> dts(1);
    dts[0] = dts_in;
    LevelsetIntegrationDomain lsetintdom(gflsets,dts,intorder,time_intorder,0,quad_dir_policy);
    return CreateMultiLevelsetCutIntegrationRule(lsetintdom,trafo,lh);
  }

} // end of namespace
