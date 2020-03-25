#include "mlsetintegration.hpp"
#include "straightcutrule.hpp"

namespace xintegration
{
  using ngfem::INT;


  tuple<const IntegrationRule *, Array<double> > CreateMultiLevelsetCutIntegrationRule(const LevelsetIntegrationDomain & lsetintdom,
                                                                                       const ElementTransformation & trafo,
                                                                                       LocalHeap & lh)
  {
    bool debug_out = false;

    Array<double> sum_wei_arr (0);
    IntegrationRule * sum_ir = new (lh) IntegrationRule(0 , lh);
      
    const Array<shared_ptr<GridFunction>> & gflsets (lsetintdom.GetLevelsetGFs());
    int intorder = lsetintdom.GetIntegrationOrder();
    int time_intorder = lsetintdom.GetTimeIntegrationOrder();
    SWAP_DIMENSIONS_POLICY quad_dir_policy = lsetintdom.GetSwapDimensionPolicy();
    
    int M = gflsets.Size();
    if (time_intorder >= 0)
      throw Exception("no space-time integration with multiple level sets, yet.");
    
    Array<DofId> dnums(0,lh);
    gflsets[0]->GetFESpace()->GetDofNrs(trafo.GetElementId(), dnums);
    FlatVector<> elvec(dnums.Size(), lh);
    FlatMatrix<> elvecs(dnums.Size(), M, lh);

    Array<DOMAIN_TYPE> lset_dts(M); //<- domain types corresponding to levelset

    // Check if element is relevant for the given integration domain
    for (int i = 0; i < M; i++)
    {
      gflsets[i]->GetVector().GetIndirect(dnums, elvec);
      elvecs.Col(i) = elvec;
      lset_dts[i] = CheckIfStraightCut(elvec);
    }

    if (debug_out)
    {
      cout << "#########################################################################" << endl;
      cout << "#########################################################################" << endl << endl;
      cout << "intorder : " << intorder << endl;
      cout << "currently on element with ID: " << trafo.GetElementId() << endl << endl;
      cout << "level set vertex values (one column per lset):\n" << elvecs << endl;
      cout << "domain types current of element (corresp. to mlset): \n" << lset_dts << endl;
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
        cout << "###################################### " << endl;
        cout << "domain types for integration: \n" << dts << endl;
      }

      Array<shared_ptr<GridFunction>> condense_gflsets;
      Array<DOMAIN_TYPE> condense_dts;// only for debugging
      Array<DOMAIN_TYPE> condense_target_dts;
      Array<int> condense_id_to_full_id(0);

      if (!compatible)
      {
        if (debug_out)
        {
          cout << "levelset configuration not compatible with integration domain: skipping" << endl;
          cout << "----------------------------------------------------------------------" << endl << endl;
        }
        continue;
      }
 
    
      if (!cut_element)
      {
        if (debug_out)
        {
          cout << "relevant, uncut element: standard integration rule" << endl;
          cout << "--------------------------------------------------" << endl << endl;
        }      
        // Uncut elements simply return the standard integration rule
        ir = & (SelectIntegrationRule (trafo.GetElementType(), intorder));
        for(int i=0; i< ir->Size(); i++)
        {
          sum_ir->Append((*ir)[i]);
          sum_wei_arr.Append((*ir)[i].Weight());
        }
        // Array<double> wei_arr (ir->Size());
        // for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
        // return make_tuple(ir, wei_arr);
      }

      else
      {
        if (debug_out)
          cout << "This element is cut and relevant ";

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
          cout << "condensed ids:\n" << condense_id_to_full_id << endl;// only for debugging
        
        // If the element is cut by only one levelset, then we can use standard cut integration rules
        if (condense_id_to_full_id.Size() == 1)
        {      
          if (debug_out)
          {
            cout << "relevant element, cut by only one levelset: standard cut integration rule" << endl;
            cout << "-------------------------------------------------------------------------" << endl << endl;
          }
          // condense_gflsets[0]->GetVector().GetIndirect(dnums, elvec);
          int j = condense_id_to_full_id[0];
          elvec = elvecs.Col(j);
          if (debug_out)
          {
            cout << elvec << endl;
          }
          const IntegrationRule * ir = StraightCutIntegrationRule(elvec,
                                                                  trafo, dts[j], intorder, quad_dir_policy, lh);
          if(ir != nullptr) 
          {
            // Array<double> wei_arr (ir->Size());
            for(int i=0; i< ir->Size(); i++)
            {
              sum_ir->Append((*ir)[i]);
              sum_wei_arr.Append((*ir)[i].Weight());
            }
          }
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
            cout << "relevant element, custom integration rule" << endl;
            cout << "-------------------------------------------------------------------------" << endl << endl;
          }
          if(ir != nullptr)
          {
            for(int i=0; i< ir->Size(); i++)
            {
              sum_ir->Append((*ir)[i]);
              sum_wei_arr.Append((*ir)[i].Weight());
            }
          }
        }
      }
    }

    if (debug_out)
    {
      double sum = 0;
      for (auto w : sum_wei_arr)
      {
        sum += w;
      }
      cout << "Sum of weigts on element: " << sum << endl << endl;
    }
    if (sum_ir->Size() == 0)
      sum_ir = nullptr;
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
