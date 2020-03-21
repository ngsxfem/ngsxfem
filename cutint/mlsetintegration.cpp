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

    Array<DOMAIN_TYPE> & dts (lsetintdom.GetDomainTypes()[0]);
    if (lsetintdom.GetDomainTypes().Size() > 1)
      throw Exception("not gathering the integration rules from the domain_type_arrays together, yet");
    // TODO
    // const Array<Array<DOMAIN_TYPE>> & list_dts (lsetintdom.GetDomainTypes()); //TODO -> should yield in a loop...
    // for (auto dts: list_dts) ... gather rules..
    

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
    bool compatible = true;
    bool cut_element = false;

    const IntegrationRule* ir = nullptr;

    Array<shared_ptr<GridFunction>> condense_gflsets;
    Array<DOMAIN_TYPE> condense_dts;// only for debugging
    Array<DOMAIN_TYPE> condense_target_dts;

    // Check if element is relevant for the given integration domain
    for (int i = 0; i < M; i++)
    {
      gflsets[i]->GetVector().GetIndirect(dnums, elvec);
      elvecs.Col(i) = elvec;
      lset_dts[i] = CheckIfStraightCut(elvec);
      if ((lset_dts[i] != IF) && (lset_dts[i] != dts[i])) compatible = false; //loop could break here now
      else if (lset_dts[i] == IF) cut_element = true;
    }

    if (debug_out)
    {
      cout << "currently on element with ID: " << trafo.GetElementId() << endl << endl;
      cout << "level set vertex values (one column per lset):\n" << elvecs << endl;
      cout << "domain types current of element (corresp. to mlset): \n" << lset_dts << endl;
      cout << "domain types for integration: \n" << dts << endl;
    }
    
    if (!compatible)
    {
      if (debug_out)
      {
        cout << "levelset configuration not compatible with integration domain: skipping" << endl;
        cout << "----------------------------------------------------------------------" << endl << endl;
      }
      return make_tuple(nullptr, Array<double>());
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
      Array<double> wei_arr (ir->Size());
      for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
      return make_tuple(ir, wei_arr);
    }

    else
    {
      if (debug_out)
        cout << "This element is cut and relevant";

      // We reduce the domain type array to the necessary domain types and corresponding level sets 
      for (int i = 0; i < M; i++)
      {
        if ((lset_dts[i] == IF) || (lset_dts[i] != dts[i]))
        {
          condense_dts.Append(lset_dts[i]); // only for debugging
          condense_target_dts.Append(dts[i]);
          condense_gflsets.Append(gflsets[i]);
        }
      }
      if (debug_out)
      {
        cout << "domain type after condensing\n" << condense_dts << endl;// only for debugging
        cout << "target type after condensing\n" << condense_target_dts << endl;
      }
      // If the element is cut by only one levelset, then we can use standard cut integration rules
      if (condense_dts.Size() == 1)
      {      
        if (debug_out)
        {
          cout << "relevant element, cut by only one levelset: standard cut integration rule" << endl;
          cout << "-------------------------------------------------------------------------" << endl << endl;
        }
        condense_gflsets[0]->GetVector().GetIndirect(dnums, elvec);
        const IntegrationRule * ir = StraightCutIntegrationRule(elvec, trafo, condense_target_dts[0], intorder, quad_dir_policy, lh);
        if(ir != nullptr) 
        {
          Array<double> wei_arr (ir->Size());
          for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
          return make_tuple(ir, wei_arr);
        }
      }
      else 
      { 
        // An element cut by multiple level sets
        // At least for debugging, this else case should also be able to handle the "if case"
        FlatMatrix<> elvecs_cond(dnums.Size(), condense_dts.Size(), lh);
        for (int i = 0; i < condense_dts.Size(); i++)
        {
          condense_gflsets[i]->GetVector().GetIndirect(dnums, elvec);
          elvecs_cond.Col(i) = elvec;
        }

        const IntegrationRule * ir = StraightCutsIntegrationRule(elvecs_cond, trafo, condense_target_dts, intorder, quad_dir_policy, lh);
        if(ir != nullptr)
        {
          Array<double> wei_arr (ir->Size());
          for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
          return make_tuple(ir, wei_arr);
        }
      }
    }

    return make_tuple(nullptr, Array<double>());
      
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
