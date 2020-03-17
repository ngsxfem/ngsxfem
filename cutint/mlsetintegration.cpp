#include "mlsetintegration.hpp"
#include "straightcutrule.hpp"

namespace xintegration
{
  using ngfem::INT;


  /// Mlset-integrationRule Factory
  tuple<const IntegrationRule *, Array<double> > CreateCutIntegrationRule(Array<shared_ptr<GridFunction>> & gflsets,
                                                                          const ElementTransformation & trafo,
                                                                          Array<DOMAIN_TYPE> & dts,
                                                                          int intorder,
                                                                          int time_intorder,
                                                                          LocalHeap & lh,
                                                                          SWAP_DIMENSIONS_POLICY quad_dir_policy)
  {
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

    for (int i = 0; i < M; i++)
    {
      gflsets[i]->GetVector().GetIndirect(dnums, elvec);
      elvecs.Col(i) = elvec;
      lset_dts[i] = CheckIfStraightCut(elvec);
      if ((lset_dts[i] != IF) && (lset_dts[i] != dts[i])) compatible = false; //loop could break here now
      else if (lset_dts[i] == IF) cut_element = true;
    }
    
    cout << "currently on element with ID: " << trafo.GetElementId() << endl << endl;
    cout << "level set vertex values (one column per lset):\n" << elvecs << endl;
    cout << "domain types current of element (corresp. to mlset): \n" << lset_dts << endl;
    cout << "domain types for integration: \n" << dts << endl;

    if (!compatible)
    {
      cout << "levelset configuration not compatible with integration domain: skipping" << endl;
      cout << "----------------------------------------------------------------------" << endl << endl;
      return make_tuple(nullptr, Array<double>());
    }
 
    if (!cut_element)
    {
      cout << "relevant, uncut element: standard integration rule" << endl;
      cout << "--------------------------------------------------" << endl << endl;
      
      ir = & (SelectIntegrationRule (trafo.GetElementType(), intorder));
      Array<double> wei_arr (ir->Size());
      for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
      return make_tuple(ir, wei_arr);
    }
    else
    {
      cout << "This element is cut and relevant";

      //condense lset_dts
      for (int i = 0; i < M; i++)
      {
        if ((lset_dts[i] == IF) || (lset_dts[i] != dts[i]))
        {
          condense_dts.Append(lset_dts[i]); // only for debugging
          condense_target_dts.Append(dts[i]);
          condense_gflsets.Append(gflsets[i]);
        }
      }
      cout << "domain type after condensing\n" << condense_dts << endl;// only for debugging
      cout << "target type after condensing\n" << condense_target_dts << endl;

      // If only one cut remaining, use standard cut integration rule
      if (condense_dts.Size() == 1)
      {      
        cout << "relevant element, cut by only one levelset: standard cut integration rule" << endl;
        cout << "-------------------------------------------------------------------------" << endl << endl;

        condense_gflsets[0]->GetVector().GetIndirect(dnums, elvec);
        const IntegrationRule * ir = StraightCutIntegrationRule( elvec, trafo, condense_target_dts[0], intorder, quad_dir_policy, lh);
        if(ir != nullptr) 
        {
          Array<double> wei_arr (ir->Size());
          for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
          return make_tuple(ir, wei_arr);
        }
      }
      else
          cout << "Corner Case !!!" << endl;
          //throw Exception("not yet implemented");
    }

    cout << "further implementation is missing!" << endl;
    // throw Exception("not yet implemented");
    return make_tuple(nullptr, Array<double>());
  }


  // Compute relevant codim 

} // end of namespace
