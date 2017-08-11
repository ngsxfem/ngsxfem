/*********************************************************************/
/* File:   lsetrefine.cpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. July 2015                                             */
/*********************************************************************/

#include "lsetrefine.hpp"
#include "calcpointshift.hpp"

namespace ngcomp
{ 

  void RefineAtLevelSet (shared_ptr<GridFunction> gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh){

    auto ma = gf_lset_p1->GetMeshAccess();
    const int D = ma->GetDimension();
    if (D == 3)
    {
      int nse = ma->GetNSE();
      for (int i = 0; i < nse; i++)
        Ng_SetSurfaceRefinementFlag (i+1, 0);
    }

    int ne=ma->GetNE();
      
    for (int elnr = 0; elnr < ne; ++elnr)
    {
      HeapReset hr(lh);
      // element only measure error if "at the interface"
      Array<int> dnums_lset_p1;
      gf_lset_p1->GetFESpace()->GetDofNrs(ElementId(VOL,elnr),dnums_lset_p1);
      FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),lh);
      gf_lset_p1->GetVector().GetIndirect(dnums_lset_p1,lset_vals_p1);

      if (ElementInRelevantBand(lset_vals_p1, lower_lset_bound, upper_lset_bound))
        Ng_SetRefinementFlag (elnr+1, 1);
      else
        Ng_SetRefinementFlag (elnr+1, 0);
    }
      
  }

}

