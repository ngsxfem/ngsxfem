/*********************************************************************/
/* File:   lsetrefine.cpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. July 2015                                             */
/*********************************************************************/

#include "lsetrefine.hpp"

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
      // element only measure error if "at the interface"
      Array<int> dnums_lset_p1;
      gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dnums_lset_p1);
      FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),lh);
      gf_lset_p1->GetVector().GetIndirect(dnums_lset_p1,lset_vals_p1);
            
      bool has_pos = (lset_vals_p1[D] > lower_lset_bound);
      bool has_neg = (lset_vals_p1[D] < upper_lset_bound);
      for (int d = 0; d < D; ++d)
      {
        if (lset_vals_p1[d] > lower_lset_bound)
          has_pos = true;
        if (lset_vals_p1[d] < upper_lset_bound)
          has_neg = true;
      }
      if (has_pos && has_neg)
        Ng_SetRefinementFlag (elnr+1, 1);
      else
        Ng_SetRefinementFlag (elnr+1, 0);
    }
      
  }



  NumProcLevelSetRefine::NumProcLevelSetRefine (shared_ptr<PDE> apde, const Flags & flags)
  {
    lower_lset_bound = flags.GetNumFlag("lower_lset_bound",0.0);
    upper_lset_bound = flags.GetNumFlag("upper_lset_bound",0.0);
    gf_lset_p1 = apde->GetGridFunction(flags.GetStringFlag("levelset","gf_lset_p1"));
  }
  
  void NumProcLevelSetRefine::Do (LocalHeap & lh)
  {
    RefineAtLevelSet(gf_lset_p1, lower_lset_bound, upper_lset_bound, lh);
  }


}

static RegisterNumProc<NumProcLevelSetRefine> nplsetrefine("levelsetrefine");

