/*********************************************************************/
/* File:   lsetrefine.hpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. Jul. 2015                                             */
/*********************************************************************/
#pragma once

#include <solve.hpp>

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

namespace ngcomp
{ 

  void RefineAtLevelSet (shared_ptr<GridFunction> gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh);
  
}
