/*********************************************************************/
/* File:   calcgeomerrors.hpp                                        */
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

  void PrintConvergenceTable(const Array<double> & tab, string label="");

  struct StatisticContainer
  {
    Array<double> ErrorL2Norm;
    Array<double> ErrorL1Norm;
    Array<double> ErrorMaxNorm;
    Array<double> ErrorMisc;
  };
  
  template <int D>
  void CalcDistances (shared_ptr<CoefficientFunction> gf_lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, StatisticContainer & cont, LocalHeap & lh, double define_threshold = -1.0, bool abs_ref_threshold = false);

  template<int D>
  void CalcDeformationError (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn, StatisticContainer & cont, LocalHeap & lh, double, double);

  
}
