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
  };
  
  template <int D>
  void CalcDistances (shared_ptr<GridFunction> gf_lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, StatisticContainer & cont, LocalHeap & lh);

  class NumProcCalcErrors : public NumProc
  {
  protected:
    // double lower_lset_bound;
    // double upper_lset_bound;
    shared_ptr<CoefficientFunction> lset;
    shared_ptr<GridFunction> gf_lset_p1;
    shared_ptr<GridFunction> deform;
    StatisticContainer error_container;
  public:
    NumProcCalcErrors (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcCalcErrors() { }

    virtual string GetClassName () const
    {
      return "NumProcCalcErrors";
    }
    
    virtual void Do (LocalHeap & lh);
  };
  
}
