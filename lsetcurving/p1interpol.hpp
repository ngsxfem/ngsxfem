/*********************************************************************/
/* File:   p1interpol.cpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   July 17th 2015                                            */
/*********************************************************************/
#pragma once

#include <ngstd.hpp> // for Array
#include <fem.hpp>   // for ScalarFiniteElement
#include <comp.hpp> 
#include <solve.hpp>

using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{ 

/* ---------------------------------------- 
   Interpolate a coefficient function 
   or an h1ho function into the space of
   piecewise linears
   ---------------------------------------- */
  class InterpolateP1
  {
  protected:
    shared_ptr<MeshAccess> ma = nullptr;
    shared_ptr<CoefficientFunction> coef;
    shared_ptr<GridFunction> gf;
    shared_ptr<GridFunction> gf_p1;
  public:
    InterpolateP1 (shared_ptr<CoefficientFunction> a_coef, shared_ptr<GridFunction> a_gf_p1);
    InterpolateP1 (shared_ptr<GridFunction> a_gf, shared_ptr<GridFunction> a_gf_p1);
    void Do (LocalHeap & lh);
  };


  class NumProcInterpolateP1 : public NumProc
  {
  protected:
    shared_ptr<InterpolateP1> interpol = nullptr;
  public:
    NumProcInterpolateP1 (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcInterpolateP1() { }
    virtual string GetClassName () const
    {
      return "NumProcInterpolateP1";
    }
    virtual void Do (LocalHeap & lh);
  };
  
}
