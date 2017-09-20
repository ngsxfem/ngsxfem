#pragma once

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
#include <comp.hpp>  // for Gridfunction, Coeff...
#include "calcpointshift.hpp"


namespace ngcomp
{

  void ProjectShift (shared_ptr<GridFunction> lset_ho, shared_ptr<GridFunction> lset_p1,
                     shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn,
                     shared_ptr<BitArray> ba,
                     shared_ptr<CoefficientFunction> blending,
                     double lower_lset_bound, double upper_lset_bound, double threshold,
                     LocalHeap & lh);

}
