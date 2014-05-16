#ifndef FILE_OUTPUT_HPP
#define FILE_OUTPUT_HPP

#include <ngstd.hpp> // for Array
#include <fem.hpp>   // for ScalarFiniteElement
#include <solve.hpp>

// #include "xintegration.hpp"
#include "../spacetime/spacetimefespace.hpp"
#include "../xfem/xfemIntegrators.hpp"
#include "../xfem/stxfemIntegrators.hpp"
#include "../xfem/setvaluesx.hpp"
#include "../utils/error.hpp"

using namespace ngsolve;

/// from ngsolve

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"

namespace ngfem
{

  template<int D>
  void DoSpecialOutput (GridFunction * gfu, 
                        SolutionCoefficients<D> & solcoef, 
                        int subdivision, 
                        LocalHeap & lh);


  void OutputMeshOnly (const MeshAccess & ma, LocalHeap & lh);

}
#endif
