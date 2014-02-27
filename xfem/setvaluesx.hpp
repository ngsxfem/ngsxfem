#ifndef FILE_SETVALUESX_HPP
#define FILE_SETVALUESX_HPP

/// from ngsolve
#include <solve.hpp>    // provides FESpace, ...
#include <comp.hpp>    // provides FESpace, ...
#include <fem.hpp>

using namespace ngsolve;
// using namespace cutinfo;

typedef std::pair<double,double> TimeInterval;

namespace ngcomp
{

  template <int D, class SCAL>
  void SetValuesX (const Array<CoefficientFunction *> & acoefs,
                   const TimeInterval & ti,
                   GridFunction & bu,
                   bool bound,
                   LocalHeap & clh);

  
}    

#endif
