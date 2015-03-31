#include "output.hpp"

#include<stdlib.h>

namespace ngfem
{

  template<int D>
  void DoSpecialOutput (shared_ptr<GridFunction> gfu, 
                        SolutionCoefficients<D> & solcoef, 
                        int subdivision, 
                        Flags & flags,
                        LocalHeap & lh)
  {
    throw Exception("empty...");
  }

  template void DoSpecialOutput<2>(shared_ptr<GridFunction> gfu, SolutionCoefficients<2> & solcoef, 
                                   int subdivision, Flags & flags, LocalHeap & lh);
  template void DoSpecialOutput<3>(shared_ptr<GridFunction> gfu, SolutionCoefficients<3> & solcoef, 
                                   int subdivision, Flags & flags,LocalHeap & lh);


  void OutputMeshOnly (shared_ptr<MeshAccess> ma, LocalHeap & lh)
  {
    throw Exception("empty...");
  }

}
