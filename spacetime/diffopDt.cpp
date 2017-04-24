#define FILE_DIFFOPDT_CPP
#include "diffopDt.hpp"
#include <diffop_impl.hpp>
#include "myElement.hpp"


namespace ngfem
{

  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDt::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {

      const SpaceTimeFE & scafe =
              dynamic_cast<const SpaceTimeFE & > (bfel);
      const int ndof = scafe.GetNDof();

      FlatVector<> dtshape (ndof,lh);
      IntegrationPoint ip(mip.IP());
      scafe.CalcDtShape(ip,dtshape);
      mat = 0.0;
      mat.Row(0) = dtshape;


    }

  template class T_DifferentialOperator<DiffOpDt>;


}


