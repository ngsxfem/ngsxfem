#define FILE_GHOSTPENALTY_CPP
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
      //dtshape = scafe.GetDtShape(mip.IP(), lh);
      dtshape = scafe.GetShape(mip.IP(), lh); // needs to be changed

      mat = 0.0;
      mat.Row(0) = dtshape;


    }


}


