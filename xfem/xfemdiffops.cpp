#define FILE_XFEMDIFFOPS_CPP
#include "xfemdiffops.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{
  template <int D, DIFFOPX DOX>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpX<D,DOX>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                       MAT & mat, LocalHeap & lh)
  {
    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&bfel);

    if (!xfe)
    {
      mat = 0.0;
    }
    else
    {
      const ScalarFiniteElement<D> & scafe =
        dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

      const int ndof = scafe.GetNDof();

      if (DOX < DIFFOPX::EXTEND_GRAD)
      {
        FlatVector<> shape (ndof,lh);
        shape = scafe.GetShape(mip.IP(), lh);

        if (DOX==DIFFOPX::RNEG || DOX==DIFFOPX::RPOS)
        {
          DOMAIN_TYPE dt_here = DOX==DIFFOPX::RNEG ? DOMAIN_TYPE::NEG : DOMAIN_TYPE::POS;
          const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
          for (int i =0; i < ndof; i++)
            if (xsign[i]==dt_here)
              mat(0,i) = shape(i);
            else
              mat(0,i) = 0.0;
        }
        else
        {
          mat.Row(0) = shape;
        }
      }
      else
      {
        FlatMatrixFixWidth<D> dshape (ndof,lh);
        scafe.CalcMappedDShape(mip, dshape);

        if (DOX==DIFFOPX::RNEG_GRAD || DOX==DIFFOPX::RPOS_GRAD)
        {
          DOMAIN_TYPE dt_here = DOX==DIFFOPX::RNEG_GRAD ? DOMAIN_TYPE::NEG : DOMAIN_TYPE::POS;
          const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
          for (int i =0; i < ndof; i++)
            if (xsign[i]==dt_here)
              mat.Col(i) = dshape.Row(i);
            else
              mat.Col(i) = 0.0;
        }
        else
        {
          mat = Trans(dshape);
        }
      } // GRAD
    }// xfe
  }// generate matrix

  template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::EXTEND>>;
  template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::RNEG>>;
  template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::RPOS>>;

  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND>>;
  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG>>;
  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS>>;
  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND_GRAD>>;
  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG_GRAD>>;
  template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS_GRAD>>;

  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND>>;
  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG>>;
  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS>>;
  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND_GRAD>>;
  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG_GRAD>>;
  template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS_GRAD>>;

}

