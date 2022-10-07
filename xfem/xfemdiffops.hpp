#ifndef FILE_XFEMDIFFOPS_HPP
#define FILE_XFEMDIFFOPS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "xfiniteelement.hpp"

namespace ngfem
{

  enum DIFFOPX { EXTEND = 0, RNEG = 1, RPOS = 2, EXTEND_GRAD = 3, RNEG_GRAD = 4, RPOS_GRAD = 5};

  template <int D, DIFFOPX DOX>
  class DiffOpX : public DiffOp<DiffOpX<D,DOX> >
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // D-dim space
    enum { DIM_ELEMENT = D };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = DOX < 3 ? 1 : D };     // D-matrix
    enum { DIFFORDER = DOX < 3 ? 0 : 1 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

#ifndef FILE_XFEMDIFFOPS_CPP
  extern template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::EXTEND>>;
  extern template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::RNEG>>;
  extern template class T_DifferentialOperator<DiffOpX<1,DIFFOPX::RPOS>>;

  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND>>;
  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG>>;
  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS>>;
  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND_GRAD>>;
  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG_GRAD>>;
  extern template class T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS_GRAD>>;

  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND>>;
  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG>>;
  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS>>;
  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND_GRAD>>;
  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG_GRAD>>;
  extern template class T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS_GRAD>>;
#endif

}

#endif

