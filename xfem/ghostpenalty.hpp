#ifndef FILE_GHOSTPENALTY_HPP
#define FILE_GHOSTPENALTY_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

namespace ngfem
{

  template <int D, int ORDER>
  class DiffOpDuDnk : public DiffOp<DiffOpDuDnk<D,ORDER> >
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // D-dim space
    enum { DIM_ELEMENT = D };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = ORDER };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };


  template <int D, int ORDER>
  class DiffOpDuDnkHDiv : public DiffOp<DiffOpDuDnkHDiv<D,ORDER> >
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // D-dim space
    enum { DIM_ELEMENT = D };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = D };     // D-matrix
    enum { DIFFORDER = ORDER };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };


#ifndef FILE_GHOSTPENALTY_CPP
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,1>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,2>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,3>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,4>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,5>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,6>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,7>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<2,8>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,1>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,2>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,3>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,4>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,5>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,6>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,7>>;
  extern template class T_DifferentialOperator<DiffOpDuDnk<3,8>>;

  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,1>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,2>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,3>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,4>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,5>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,6>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,7>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,8>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,1>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,2>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,3>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,4>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,5>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,6>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,7>>;
  extern template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,8>>;
#endif

}

#endif
