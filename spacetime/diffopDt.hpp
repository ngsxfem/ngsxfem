#ifndef FILE_DIFFOPDT_HPP
#define FILE_DIFFOPDT_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

namespace ngfem
{

template<int SpaceD, int DerivOrder>
  class DiffOpDt : public DiffOp<DiffOpDt<SpaceD, DerivOrder>>
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = SpaceD };    // D-dim space
    enum { DIM_ELEMENT = SpaceD };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    static string Name() { return DerivOrder == 0? "id" : (DerivOrder == 1 ? "dt" : "ddt"); }
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

  template <int SpaceD, int D, int DerivOrder>
  class DiffOpDtVec : public DiffOp<DiffOpDtVec<SpaceD, D, DerivOrder>>
  {

  public:
    enum { DIM = D };          // two copies of the spaces
    enum { DIM_SPACE = SpaceD };    // D-dim space
    enum { DIM_ELEMENT = SpaceD };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = D };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    static string Name() { return DerivOrder == 0? "id" : (DerivOrder == 1 ? "dt" : "ddt"); }
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

  template <int SpaceD, int DerivOrder, VorB VB = VOL>
  class DiffOpDtVectorH1 : public DiffOp<DiffOpDtVectorH1<SpaceD, DerivOrder, VB> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = SpaceD };
    enum { DIM_ELEMENT = SpaceD-VB };
    enum { DIM_DMAT = SpaceD };
    enum { DIFFORDER = 0 };

    static string Name() { return DerivOrder == 0? "id" : (DerivOrder == 1 ? "dt" : "ddt"); }
    static constexpr bool SUPPORT_PML = false;
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh);

  };

  template <int SpaceD, int DerivOrder, VorB VB = VOL>
  class DiffOpDtGradVectorH1 : public DiffOp<DiffOpDtGradVectorH1<SpaceD, DerivOrder, VB> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = SpaceD };
    enum { DIM_ELEMENT = SpaceD-VB };
    enum { DIM_DMAT = SpaceD*SpaceD };
    enum { DIFFORDER = 1 };
    
    static Array<int> GetDimensions() { return Array<int> ( { SpaceD, SpaceD } ); }

    static string Name() { return DerivOrder == 0? "grad" : (DerivOrder == 1 ? "dxt" : "dxtt"); }
    static constexpr bool SUPPORT_PML = false;
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh);

  };

  template <int SpaceD, int DerivOrder, VorB VB = VOL>
  class DiffOpDtDivVectorH1 : public DiffOp<DiffOpDtDivVectorH1<SpaceD, DerivOrder, VB> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = SpaceD };
    enum { DIM_ELEMENT = SpaceD-VB };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };
    
    static string Name() { return DerivOrder == 0? "div" : (DerivOrder == 1 ? "dtdiv" : "dttdiv"); }
    static constexpr bool SUPPORT_PML = false;
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh);

  };


#ifndef FILE_DIFFOPDT_CPP
  extern template class T_DifferentialOperator<DiffOpDt<1,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<2,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<3,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<1,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<2,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<3,1>>;
  extern template class T_DifferentialOperator<DiffOpDt<1,2>>;
  extern template class T_DifferentialOperator<DiffOpDt<2,2>>;
  extern template class T_DifferentialOperator<DiffOpDt<3,2>>;
  extern template class T_DifferentialOperator<DiffOpDt<1,2>>;
  extern template class T_DifferentialOperator<DiffOpDt<2,2>>;
  extern template class T_DifferentialOperator<DiffOpDt<3,2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 1, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 2, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 3, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<1, 1, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<1, 2, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 1, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 2, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 3, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 1, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 2, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 3, 1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 1, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 2, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<0, 3, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<1, 1, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<1, 2, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 1, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 2, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2, 3, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 1, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 2, 2>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<3, 3, 2>>;

#endif

}
#endif
