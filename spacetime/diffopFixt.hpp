#ifndef FILE_DIFFOPFIXT_HPP
#define FILE_DIFFOPFIXT_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

namespace ngfem
{

  template <int SpaceD, int DerivOrder, int time, VorB VB = VOL>
  class DiffOpDtFixtVectorH1 : public DiffOp<DiffOpDtFixtVectorH1<SpaceD, DerivOrder, time, VB> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = SpaceD };
    enum { DIM_ELEMENT = SpaceD-VB };
    enum { DIM_DMAT = SpaceD };
    enum { DIFFORDER = 0 };
    
    static string Name() { return DerivOrder == 0? "fixt" : (DerivOrder == 1 ? "fixtdt" : "fixtdtt"); }
    static constexpr bool SUPPORT_PML = false;
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh);

  };


  template <int SpaceD, int time>
  class DiffOpFixt : public DiffOp<DiffOpFixt<SpaceD, time>>
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = SpaceD };    // D-dim space
    enum { DIM_ELEMENT = SpaceD };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };


  template<int SpaceD>
  class DiffOpFixAnyTime : public DifferentialOperator
  {
    double time;

  public:

    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = SpaceD };    // D-dim space
    enum { DIM_ELEMENT = SpaceD };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    DiffOpFixAnyTime(double atime)
      : DifferentialOperator(DIM_DMAT, 1, VorB(int(DIM_SPACE)-int(DIM_ELEMENT)), DIFFORDER)
    {
      SetDimensions(Array<int> ( { DIM_DMAT } ));
      time = atime;
    }
    /*
    virtual int Dim() const { return DIM_DMAT; }
    virtual bool Boundary() const { return int(DIM_SPACE) > int(DIM_ELEMENT); }
    virtual int DiffOrder() const { return DIFFORDER; }
    */
    virtual string Name() const { return "Fix_time"; }

    virtual bool operator== (const DifferentialOperator & diffop2) const
    { return typeid(*this) == typeid(diffop2); }


    virtual void
    CalcMatrix (const FiniteElement & bfel,
        const BaseMappedIntegrationPoint & bmip,
        BareSliceMatrix<double,ColMajor> mat,
        LocalHeap & lh) const;

    virtual void
    ApplyTrans (const FiniteElement & fel,
        const BaseMappedIntegrationPoint & mip,
        FlatVector<double> flux,
        FlatVector<double> x,
        LocalHeap & lh) const;

  };



#ifndef FILE_DIFFOPFIXT_CPP
  extern template class T_DifferentialOperator<DiffOpFixt<1, 0>>;
  extern template class T_DifferentialOperator<DiffOpFixt<1, 1>>;
  extern template class T_DifferentialOperator<DiffOpFixt<2, 0>>;
  extern template class T_DifferentialOperator<DiffOpFixt<2, 1>>;
  extern template class T_DifferentialOperator<DiffOpFixt<3, 0>>;
  extern template class T_DifferentialOperator<DiffOpFixt<3, 1>>;
  extern template class DiffOpFixAnyTime<2>;
  extern template class DiffOpFixAnyTime<3>;

#endif

}
#endif
