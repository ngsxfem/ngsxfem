#ifndef FILE_DIFFOPDT_HPP
#define FILE_DIFFOPDT_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

namespace ngfem
{

  class DiffOpDt : public DiffOp<DiffOpDt>
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = 2 };    // D-dim space
    enum { DIM_ELEMENT = 2 };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

  template <int D>
  class DiffOpDtVec : public DiffOp<DiffOpDtVec<D>>
  {

  public:
    enum { DIM = D };          // two copies of the spaces
    enum { DIM_SPACE = 2 };    // D-dim space
    enum { DIM_ELEMENT = 2 };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = D };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };


  template <int time>
  class DiffOpFixt : public DiffOp<DiffOpFixt<time>>
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = 2 };    // D-dim space
    enum { DIM_ELEMENT = 2 };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };



  class DiffOpFixAnyTime : public DifferentialOperator
  {
    double time;

  public:

    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = 2 };    // D-dim space
    enum { DIM_ELEMENT = 2 };  // D-dim elements (in contrast to boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    DiffOpFixAnyTime(double atime)
      : DifferentialOperator(DIM_DMAT, 1, VorB(int(DIM_SPACE)-int(DIM_ELEMENT)), DIFFORDER)
    {
      //dimensions = DIFFOP::GetDimensions();
      dimensions = Array<int> ( { DIM_DMAT } );
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
        SliceMatrix<double,ColMajor> mat,
        LocalHeap & lh) const;

    virtual void
    ApplyTrans (const FiniteElement & fel,
        const BaseMappedIntegrationPoint & mip,
        FlatVector<double> flux,
        FlatVector<double> x,
        LocalHeap & lh) const;

  };



#ifndef FILE_DIFFOPDT_CPP
  extern template class T_DifferentialOperator<DiffOpDt>;
  extern template class T_DifferentialOperator<DiffOpDtVec<1>>;
  extern template class T_DifferentialOperator<DiffOpDtVec<2>>;
  extern template class T_DifferentialOperator<DiffOpFixt<0>>;
  extern template class T_DifferentialOperator<DiffOpFixt<1>>;
#endif

}
#endif
