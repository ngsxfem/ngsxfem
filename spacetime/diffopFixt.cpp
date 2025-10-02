#define FILE_DIFFOPFIXT_CPP
#include "diffopFixt.hpp"
#include <diffop_impl.hpp>
#include "SpaceTimeFE.hpp"
#include "../utils/ngsxstd.hpp"


namespace ngfem
{

  template <int SpaceD, int time>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpFixt<SpaceD, time>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {

      IntegrationPoint ip(mip.IP()(0),mip.IP()(1), mip.IP()(2), time);
      MarkAsSpaceTimeIntegrationPoint(ip);
      mat = 0.0;
      const SpaceTimeFE<SpaceD>& scafed = dynamic_cast<const SpaceTimeFE<SpaceD> & > (bfel);

      FlatVector<> shape (scafed.GetNDof(),lh);
      scafed.CalcShape(ip,shape);
      mat.Row(0) = shape;
   }

  template<int SpaceD>
  void DiffOpFixAnyTime<SpaceD> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {

      //mat = 0.0;
      const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

      IntegrationPoint ip(mip.IP()(0),mip.IP()(1),mip.IP()(2), time);
      MarkAsSpaceTimeIntegrationPoint(ip);

      const SpaceTimeFE<SpaceD>& scafed = dynamic_cast<const SpaceTimeFE<SpaceD> &> (bfel);

      FlatVector<> shape (scafed.GetNDof(),lh);
      scafed.CalcShape(ip,shape);
      mat.Row(0) = shape;
  }

  template<int SpaceD>
  void DiffOpFixAnyTime<SpaceD> ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x,
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }


  template<int SpaceD>
  void DiffOpFixAnyTimeVectorH1<SpaceD> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {
    DiffOpDtFixtVectorH1_GenerateMatrix_Impl<SpaceD, 0, FiniteElement, BaseMappedIntegrationPoint, BareSliceMatrix<double,ColMajor>> 
                            (bfel, bmip, mat, DIM_DMAT, time, lh);
  }

  template<int SpaceD>
  void DiffOpFixAnyTimeVectorH1<SpaceD> ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x,
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }


  template class T_DifferentialOperator<DiffOpFixt<1, 0>>;
  template class T_DifferentialOperator<DiffOpFixt<1, 1>>;

  template class T_DifferentialOperator<DiffOpFixt<2, 0>>;
  template class T_DifferentialOperator<DiffOpFixt<2, 1>>;

  template class T_DifferentialOperator<DiffOpFixt<3, 0>>;
  template class T_DifferentialOperator<DiffOpFixt<3, 1>>;

  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 0, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 0, 1, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 1, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 1, 1, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 0, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 0, 1, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 1, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 1, 1, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 0, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 0, 1, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 1, 0, VOL>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 1, 1, VOL>>;

  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 0, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 0, 1, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 1, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<1, 1, 1, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 0, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 0, 1, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 1, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<2, 1, 1, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 0, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 0, 1, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 1, 0, BND>>;
  template class T_DifferentialOperator<DiffOpDtFixtVectorH1<3, 1, 1, BND>>;

  template class DiffOpFixAnyTime<1>;
  template class DiffOpFixAnyTime<2>;
  template class DiffOpFixAnyTime<3>;

  template class DiffOpFixAnyTimeVectorH1<1>;
  template class DiffOpFixAnyTimeVectorH1<2>;
  template class DiffOpFixAnyTimeVectorH1<3>;


  
}


