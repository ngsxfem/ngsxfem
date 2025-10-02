#define FILE_DIFFOPFIXT_CPP
#include "diffopFixt.hpp"
#include <diffop_impl.hpp>
#include "SpaceTimeFE.hpp"
#include "../utils/ngsxstd.hpp"


namespace ngfem
{

  template <int SpaceD, int DerivOrder, int time, VorB VB>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDtFixtVectorH1<SpaceD, DerivOrder, time, VB>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                                                           MAT && mat, LocalHeap & lh)
  {
    auto & fel = static_cast<const VectorFiniteElement&> (bfel);
    mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;

    auto & feli = static_cast<const SpaceTimeFE<SpaceD>&> (fel[0]);
    const int tndof = feli.GetTimeFE()->GetNDof();
    const int sndof = feli.GetNDof() / tndof;
    HeapReset hr(lh);
    FlatVector<> space_shape (sndof,lh); 
    feli.GetSpaceFE()->CalcShape(mip.IP(),space_shape);
    FlatVector<> time_shape (tndof,lh); 


    IntegrationPoint ip(mip.IP()(0),mip.IP()(1), mip.IP()(2), time);
    MarkAsSpaceTimeIntegrationPoint(ip);    
    feli.CalcTimeShape(ip,time_shape,DerivOrder);

    for (int i = 0; i < tndof; i++)
      for (int d = 0; d < SpaceD; d++)
        for (int j = 0; j < sndof; j++)
          mat(d,i*SpaceD*sndof+d*sndof+j) = space_shape(j)*time_shape(i);
  }        



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


#define INSTANTIATE_ONE(CLASS, SPACEDIM, TD, VB) \
  template class T_DifferentialOperator<CLASS<SPACEDIM, TD, VB>>;

#define INSTANTIATE_ALL_VB(CLASS, SPACEDIM, TD) \
  INSTANTIATE_ONE(CLASS, SPACEDIM, TD, VOL) 
//  \
//  INSTANTIATE_ONE(CLASS, SPACEDIM, TD, BND)

#define INSTANTIATE_ALL_TIME_DERIVS(CLASS, SPACEDIM) \
  INSTANTIATE_ALL_VB(CLASS, SPACEDIM, 0) \
  INSTANTIATE_ALL_VB(CLASS, SPACEDIM, 1) \
  INSTANTIATE_ALL_VB(CLASS, SPACEDIM, 2)
 
#define INSTANTIATE_ALL_SPACEDIMS(CLASS) \
  INSTANTIATE_ALL_TIME_DERIVS(CLASS, 0) \
  INSTANTIATE_ALL_TIME_DERIVS(CLASS, 1) \
  INSTANTIATE_ALL_TIME_DERIVS(CLASS, 2) \
  INSTANTIATE_ALL_TIME_DERIVS(CLASS, 3)

  INSTANTIATE_ALL_SPACEDIMS(DiffOpDtVectorH1)
  INSTANTIATE_ALL_SPACEDIMS(DiffOpDtGradVectorH1)
  INSTANTIATE_ALL_SPACEDIMS(DiffOpDtDivVectorH1)
  

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

  template class DiffOpFixAnyTime<1>;
  template class DiffOpFixAnyTime<2>;
  template class DiffOpFixAnyTime<3>;

}


