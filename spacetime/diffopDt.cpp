#define FILE_DIFFOPDT_CPP
#include "diffopDt.hpp"
#include <diffop_impl.hpp>
#include "SpaceTimeFE.hpp"
#include "../utils/ngsxstd.hpp"


namespace ngfem
{

  template <int SpaceD, int DerivOrder>
    template<typename FEL, typename MIP, typename MAT>
  void DiffOpDt<SpaceD, DerivOrder>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {
      IntegrationPoint ip(mip.IP());
      mat = 0.0;

      const SpaceTimeFE<SpaceD>* scafed = dynamic_cast<const SpaceTimeFE<SpaceD> * > (& bfel);
      FlatVector<> dtshape (scafed->GetNDof(),lh);
      if (DerivOrder == 1)
        scafed->CalcDtShape(ip,dtshape);
      else
        scafed->CalcDDtShape(ip,dtshape);
      mat.Row(0) = dtshape;
    }

// explicit instantiation (needed for clang starting from NGSolve v6.2.2103 - not sure why exactly though)
  template void DiffOpDt<1,1>::GenerateMatrix(const FiniteElement & bfel, const ngfem::MappedIntegrationPoint<1, 1, double> & mip,
                                              ngbla::FlatMatrixFixHeight<1, double, 1> & mat, LocalHeap & lh);
  template void DiffOpDt<1,2>::GenerateMatrix(const FiniteElement & bfel, const ngfem::MappedIntegrationPoint<1, 1, double> & mip,
                                              ngbla::FlatMatrixFixHeight<1, double, 1> & mat, LocalHeap & lh);

  template class T_DifferentialOperator<DiffOpDt<1,1>>;
  template class T_DifferentialOperator<DiffOpDt<2,1>>;
  template class T_DifferentialOperator<DiffOpDt<3,1>>;
  template class T_DifferentialOperator<DiffOpDt<1,2>>;
  template class T_DifferentialOperator<DiffOpDt<2,2>>;
  template class T_DifferentialOperator<DiffOpDt<3,2>>;

  template <int SpaceD, int D, int DerivOrder>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDtVec<SpaceD, D, DerivOrder>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {
      IntegrationPoint ip(mip.IP());
      mat = 0.0;

      const SpaceTimeFE<SpaceD>& scafed = dynamic_cast<const SpaceTimeFE<SpaceD>& > (bfel);
      FlatVector<> dtshape (scafed.GetNDof(),lh);

      if (DerivOrder == 1)  
        scafed.CalcDtShape(ip,dtshape);
      else  
        scafed.CalcDDtShape(ip,dtshape);
      for (int j = 0; j < D; j++)
        for (int k = 0; k < dtshape.Size(); k++)
            mat(j,k*D+j) = dtshape(k);
    }

  template class T_DifferentialOperator<DiffOpDtVec<0, 1, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<0, 2, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<0, 3, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 1, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 2, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 3, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 1, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 2, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 3, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 1, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 2, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 3, 1>>;
  template class T_DifferentialOperator<DiffOpDtVec<0, 1, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<0, 2, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<0, 3, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 1, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 2, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<1, 3, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 1, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 2, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<2, 3, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 1, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 2, 2>>;
  template class T_DifferentialOperator<DiffOpDtVec<3, 3, 2>>;

  template <int SpaceD, int DerivOrder, VorB VB>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDtVectorH1<SpaceD, DerivOrder, VB>::GenerateMatrix (const FEL & bfel, const MIP & mip,
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
    feli.CalcTimeShape(mip.IP(),time_shape,DerivOrder);

    for (int i = 0; i < tndof; i++)
      for (int d = 0; d < SpaceD; d++)
        for (int j = 0; j < sndof; j++)
          mat(d,i*SpaceD*sndof+d*sndof+j) = space_shape(j)*time_shape(i);
  }        

  template <int SpaceD, int DerivOrder, VorB VB>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDtGradVectorH1<SpaceD, DerivOrder, VB>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                                                 MAT && mat, LocalHeap & lh)
  {
    auto & fel = static_cast<const VectorFiniteElement&> (bfel);
    mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;

    auto & feli = static_cast<const SpaceTimeFE<SpaceD>&> (fel[0]);
    const int tndof = feli.GetTimeFE()->GetNDof();
    const int sndof = feli.GetNDof() / tndof;
    HeapReset hr(lh);
    FlatMatrix<> space_shape (sndof,SpaceD,lh); 
    feli.GetSpaceFE()->CalcMappedDShape(mip,space_shape);
    FlatVector<> time_shape (tndof,lh); 
    feli.CalcTimeShape(mip.IP(),time_shape,DerivOrder);

    for (int i = 0; i < tndof; i++)
      for (int d = 0; d < SpaceD; d++)
        for (int j = 0; j < sndof; j++)
          for (int l = 0; l < SpaceD; l++)
            mat(d*SpaceD+l,i*SpaceD*sndof+d*sndof+j) = space_shape(j,l)*time_shape(i);
  }        

  template <int SpaceD, int DerivOrder, VorB VB>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDtDivVectorH1<SpaceD, DerivOrder, VB>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                                                 MAT && mat, LocalHeap & lh)
  {
    auto & fel = static_cast<const VectorFiniteElement&> (bfel);
    mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;

    auto & feli = static_cast<const SpaceTimeFE<SpaceD>&> (fel[0]);
    const int tndof = feli.GetTimeFE()->GetNDof();
    const int sndof = feli.GetNDof() / tndof;
    HeapReset hr(lh);
    FlatMatrix<> space_shape (sndof,SpaceD,lh); 
    feli.GetSpaceFE()->CalcMappedDShape(mip,space_shape);
    FlatVector<> time_shape (tndof,lh); 
    feli.CalcTimeShape(mip.IP(),time_shape,DerivOrder);

    for (int i = 0; i < tndof; i++)
      for (int d = 0; d < SpaceD; d++)
        for (int j = 0; j < sndof; j++)
          mat(0,i*SpaceD*sndof+d*sndof+j) = space_shape(j,d)*time_shape(i);
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

template class T_DifferentialOperator<DiffOpDtVectorH1<1, 0, BND>>;
template class T_DifferentialOperator<DiffOpDtVectorH1<2, 0, BND>>;
template class T_DifferentialOperator<DiffOpDtVectorH1<3, 0, BND>>;

}


