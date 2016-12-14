#define FILE_XFEMVISINTS_CPP
#include "xfemVisInts.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{



  
/*  
  template <int D, TIME t>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalSTX<D,t>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarSpaceTimeFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> & > (cfel[0]);

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);

    scafe.CalcShapeSpaceTime(mip.IP(),t == PAST ? 0.0 : 1.0,shape,lh);

    IntRange range = cfel.GetRange(0);
    mat.Row(0).Range(range) = shape;
    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&cfel[1]);

    if (xfe)
    {
      const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
      const double lsetval = xgeom.EvaluateLsetAtPoint<D,D+1>(mip.IP(),t == PAST ? 0.0 : 1.0);
      
      DOMAIN_TYPE dt_here = lsetval > 0 ? POS : NEG;
      const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
      for (int i =0; i < ndof; i++)
        if (xsign[i]==dt_here)
          mat(0,ndof+i) = shape(i);
        else
          mat(0,ndof+i) = 0.0;
    } 
  }


  template <int D, TIME t>  STXVisIntegrator<D,t> :: STXVisIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalSTX<D,t>, DiagDMat<1>, CompoundFiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D, TIME t>  STXVisIntegrator<D,t> :: STXVisIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSTX<D,t>, DiagDMat<1>, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D, TIME t>  STXVisIntegrator<D,t> :: ~STXVisIntegrator () { ; }

  template class STXVisIntegrator<2,PAST>;
  template class STXVisIntegrator<3,PAST>;
  template class STXVisIntegrator<2,FUTURE>;
  template class STXVisIntegrator<3,FUTURE>;

  static RegisterBilinearFormIntegrator<STXVisIntegrator<2,PAST> > initxwmass0p ("stxvis_past", 2, 1);
  static RegisterBilinearFormIntegrator<STXVisIntegrator<3,PAST> > initxwmass1p ("stxvis_past", 3, 1);

  static RegisterBilinearFormIntegrator<STXVisIntegrator<2,FUTURE> > initxwmass0f ("stxvis_future", 2, 1);
  static RegisterBilinearFormIntegrator<STXVisIntegrator<3,FUTURE> > initxwmass1f ("stxvis_future", 3, 1);


  template <int D, TIME t>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalSTNegPos<D,t>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarSpaceTimeFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> & > (cfel[0]);

    const LevelsetContainerFE & lsetcontfe = 
      dynamic_cast<const LevelsetContainerFE & > (cfel[2]);

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);

    scafe.CalcShapeSpaceTime(mip.IP(),t == PAST ? 0.0 : 1.0,shape,lh);

    IntRange range0 = cfel.GetRange(0);
    IntRange range1 = cfel.GetRange(1);

    const shared_ptr<CoefficientFunction> coef_lset = lsetcontfe.GetLevelsetCoefficient();

    DimMappedIntegrationPoint<D+1> mipp(mip.IP(),mip.GetTransformation());
    mipp.Point().Range(0,D) = mip.GetPoint();
    mipp.Point()[D] = t == PAST ? lsetcontfe.told : lsetcontfe.tnew;

    const double lsetval = coef_lset->Evaluate(mipp);

    if (lsetval < 0)
    {
      mat.Row(0).Range(range0) = shape;
      mat.Row(0).Range(range1) = 0.0;
    }
    else
    {
      mat.Row(0).Range(range0) = 0.0;
      mat.Row(0).Range(range1) = shape;
    }
  }


  template <int D, TIME t>  STNegPosVisIntegrator<D,t> :: STNegPosVisIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalSTNegPos<D,t>, DiagDMat<1>, CompoundFiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D, TIME t>  STNegPosVisIntegrator<D,t> :: STNegPosVisIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSTNegPos<D,t>, DiagDMat<1>, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D, TIME t>  STNegPosVisIntegrator<D,t> :: ~STNegPosVisIntegrator () { ; }

  template class STNegPosVisIntegrator<2,PAST>;
  template class STNegPosVisIntegrator<3,PAST>;
  template class STNegPosVisIntegrator<2,FUTURE>;
  template class STNegPosVisIntegrator<3,FUTURE>;

  static RegisterBilinearFormIntegrator<STNegPosVisIntegrator<2,PAST> > initnegposxwmass0p ("st_np_vis_past", 2, 1);
  static RegisterBilinearFormIntegrator<STNegPosVisIntegrator<3,PAST> > initnegposxwmass1p ("st_np_vis_past", 3, 1);

  static RegisterBilinearFormIntegrator<STNegPosVisIntegrator<2,FUTURE> > initnegposxwmass0f ("st_np_vis_future", 2, 1);
  static RegisterBilinearFormIntegrator<STNegPosVisIntegrator<3,FUTURE> > initnegposxwmass1f ("st_np_vis_future", 3, 1);

*/

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
  
  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalX<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarFiniteElement<D> & > (cfel[0]);

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);

    shape = scafe.GetShape(mip.IP(), lh);
    IntRange range = cfel.GetRange(0);
    mat.Row(0).Range(range) = shape;
    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&cfel[1]);

    if (xfe)
    {
      const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
      const double lsetval = xgeom.EvaluateLsetAtPoint<D,D>(mip.IP(),0.0);

      DOMAIN_TYPE dt_here = lsetval > 0 ? POS : NEG;
      const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
      for (int i =0; i < ndof; i++)
        if (xsign[i]==dt_here)
          mat(0,ndof+i) = shape(i);
        else
          mat(0,ndof+i) = 0.0;
    } 
  }

  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpGradX<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                       MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarFiniteElement<D> & > (cfel[0]);

    const int ndof = scafe.GetNDof();
    FlatMatrix<> dshape (ndof,D,lh);

    scafe.CalcMappedDShape(mip, dshape);
    IntRange range = cfel.GetRange(0);
    mat.Cols(range) = Trans(dshape);
    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&cfel[1]);

    if (xfe)
    {
      const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
      const double lsetval = xgeom.EvaluateLsetAtPoint<D,D>(mip.IP(),0.0);

      DOMAIN_TYPE dt_here = lsetval > 0 ? POS : NEG;
      const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
      for (int i =0; i < ndof; i++)
        if (xsign[i]==dt_here)
          mat.Col(ndof+i) = dshape.Row(i);
        else
          mat.Col(ndof+i) = Vec<D>(0.0);
    } 
    
  }
  
  template <int D>  XVisIntegrator<D> :: XVisIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalX<D>, DiagDMat<1>, CompoundFiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D>  XVisIntegrator<D> :: XVisIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalX<D>, DiagDMat<1>, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D>  XVisIntegrator<D> :: ~XVisIntegrator () { ; }

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

  
  template class T_DifferentialOperator<DiffOpGradX<2>>;
  template class T_DifferentialOperator<DiffOpGradX<3>>;

  template class T_DifferentialOperator<DiffOpEvalX<2>>;
  template class T_DifferentialOperator<DiffOpEvalX<3>>;

  template class XVisIntegrator<2>;
  template class XVisIntegrator<3>;

  static RegisterBilinearFormIntegrator<XVisIntegrator<2> > initxwmass0 ("xvis", 2, 1);
  static RegisterBilinearFormIntegrator<XVisIntegrator<3> > initxwmass1 ("xvis", 3, 1);



  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalSigned<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarFiniteElement<D> & > (cfel[0]);

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);
    shape = scafe.GetShape(mip.IP(), lh);
    IntRange range = cfel.GetRange(0);
    mat.Row(0).Range(range) = shape;
    mat.Row(1).Range(range) = shape;

    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&cfel[1]);

    if (xfe)
    {
      const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
      for (int i =0; i < xfe->GetNDof(); i++)
        if (xsign[i]==POS){
          mat(0,ndof+i) = shape(i);
          mat(1,ndof+i) = 0.0;
        }
        else
        {
          mat(0,ndof+i) = 0.0;
          mat(1,ndof+i) = shape(i);
        }
    } 
  }


  /*
  template <int D, TIME t>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalSpaceTimeSigned<D,t>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                              MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    const ScalarSpaceTimeFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> & > (cfel[0]);

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);

    if (t == PAST)
      scafe.CalcShapeSpaceTime(mip.IP(),0.0,shape,lh);
    else if (t==FUTURE)
      scafe.CalcShapeSpaceTime(mip.IP(),1.0,shape,lh);

    IntRange range = cfel.GetRange(0);
    mat.Row(0).Range(range) = shape;
    mat.Row(1).Range(range) = shape;

    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&cfel[1]);

    if (xfe)
    {
      const FlatArray<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
      for (int i =0; i < xfe->GetNDof(); i++)
        if (xsign[i]==POS){
          mat(0,ndof+i) = shape(i);
          mat(1,ndof+i) = 0.0;
        }
        else
        {
          mat(0,ndof+i) = 0.0;
          mat(1,ndof+i) = shape(i);
        }
    } 
  }

  */

  template <typename FEL, typename MIP, typename MAT>
  void XHeavisideDMat::GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0.0;
    double val = lvlset -> Evaluate (mip);
    if (val>=0.0)
      mat(0, 0) = 1.0;
    else
      mat(1, 1) = 1.0;
  }  

  template <int D>  SignedXMassIntegrator<D> :: 
  SignedXMassIntegrator (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalSigned<D>, XHeavisideDMat, CompoundFiniteElement > (XHeavisideDMat (coeff))
  { ; }

  template <int D>  SignedXMassIntegrator<D> :: 
  SignedXMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSigned<D>, XHeavisideDMat, CompoundFiniteElement > (coeffs)
  { ; }

  template class SignedXMassIntegrator<2>;
  template class SignedXMassIntegrator<3>;

  static RegisterBilinearFormIntegrator<SignedXMassIntegrator<2> > initxvmass0 ("xvis_sign", 2, 1);
  static RegisterBilinearFormIntegrator<SignedXMassIntegrator<3> > initxvmass1 ("xvis_sign", 3, 1);

  template <int D>  SignedXMassIntegrator<D> :: ~SignedXMassIntegrator () { ; }
/*
  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: 
  SignedSpaceTimeXMassIntegrator (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalSpaceTimeSigned<D,t>, XHeavisideDMat, CompoundFiniteElement > (XHeavisideDMat (coeff))
  { ; }

  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: 
  SignedSpaceTimeXMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSpaceTimeSigned<D,t>, XHeavisideDMat, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: ~SignedSpaceTimeXMassIntegrator () { ; }

  template class SignedSpaceTimeXMassIntegrator<2,PAST>;
  template class SignedSpaceTimeXMassIntegrator<3,PAST>;
  template class SignedSpaceTimeXMassIntegrator<2,FUTURE>;
  template class SignedSpaceTimeXMassIntegrator<3,FUTURE>;

  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<2,PAST> > initxv_stp_mass0 ("xvis_st_past", 2, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<3,PAST> > initxv_stp_mass1 ("xvis_st_past", 3, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<2,FUTURE> > initxv_stf_mass0 ("xvis_st_future", 2, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<3,FUTURE> > initxv_stf_mass1 ("xvis_st_future", 3, 1);
*/



  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalFict<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const XFiniteElement * xfe = 
      dynamic_cast<const XFiniteElement *> (&bfel);

    if (!xfe)
    {
      mat = 0.0;
      return;
    }

    const ScalarFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    const int ndof = scafe.GetNDof();
    FlatVector<> shape (ndof,lh);

    shape = scafe.GetShape(mip.IP(), lh);
    mat.Row(0) = shape;
  }


  template <int D>  FictVisIntegrator<D> :: FictVisIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalFict<D>, DiagDMat<1>, FiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D>  FictVisIntegrator<D> :: FictVisIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalFict<D>, DiagDMat<1>, FiniteElement > (coeffs)
  { ; }

  template <int D>  FictVisIntegrator<D> :: ~FictVisIntegrator () { ; }

  template class FictVisIntegrator<2>;
  template class FictVisIntegrator<3>;

  static RegisterBilinearFormIntegrator<FictVisIntegrator<2> > initfictwmass0 ("fictvis", 2, 1);
  static RegisterBilinearFormIntegrator<FictVisIntegrator<3> > initfictwmass1 ("fictvis", 3, 1);




// ---------------------

  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalExtTrace<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const XFiniteElement * xfe = dynamic_cast<const XFiniteElement *> (&bfel);
    if (!xfe)
      return;
    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());
    FlatVector<> shape (scafe.GetNDof(), &mat(0,0));
    shape = scafe.GetShape(mip.IP(), lh);
  }


  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpGradExtTrace<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const XFiniteElement * xfe = dynamic_cast<const XFiniteElement *> (&bfel);
    if (!xfe)
      return;
    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());
    FlatMatrix<> dshape (scafe.GetNDof(), D, &mat(0,0));
    scafe.CalcMappedDShape(mip, dshape);
  }

  
  template <int D>  ExtTraceIntegrator<D> :: ExtTraceIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalExtTrace<D>, DiagDMat<1>, CompoundFiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D>  ExtTraceIntegrator<D> :: ExtTraceIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalExtTrace<D>, DiagDMat<1>, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D>  ExtTraceIntegrator<D> :: ~ExtTraceIntegrator () { ; }

  template class T_DifferentialOperator<DiffOpEvalExtTrace<2>>;
  template class T_DifferentialOperator<DiffOpEvalExtTrace<3>>;
  template class DiffOpEvalExtTrace<2>;
  template class DiffOpEvalExtTrace<3>;

  template class T_DifferentialOperator<DiffOpGradExtTrace<2>>;
  template class T_DifferentialOperator<DiffOpGradExtTrace<3>>;
  template class DiffOpGradExtTrace<2>;
  template class DiffOpGradExtTrace<3>;
  
  static RegisterBilinearFormIntegrator<ExtTraceIntegrator<2> > initexttracemass0 ("exttrace", 2, 1);
  static RegisterBilinearFormIntegrator<ExtTraceIntegrator<3> > initexttracemass1 ("exttrace", 3, 1);

// ---------------------



  


}

