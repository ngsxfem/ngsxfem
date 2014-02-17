#include "xfemVisInts.hpp"

namespace ngfem
{

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
      const Array<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
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
      const Array<DOMAIN_TYPE> & xsign = xfe->GetSignsOfDof();
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
  SignedXMassIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpEvalSigned<D>, XHeavisideDMat, CompoundFiniteElement > (XHeavisideDMat (coeff))
  { ; }

  template <int D>  SignedXMassIntegrator<D> :: 
  SignedXMassIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSigned<D>, XHeavisideDMat, CompoundFiniteElement > (coeffs)
  { ; }

  template class SignedXMassIntegrator<2>;
  template class SignedXMassIntegrator<3>;

  static RegisterBilinearFormIntegrator<SignedXMassIntegrator<2> > initxvmass0 ("visx", 2, 1);
  static RegisterBilinearFormIntegrator<SignedXMassIntegrator<3> > initxvmass1 ("visx", 3, 1);

  template <int D>  SignedXMassIntegrator<D> :: ~SignedXMassIntegrator () { ; }

  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: 
  SignedSpaceTimeXMassIntegrator (CoefficientFunction * coeff)
    : T_BDBIntegrator<DiffOpEvalSpaceTimeSigned<D,t>, XHeavisideDMat, CompoundFiniteElement > (XHeavisideDMat (coeff))
  { ; }

  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: 
  SignedSpaceTimeXMassIntegrator (Array<CoefficientFunction*> & coeffs)
    : T_BDBIntegrator<DiffOpEvalSpaceTimeSigned<D,t>, XHeavisideDMat, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D, TIME t>  SignedSpaceTimeXMassIntegrator<D,t> :: ~SignedSpaceTimeXMassIntegrator () { ; }

  template class SignedSpaceTimeXMassIntegrator<2,PAST>;
  template class SignedSpaceTimeXMassIntegrator<3,PAST>;
  template class SignedSpaceTimeXMassIntegrator<2,FUTURE>;
  template class SignedSpaceTimeXMassIntegrator<3,FUTURE>;

  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<2,PAST> > initxv_stp_mass0 ("visx_st_past", 2, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<3,PAST> > initxv_stp_mass1 ("visx_st_past", 3, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<2,FUTURE> > initxv_stf_mass0 ("visx_st_future", 2, 1);
  static RegisterBilinearFormIntegrator<SignedSpaceTimeXMassIntegrator<3,FUTURE> > initxv_stf_mass1 ("visx_st_future", 3, 1);


}

