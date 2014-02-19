#include "xfemNitsche.hpp"

namespace ngfem
{

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void XNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarFiniteElement<D> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel[i]);
    }

    elmat = 0.0;
    
    if (D==3)
      throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

    if (!xfe) 
    {
      if(dummfe)
        return;
      else
        throw Exception(" not containing X-elements?");
    }

    int ndof_x = xfe->GetNDof();
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> jump(ndof,lh);
    FlatVector<> shape(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);
    const Array<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    int p = scafe->Order();

    const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
    if (lset_eval_p == NULL)
      throw Exception(" no local geometry");
    const CompositeQuadratureRule<D> * compr (lset_eval_p->GetCompositeRule<D>());
    const QuadratureRuleCoDim1<D> & quad(compr->GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);
    const double b_t_neg = beta_neg->Evaluate(mipc);
    const double b_t_pos = beta_pos->Evaluate(mipc);
    const double a_t_neg = alpha_neg->Evaluate(mipc);
    const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;
    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HANSBOBETA:
      {
        throw Exception("No Kappa yet");
        // Vec<2> kappa = masterel.CalcKappa();
        // double sum = kappa(0) * b_t_neg + kappa(1) * b_t_pos;
        // kappa_neg = kappa(0) * b_t_neg / sum;
        // kappa_pos = kappa(1) * b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::BETA:
      {
        double sum = b_t_neg + b_t_pos;
        kappa_neg = b_t_neg / sum;
        kappa_pos = b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::ALPHA:
      {
        double sum = a_t_neg + a_t_pos;
        kappa_neg = a_t_neg / sum;
        kappa_pos = a_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::ALPHABETA:
      {
        double sum = a_t_neg / b_t_neg + a_t_pos / b_t_pos;
        kappa_neg = a_t_neg / b_t_neg / sum;
        kappa_pos = a_t_pos / b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::HANSBO:
    default:
      {
        throw Exception("No Kappa yet");
        // Vec<2> kappa = masterel.CalcKappa();
        // kappa_neg = kappa(0);
        // kappa_pos = kappa(1);
        break;	      
      }
    }

    for (int i = 0; i < quad.Size(); ++i)
    {
      IntegrationPoint ip(quad.points[i]);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      const double h = D == 2 ? sqrt(absdet) : cbrt(absdet);

      Vec<D> nref = quad.normals[i];
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;

      const double weight = quad.weights[i] * len;
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);
      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
      const double lam = lambda->Evaluate(mip);
        
      shape = scafe->GetShape(mip.IP(), lh);
      jump.Range(0,ndof_h1) = (b_pos-b_neg) * shape;
      jump.Range(ndof_h1,ndof) = shape;

      scafe->CalcMappedDShape (mip,dshape_h1);

      dshape.Range(0,ndof_h1) = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape_h1 * normal);
      dshape.Range(ndof_h1,ndof) = dshape_h1 * normal;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l] == POS){
          jump(ndof_h1+l) *= -b_neg;
          dshape(ndof_h1+l) *= kappa_neg * a_neg;
        }
        else{
          jump(ndof_h1+l) *= b_pos;
          dshape(ndof_h1+l) *= kappa_pos * a_pos;
        }
      }

      double ava = a_pos;

      switch (kappa_choice){
      case NITSCHE_VARIANTS::HALFHALF:
        {
        ava = a_pos*0.5+a_neg*0.5;
        break;
        }
      case NITSCHE_VARIANTS::BETA:
      case NITSCHE_VARIANTS::ALPHA:
        {
          ava = 2*a_pos*a_neg/(a_neg+a_pos);
        }
      case NITSCHE_VARIANTS::ALPHABETA:
        {
          ava = a_pos*kappa_pos+a_neg*kappa_neg;
        }
      case NITSCHE_VARIANTS::HANSBO:
      default:
        {
          ava = a_pos*kappa_pos+a_neg*kappa_neg;
        }
      }


      elmat -= weight * jump * Trans(dshape);
      elmat -= weight * dshape * Trans(jump);
      elmat += lam * (p+1)*p/h * ava * weight * jump * Trans(jump);

    }

  }

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void SpaceTimeXNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarSpaceTimeFiniteElement<D> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel[i]);
    }

    elmat = 0.0;
    
    if (D==3)
      throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

    if (!xfe) 
    {
      if(dummfe)
        return;
      else
        throw Exception(" not containing X-elements?");
    }

    const double t0 = t_old->EvaluateConst();
    const double t1 = t_new->EvaluateConst();
    const double dt = t1-t0;

    int ndof_x = xfe->GetNDof();
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> jump(ndof,lh);
    FlatVector<> shape(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);
    const Array<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    int ps = scafe->OrderSpace();
    // int pt = scafe->OrderTime();

    const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
    if (lset_eval_p == NULL)
      throw Exception(" no local geometry");
    const CompositeQuadratureRule<D+1> * compr (lset_eval_p->GetCompositeRule<D+1>());
    const QuadratureRuleCoDim1<D+1> & quad(compr->GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);
    const double b_t_neg = beta_neg->Evaluate(mipc);
    const double b_t_pos = beta_pos->Evaluate(mipc);
    const double a_t_neg = alpha_neg->Evaluate(mipc);
    const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;
    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HANSBOBETA:
      {
        throw Exception("No Kappa yet");
        // Vec<2> kappa = masterel.CalcKappa();
        // double sum = kappa(0) * b_t_neg + kappa(1) * b_t_pos;
        // kappa_neg = kappa(0) * b_t_neg / sum;
        // kappa_pos = kappa(1) * b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::BETA:
      {
        double sum = b_t_neg + b_t_pos;
        kappa_neg = b_t_neg / sum;
        kappa_pos = b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::ALPHA:
      {
        double sum = a_t_neg + a_t_pos;
        kappa_neg = a_t_neg / sum;
        kappa_pos = a_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::ALPHABETA:
      {
        double sum = a_t_neg / b_t_neg + a_t_pos / b_t_pos;
        kappa_neg = a_t_neg / b_t_neg / sum;
        kappa_pos = a_t_pos / b_t_pos / sum;
        break;
      }
    case NITSCHE_VARIANTS::HANSBO:
    default:
      {
        throw Exception("No Kappa yet");
        // Vec<2> kappa = masterel.CalcKappa();
        // kappa_neg = kappa(0);
        // kappa_pos = kappa(1);
        break;	      
      }
    }

    for (int i = 0; i < quad.Size(); ++i)
    {
      IntegrationPoint ip(quad.points[i]);
      const double time = quad.points[i][D];
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      const double h = D == 2 ? sqrt(absdet) : cbrt(absdet);

      Vec<D> nref_space; 
      for (int d = 0; d < D; ++d) 
        nref_space[d] = quad.normals[i][d];
      Vec<D> normal_space = dt * absdet * Trans(Finv) * nref_space;
      double n_t = quad.normals[i][D] * absdet;

      Vec<D+1> normal_st; 
      for (int d = 0; d < D; ++d) 
        normal_st[d] = normal_space[d];
      normal_st[D] = n_t;
      
      double len = L2Norm(normal_st);
      normal_st /= len;

      for (int d = 0; d < D; ++d) 
        normal_space[d] = normal_st[d];

      const double nu = L2Norm(normal_space);
      normal_space /= nu;

      const double weight = quad.weights[i] * len * nu;
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);
      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
      const double lam = lambda->Evaluate(mip);
        
      scafe->CalcShapeSpaceTime(mip.IP(), time, shape, lh);
      jump.Range(0,ndof_h1) = (b_pos-b_neg) * shape;
      jump.Range(ndof_h1,ndof) = shape;

      scafe->CalcMappedDxShapeSpaceTime(mip, time, dshape_h1, lh);

      dshape.Range(0,ndof_h1) = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape_h1 * normal_space);
      dshape.Range(ndof_h1,ndof) = dshape_h1 * normal_space;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l] == POS){
          jump(ndof_h1+l) *= -b_neg;
          dshape(ndof_h1+l) *= kappa_neg * a_neg;
        }
        else{
          jump(ndof_h1+l) *= b_pos;
          dshape(ndof_h1+l) *= kappa_pos * a_pos;
        }
      }

      double ava = a_pos;

      switch (kappa_choice){
      case NITSCHE_VARIANTS::HALFHALF:
        {
        ava = a_pos*0.5+a_neg*0.5;
        break;
        }
      case NITSCHE_VARIANTS::BETA:
      case NITSCHE_VARIANTS::ALPHA:
        {
          ava = 2*a_pos*a_neg/(a_neg+a_pos);
        }
      case NITSCHE_VARIANTS::ALPHABETA:
        {
          ava = a_pos*kappa_pos+a_neg*kappa_neg;
        }
      case NITSCHE_VARIANTS::HANSBO:
      default:
        {
          ava = a_pos*kappa_pos+a_neg*kappa_neg;
        }
      }


      elmat -= weight * jump * Trans(dshape);
      elmat -= weight * dshape * Trans(jump);
      elmat += lam * (ps+1)*ps/h * ava * weight * jump * Trans(jump);

    }

  }

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitsche2d_1 ("xnitsche_halfhalf", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_2 ("xnitsche_hansbo", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_2b ("xnitsche", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBOBETA> > initxnitsche2d_3 ("xnitsche_hansbobeta", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::BETA> > initxnitsche2d_4 ("xnitsche_beta", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHA> > initxnitsche2d_5 ("xnitsche_alpha", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHABETA> > initxnitsche2d_6 ("xnitsche_alphabeta", 2, 5);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initxnitsche3d_1 ("xnitsche_halfhalf", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_2 ("xnitsche_hansbo", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_2b ("xnitsche", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBOBETA> > initxnitsche3d_3 ("xnitsche_hansbobeta", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::BETA> > initxnitsche3d_4 ("xnitsche_beta", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHA> > initxnitsche3d_5 ("xnitsche_alpha", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHABETA> > initxnitsche3d_6 ("xnitsche_alphabeta", 3, 5);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitsche2d_st_1 ("stx_nitsche_halfhalf", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_st_2 ("stx_nitsche_hansbo", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_st_2b ("stx_nitsche", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBOBETA> > initxnitsche2d_st_3 ("stx_nitsche_hansbobeta", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::BETA> > initxnitsche2d_st_4 ("stx_nitsche_beta", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHA> > initxnitsche2d_st_5 ("stx_nitsche_alpha", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHABETA> > initxnitsche2d_st_6 ("stx_nitsche_alphabeta", 2, 7);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initxnitsche3d_st_1 ("stx_nitsche_halfhalf", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_st_2 ("stx_nitsche_hansbo", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_st_2b ("stx_nitsche", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBOBETA> > initxnitsche3d_st_3 ("stx_nitsche_hansbobeta", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::BETA> > initxnitsche3d_st_4 ("stx_nitsche_beta", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHA> > initxnitsche3d_st_5 ("stx_nitsche_alpha", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHABETA> > initxnitsche3d_st_6 ("stx_nitsche_alphabeta", 3, 7);

}

/// coefficientfunction statt function-pointer
