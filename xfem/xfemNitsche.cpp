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
    static Timer timer ("XNitscheIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

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
    FlatVector<> shape_total(ndof,lh);
    FlatVector<> shape(ndof_h1,&shape_total(0));
    FlatVector<> shape_x(ndof_x,&shape_total(ndof_h1));
    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);

    FlatMatrix<> Nc(ndof,ndof,lh);
    FlatMatrix<> Ns(ndof,ndof,lh);
    FlatMatrix<> A(ndof,ndof,lh);

    FlatMatrix<> Lsys(ndof+2,ndof+2,lh);
    FlatMatrix<> L(ndof,ndof,lh); // lifting matrix

    ablockintegrator-> CalcElementMatrix (base_fel,
                                          eltrans, 
                                          A, lh);
    Ns = 0.0;
    Nc = 0.0;

    FlatMatrixFixWidth<2> constrb(ndof,lh);
    constrb = 0.0;

    const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    int p = scafe->Order();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    { // calc constrb
      DOMAIN_TYPE dt = POS;
      for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
      {
        const FlatQuadratureRule<D> & fdomquad(fcompr.GetRule(dt));
        for (int i = 0; i < fdomquad.Size(); ++i)
        {
          IntegrationPoint ip(&fdomquad.points(i,0),fdomquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);

          shape = scafe->GetShape(ip, lh);
          shape_x = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shape_x(l) = 0.0;
          }

          constrb.Col(dt) += mip.GetWeight() * shape_total;
        } // quad rule
      }
    }

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

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
        double sum = xgeom.kappa[NEG] * b_t_neg + xgeom.kappa[POS] * b_t_pos;
        kappa_neg = xgeom.kappa[NEG] * b_t_neg / sum;
        kappa_pos = xgeom.kappa[POS] * b_t_pos / sum;
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
        kappa_neg = xgeom.kappa[NEG];
        kappa_pos = xgeom.kappa[POS];
        break;	      
      }
    }

    const double lam = minimal_stabilization ? 0.0 : lambda->EvaluateConst();


    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;

      const double weight = fquad.weights(i) * len; 
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);
      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
        
      shape = scafe->GetShape(mip.IP(), lh);
      jump.Range(0,ndof_h1) = (b_pos-b_neg) * shape;
      jump.Range(ndof_h1,ndof) = shape;

      scafe->CalcMappedDShape (mip,dshape_h1);

      dshape.Range(0,ndof_h1) = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape_h1 * normal);
      dshape.Range(ndof_h1,ndof) = dshape_h1 * normal;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l] == NEG){
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

      Nc -= weight * jump * Trans(dshape);
      Ns += ava * weight * jump * Trans(jump);

    }

    if (minimal_stabilization)
    {
      Lsys = 0.0;

      Lsys.Cols(0,ndof).Rows(0,ndof) = A;
      Lsys.Cols(ndof,ndof+2).Rows(0,ndof) = constrb;
      Lsys.Rows(ndof,ndof+2).Cols(0,ndof) = Trans(constrb);
      LapackInverse (Lsys);

      L = Lsys.Cols(0,ndof).Rows(0,ndof) * Trans(Nc);

      elmat = Nc + Trans(Nc) + 1.0 * /*lam*(p+1)**/ p/h * Ns; 

      elmat += 1.5 * Trans(L) * A * L;
    }
    else
      elmat = Nc + Trans(Nc) + lam*(p+1)/p/h * Ns; 

  }

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void SpaceTimeXNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXNitscheIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

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

    int ndof_x = xfe->GetNDof();
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> jump(ndof,lh);

    FlatVector<> shape_total(ndof,lh);
    FlatVector<> shape(ndof_h1,&shape_total(0));
    FlatVector<> shape_x(ndof_x,&shape_total(ndof_h1));

    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<D> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);

    FlatMatrix<> Nc(ndof,ndof,lh);
    FlatMatrix<> Ns(ndof,ndof,lh);
    FlatMatrix<> A(ndof,ndof,lh);

    FlatMatrix<> Lsys(ndof+4,ndof+4,lh);
    FlatMatrix<> L(ndof,ndof,lh); // lifting matrix

    ablockintegrator-> CalcElementMatrix (base_fel,
                                          eltrans, 
                                          A, lh);
    Ns = 0.0;
    Nc = 0.0;

    FlatMatrixFixWidth<4> constrb(ndof,lh);
    constrb = 0.0;

    const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    int ps = scafe->OrderSpace();
    // int pt = scafe->OrderTime();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
    const FlatQuadratureRuleCoDim1<D+1> & fquad(fcompr.GetInterfaceRule());

    { // calc constrb
      DOMAIN_TYPE dt = POS;
      for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
      {
        const FlatQuadratureRule<D+1> & fdomquad(fcompr.GetRule(dt));
        for (int i = 0; i < fdomquad.Size(); ++i)
        {
          IntegrationPoint ip(&fdomquad.points(i,0),fdomquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);

          scafe->CalcShapeSpaceTime(ip, fdomquad.points(i,D), shape, lh);

          shape_x = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shape_x(l) = 0.0;
          }

          //constant in time, constant in space constraint
          constrb.Col(dt) += tau * mip.GetMeasure() * fdomquad.weights(i) * shape_total;
          //linear in time, constant in space constraint
          constrb.Col(dt+2) += tau * mip.GetMeasure() * fdomquad.weights(i) * fdomquad.points(i,D) * shape_total ;
          
        } // quad rule
      }
    }

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

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
        double sum = xgeom.kappa[NEG] * b_t_neg + xgeom.kappa[POS] * b_t_pos;
        kappa_neg = xgeom.kappa[NEG] * b_t_neg / sum;
        kappa_pos = xgeom.kappa[POS] * b_t_pos / sum;
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
        kappa_neg = xgeom.kappa[NEG];
        kappa_pos = xgeom.kappa[POS];
        break;	      
      }
    }

    const double lam = minimal_stabilization ? 0.0 : lambda->EvaluateConst();

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      const double time = fquad.points(i,D);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      // const double h = D == 2 ? sqrt(absdet) : cbrt(absdet);

      Vec<D> nref_space; 
      for (int d = 0; d < D; ++d) 
        nref_space[d] = fquad.normals(i,d);
      Vec<D> normal_space = tau * absdet * Trans(Finv) * nref_space;
      double n_t = fquad.normals(i,D) * absdet;

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

      const double weight = fquad.weights(i) * len * nu;
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);
      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
        
      scafe->CalcShapeSpaceTime(mip.IP(), time, shape, lh);
      jump.Range(0,ndof_h1) = (b_pos-b_neg) * shape;
      jump.Range(ndof_h1,ndof) = shape;

      scafe->CalcMappedDxShapeSpaceTime(mip, time, dshape_h1, lh);

      dshape.Range(0,ndof_h1) = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape_h1 * normal_space);
      dshape.Range(ndof_h1,ndof) = dshape_h1 * normal_space;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l] == NEG){
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

      Nc -= weight * jump * Trans(dshape);
      Ns += ava * weight * jump * Trans(jump);

    }

    if (minimal_stabilization)
    {
      Lsys = 0.0;

      Lsys.Cols(0,ndof).Rows(0,ndof) = A;
      Lsys.Cols(ndof,ndof+4).Rows(0,ndof) = constrb;
      Lsys.Rows(ndof,ndof+4).Cols(0,ndof) = Trans(constrb);
      LapackInverse (Lsys);

      L = Lsys.Cols(0,ndof).Rows(0,ndof) * Trans(Nc);

      elmat = Nc + Trans(Nc) + 1.0 * /*lam*(p+1)**/ ps/h * Ns; 

      elmat += 1.5 * Trans(L) * A * L;
    }
    else
      elmat = Nc + Trans(Nc) + lam * (ps+1)*ps/h * Ns; 

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

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initx_min_stab_nitsche2d_1 ("xnitsche_minstab_halfhalf", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche2d_2 ("xnitsche_minstab_hansbo", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche2d_2b ("xnitsche_minstab", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBOBETA> > initx_min_stab_nitsche2d_3 ("xnitsche_minstab_hansbobeta", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::BETA> > initx_min_stab_nitsche2d_4 ("xnitsche_minstab_beta", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHA> > initx_min_stab_nitsche2d_5 ("xnitsche_minstab_alpha", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHABETA> > initx_min_stab_nitsche2d_6 ("xnitsche_minstab_alphabeta", 2, 4);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initx_min_stab_nitsche3d_1 ("xnitsche_minstab_halfhalf", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche3d_2 ("xnitsche_minstab_hansbo", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche3d_2b ("xnitsche_minstab", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBOBETA> > initx_min_stab_nitsche3d_3 ("xnitsche_minstab_hansbobeta", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::BETA> > initx_min_stab_nitsche3d_4 ("xnitsche_minstab_beta", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHA> > initx_min_stab_nitsche3d_5 ("xnitsche_minstab_alpha", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHABETA> > initx_min_stab_nitsche3d_6 ("xnitsche_minstab_alphabeta", 3, 4);

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

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > init_min_stab_xnitsche2d_st_1 ("stx_nitsche_min_stab_halfhalf", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche2d_st_2 ("stx_nitsche_min_stab_hansbo", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche2d_st_2b ("stx_nitsche_min_stab", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBOBETA> > init_min_stab_xnitsche2d_st_3 ("stx_nitsche_min_stab_hansbobeta", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::BETA> > init_min_stab_xnitsche2d_st_4 ("stx_nitsche_min_stab_beta", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHA> > init_min_stab_xnitsche2d_st_5 ("stx_nitsche_min_stab_alpha", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHABETA> > init_min_stab_xnitsche2d_st_6 ("stx_nitsche_min_stab_alphabeta", 2, 6);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > init_min_stab_xnitsche3d_st_1 ("stx_nitsche_min_stab_halfhalf", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche3d_st_2 ("stx_nitsche_min_stab_hansbo", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche3d_st_2b ("stx_nitsche_min_stab", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBOBETA> > init_min_stab_xnitsche3d_st_3 ("stx_nitsche_min_stab_hansbobeta", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::BETA> > init_min_stab_xnitsche3d_st_4 ("stx_nitsche_min_stab_beta", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHA> > init_min_stab_xnitsche3d_st_5 ("stx_nitsche_min_stab_alpha", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::ALPHABETA> > init_min_stab_xnitsche3d_st_6 ("stx_nitsche_min_stab_alphabeta", 3, 6);

}

/// coefficientfunction statt function-pointer
