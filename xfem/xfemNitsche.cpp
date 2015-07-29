#include "xfemNitsche.hpp"

namespace ngfem
{

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice, NITSCHE_VARIANTS::SCALING_CHOICE scale_choice>
  void XNitscheIntegrator<D, kappa_choice, scale_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XNitscheIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;
    if (!xfe) return;

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

    ablockintegrator-> CalcElementMatrix (base_fel,eltrans,A, lh);
    Ns = 0.0;
    Nc = 0.0;

    FlatMatrixFixWidth<2> constrb(ndof,lh);
    constrb = 0.0;

    const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    int p = scafe->Order();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    double convmax=0.0;
    if (scale_choice == NITSCHE_VARIANTS::CONVECTIVE)
    { // START estimate convmax;
      for (auto dt : {POS,NEG})
      {
        const FlatQuadratureRule<D> & fquaddt(fcompr.GetRule(dt));
        for (int i = 0; i < fquaddt.Size(); ++i)
        {
          IntegrationPoint ip(&fquaddt.points(i,0),fquaddt.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          Vec<D> conv;
          if (dt == POS)
            conv_pos->Evaluate(mip,conv);
          else
            conv_neg->Evaluate(mip,conv);
          convmax = max( convmax, max(abs(conv(0)),abs(conv(1))) );
        }
      }
    } // END estimate convmax;


    { // calc constrb
      for (auto dt : {POS,NEG})
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

    const double a_t_neg = alpha_neg->Evaluate(mipc);
    const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;

    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HEAVISIDE:
      {
        if (xgeom.kappa[NEG] >= 0.5)
            kappa_neg = 1.0;
        else
            kappa_neg = 0.0;
        kappa_pos = 1.0 - kappa_neg;
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

    double ava = a_t_pos*0.5+a_t_neg*0.5;

    const double Pe = 0.5 * h/p * convmax / ava;

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

      Nc -= weight * jump * Trans(dshape);
      Ns += weight * jump * Trans(jump);

    }

    double avah = ava/h;
    
    if (scale_choice == NITSCHE_VARIANTS::CONVECTIVE)
    {
      avah *= max(1.0, Pe);
    }

    if (minimal_stabilization)
    {
      Lsys = 0.0;

      Lsys.Cols(0,ndof).Rows(0,ndof) = A;
      Lsys.Cols(ndof,ndof+2).Rows(0,ndof) = constrb;
      Lsys.Rows(ndof,ndof+2).Cols(0,ndof) = Trans(constrb);
      LapackInverse (Lsys);

      L = Lsys.Cols(0,ndof).Rows(0,ndof) * Trans(Nc);

      elmat = Nc + Trans(Nc) + avah * 1.0 * /*lam*(p+1)**/ p * Ns; 

      elmat += 1.5 * Trans(L) * A * L;
    }
    else
      elmat = Nc + Trans(Nc) + lam*(p+1)*p * avah * Ns; 
  }


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void XNitscheRhsJumpIntegrator<D, kappa_choice> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
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

    elvec = 0.0;
    
    // if (D==3)
    //   throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

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

    FlatVector<> Nc(ndof,lh);
    FlatVector<> Ns(ndof,lh);
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

    // const double b_t_neg = beta_neg->Evaluate(mipc);
    // const double b_t_pos = beta_pos->Evaluate(mipc);
    // const double a_t_neg = alpha_neg->Evaluate(mipc);
    // const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;

    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HEAVISIDE:
      {
        if (xgeom.kappa[NEG] >= 0.5)
          { 
            kappa_neg = 1.0;
            kappa_pos = 0.0;
          }
        else
          { 
            kappa_neg = 0.0;
            kappa_pos = 1.0;
          }
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
      
      double rhs = coef_rhs->Evaluate(mip);

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


      double ava = a_pos*0.5+a_neg*0.5;

      Nc -= weight * rhs * dshape;
      Ns += ava * weight * rhs * jump;

    }

    if (minimal_stabilization)
    {
      throw Exception("not yet minstabable");
      /*
      Lsys = 0.0;

      Lsys.Cols(0,ndof).Rows(0,ndof) = A;
      Lsys.Cols(ndof,ndof+2).Rows(0,ndof) = constrb;
      Lsys.Rows(ndof,ndof+2).Cols(0,ndof) = Trans(constrb);
      LapackInverse (Lsys);

      L = Lsys.Cols(0,ndof).Rows(0,ndof) * Trans(Nc);
      */
      // elvec = Nc + Trans(Nc) + 1.0 * /*lam*(p+1)**/ p/h * Ns; 

      // elvec += 1.5 * Trans(L) * A * L;
    }
    else
      elvec = Nc + lam*(p+1)*p/h * Ns;
      // elvec = Nc + Trans(Nc) + lam*(p+1)/p/h * Ns; 
      // elmat = lam*(p+1)/p/h * Ns; 

  }

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void XNitscheRhsFluxJumpIntegrator<D, kappa_choice> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
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

    elvec = 0.0;
    
    // if (D==3)
    //   throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

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

    FlatVector<> Ns(ndof,lh);
    Ns = 0.0;

    const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();
    // int p = scafe->Order();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    // const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

    // const double b_t_neg = beta_neg->Evaluate(mipc);
    // const double b_t_pos = beta_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;

    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HEAVISIDE:
      {
        if (xgeom.kappa[NEG] >= 0.5)
          { 
            kappa_neg = 1.0;
            kappa_pos = 0.0;
          }
        else
          { 
            kappa_neg = 0.0;
            kappa_pos = 1.0;
          }
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

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
      double rhs = coef_rhs->Evaluate(mip);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;

      const double weight = fquad.weights(i) * len; 

      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
        
      shape = scafe->GetShape(mip.IP(), lh);
      jump.Range(0,ndof_h1) = (b_pos*kappa_neg+b_neg*kappa_pos) * shape;
      jump.Range(ndof_h1,ndof) = shape;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l] == NEG){
          jump(ndof_h1+l) *= b_neg*kappa_pos;
        }
        else{
          jump(ndof_h1+l) *= b_pos*kappa_neg;
        }
      }

      Ns += weight * rhs * jump;
    }
    elvec = Ns;
  }


  template <int D>
  void CalcSTCrazyKappaCoeffs( const FlatXLocalGeometryInformation & xgeom, const ElementTransformation & eltrans, const double tau, double & kappa_neg, double & kappa_pos, LocalHeap & lh)

  {
    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);
    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

    Mat<2> mass_i[2];
    mass_i[NEG] *= 0.0;
    mass_i[POS] *= 0.0;

    double gamma[2];

    double int_dom_1[2];
    double int_dom_2[2];
    double int_if_1 = 0;
    double int_if_2 = 0;
    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      int_dom_1[(int)dt] = 0;
      int_dom_2[(int)dt] = 0;
      const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
      const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));
      for (int i = 0; i < fquad.Size(); ++i)
      {
        IntegrationPoint ips;
        for (int d = 0; d < D; ++d)
          ips(d) = fquad.points(i,d);
        MappedIntegrationPoint<D,D> mip(ips, eltrans);

        const double psi1sq = sqr(fquad.points(i,D));
        const double psi2sq = sqr(1-fquad.points(i,D));
        int_dom_1[(int)dt] += mip.GetMeasure() * fquad.weights(i) * psi1sq;
        int_dom_2[(int)dt] += mip.GetMeasure() * fquad.weights(i) * psi2sq;
        mass_i[dt](0,0) += 2.0 * fquad.weights(i) * (1-fquad.points(i,D)) * (1-fquad.points(i,D));
        mass_i[dt](1,0) += 2.0 * fquad.weights(i) * fquad.points(i,D) * (1-fquad.points(i,D));
        mass_i[dt](0,1) += 2.0 * fquad.weights(i) * fquad.points(i,D) * (1-fquad.points(i,D));
        mass_i[dt](1,1) += 2.0 * fquad.weights(i) * fquad.points(i,D) * fquad.points(i,D);
      }

      FlatMatrix<> Msys(2,2,&mass_i[dt](0,0));

      // std::cout << " Msys = " << Msys << std::endl;
      LapackInverse(Msys);
      // std::cout << " Msys = " << Msys << std::endl;
      
      const double alpha = 4*Msys(0,0) + 2*Msys(0,1) + 2*Msys(1,0) + Msys(1,1);
      gamma[dt] = 4.0 /alpha;
    }

    if (false) {
      std::cout << " gamma[NEG] = " << gamma[NEG] << std::endl;
      std::cout << " gamma[POS] = " << gamma[POS] << std::endl;
      // getchar();
    }


    if (false) {
      std::cout << " gamma[NEG] + gamma[POS] = " << gamma[NEG] + gamma[POS] << std::endl;
      // getchar();
    }

    // if (gamma[NEG] < gamma[POS])
    // {
    //   kappa_neg = gamma[NEG];
    //   kappa_pos = 1.0 - kappa_neg;
    // }
    // else
    // {
    //   kappa_pos = gamma[POS];
    //   kappa_neg = 1.0 - kappa_pos;
    // }

    // kappa_pos = 0.5;
    // kappa_neg = 0.5;

    return;


    const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
    const FlatQuadratureRuleCoDim1<D+1> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      // const double time = fquad.points(i,D);
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

      const double psi1sq = sqr(fquad.points(i,D));
      const double psi2sq = sqr(1-fquad.points(i,D));

      int_if_1 += weight * psi1sq;
      int_if_2 += weight * psi2sq;

    }

    // std::cout << " int_if_1 = " << int_if_1 << std::endl;
    // std::cout << " int_if_2 = " << int_if_2 << std::endl;

    // getchar();

    const double qneg_1 = int_dom_1[NEG] / (h * int_if_1);
    const double qneg_2 = int_dom_2[NEG] / (h * int_if_2);
    const double qneg = max(qneg_1,qneg_2);
    const double qpos_1 = int_dom_1[POS] / (h * int_if_1);
    const double qpos_2 = int_dom_2[POS] / (h * int_if_2);
    const double qpos = max(qpos_1,qpos_2);


    // if (qneg > qpos)
    // {
    //   kappa_neg = 1.0;
    //   kappa_pos = 0.0;
    // }
    // else
    // {
    //   kappa_neg = 1.0;
    //   kappa_pos = 0.0;
    // }

    if (qneg < 1 || qpos < 1)
    {
      std::cout << " int_dom_1[POS] = " << int_dom_1[POS] << std::endl;
      std::cout << " int_dom_1[NEG] = " << int_dom_1[NEG] << std::endl;
      std::cout << " int_dom_2[POS] = " << int_dom_2[POS] << std::endl;
      std::cout << " int_dom_2[NEG] = " << int_dom_2[NEG] << std::endl;

      // getchar();

      std::cout << " qneg = " << qneg << std::endl;
      std::cout << " qpos = " << qpos << std::endl;

      if (qneg < 1)
      {
        kappa_neg = 0.01*qneg;
        kappa_pos = 1 - kappa_neg;
      }
      else
      {
        kappa_pos = 0.01*qpos;
        kappa_neg = 1 - kappa_pos;
      }
      getchar();
    }

    // kappa_neg = 1.0 / (1.0 + sqrt(qpos/qneg));
    // kappa_pos = 1.0 / (1.0 + sqrt(qneg/qpos));

    // kappa_neg = qneg;
    // kappa_pos = 1.0 / (1.0 + sqrt(qneg/qpos));

    // std::cout << " kappa_neg = " << kappa_neg << std::endl;
    // std::cout << " kappa_pos = " << kappa_pos << std::endl;

    // std::cout << " kappa_neg * kappa_neg / qneg = " << kappa_neg * kappa_neg / qneg << std::endl;
    // std::cout << " kappa_pos * kappa_pos / qpos = " << kappa_pos * kappa_pos / qpos << std::endl;

    // getchar();

  }


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void SpaceTimeXNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
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
    
    // if (D==3)
    //   throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

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

    // const double b_t_neg = beta_neg->Evaluate(mipc);
    // const double b_t_pos = beta_pos->Evaluate(mipc);
    // const double a_t_neg = alpha_neg->Evaluate(mipc);
    // const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;

    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HEAVISIDE:
      {
        if (xgeom.kappa[NEG] >= 0.5)
          { 
            kappa_neg = 1.0;
            kappa_pos = 0.0;
          }
        else
          { 
            kappa_neg = 0.0;
            kappa_pos = 1.0;
          }
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

    CalcSTCrazyKappaCoeffs<D>(xgeom,eltrans,tau, kappa_neg, kappa_pos, lh);

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

      double ava = a_pos*0.5+a_neg*0.5;

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
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > initxnitsche2d_3 ("xnitsche_heaviside", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_2b ("xnitsche", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO,NITSCHE_VARIANTS::CONVECTIVE> > initxnitscheconv2d_2b ("xnitsche_conv", 2, 7);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initxnitsche3d_1 ("xnitsche_halfhalf", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_2 ("xnitsche_hansbo", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HEAVISIDE> > initxnitsche3d_3 ("xnitsche_heaviside", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_2b ("xnitsche", 3, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO,NITSCHE_VARIANTS::CONVECTIVE> > initxnitscheconv3d_2b ("xnitsche_conv", 3, 7);

  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitscherhsjump2d_1 ("xnitscherhsjump_halfhalf", 2, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsjump2d_2 ("xnitscherhsjump_hansbo", 2, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > initxnitscherhsjump2d_3 ("xnitscherhsjump_heaviside", 2, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsjump2d_2b ("xnitscherhsjump", 2, 6);

  static RegisterLinearFormIntegrator<XNitscheRhsFluxJumpIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitscherhsfluxjump2d_1 ("xnitscherhsfluxjump_halfhalf", 2, 3);
  static RegisterLinearFormIntegrator<XNitscheRhsFluxJumpIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsfluxjump2d_2 ("xnitscherhsfluxjump_hansbo", 2, 3);
  static RegisterLinearFormIntegrator<XNitscheRhsFluxJumpIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > initxnitscherhsfluxjump2d_3 ("xnitscherhsfluxjump_heaviside", 2, 3);
  static RegisterLinearFormIntegrator<XNitscheRhsFluxJumpIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsfluxjump2d_2b ("xnitscherhsfluxjump", 2, 3);
                                              
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initxnitscherhsjump3d_1 ("xnitscherhsjump_halfhalf", 3, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsjump3d_2 ("xnitscherhsjump_hansbo", 3, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<3,NITSCHE_VARIANTS::HEAVISIDE> > initxnitscherhsjump3d_3 ("xnitscherhsjump_heaviside", 3, 6);
  static RegisterLinearFormIntegrator<XNitscheRhsJumpIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitscherhsjump3d_2b ("xnitscherhsjump", 3, 6);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initx_min_stab_nitsche2d_1 ("xnitsche_minstab_halfhalf", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche2d_2 ("xnitsche_minstab_hansbo", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > initx_min_stab_nitsche2d_3 ("xnitsche_minstab_heaviside", 2, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche2d_2b ("xnitsche_minstab", 2, 4);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initx_min_stab_nitsche3d_1 ("xnitsche_minstab_halfhalf", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche3d_2 ("xnitsche_minstab_hansbo", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HEAVISIDE> > initx_min_stab_nitsche3d_3 ("xnitsche_minstab_heaviside", 3, 4);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initx_min_stab_nitsche3d_2b ("xnitsche_minstab", 3, 4);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitsche2d_st_1 ("stx_nitsche_halfhalf", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_st_2 ("stx_nitsche_hansbo", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > initxnitsche2d_st_3 ("stx_nitsche_heaviside", 2, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_st_2b ("stx_nitsche", 2, 7);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > initxnitsche3d_st_1 ("stx_nitsche_halfhalf", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_st_2 ("stx_nitsche_hansbo", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HEAVISIDE> > initxnitsche3d_st_3 ("stx_nitsche_heaviside", 3, 7);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > initxnitsche3d_st_2b ("stx_nitsche", 3, 7);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > init_min_stab_xnitsche2d_st_1 ("stx_nitsche_min_stab_halfhalf", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche2d_st_2 ("stx_nitsche_min_stab_hansbo", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HEAVISIDE> > init_min_stab_xnitsche2d_st_3 ("stx_nitsche_min_stab_heaviside", 2, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche2d_st_2b ("stx_nitsche_min_stab", 2, 6);

  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HALFHALF> > init_min_stab_xnitsche3d_st_1 ("stx_nitsche_min_stab_halfhalf", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche3d_st_2 ("stx_nitsche_min_stab_hansbo", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HEAVISIDE> > init_min_stab_xnitsche3d_st_3 ("stx_nitsche_min_stab_heaviside", 3, 6);
  static RegisterBilinearFormIntegrator<SpaceTimeXNitscheIntegrator<3,NITSCHE_VARIANTS::HANSBO> > init_min_stab_xnitsche3d_st_2b ("stx_nitsche_min_stab", 3, 6);

  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void FictXNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("FictXNitscheIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe[2]; xfe[0] = NULL; xfe[1] = NULL;
    // const XDummyFE * dummfe[2]; dummfe[0] = NULL; dummfe[1] = NULL;

    const ScalarFiniteElement<D> * scafe = NULL;

    xfe[0] = dynamic_cast<const XFiniteElement* >(&cfel[0]);
    xfe[1] = dynamic_cast<const XFiniteElement* >(&cfel[1]);
    // dummfe[0] = dynamic_cast<const XDummyFE* >(&cfel[0]);
    // dummfe[1] = dynamic_cast<const XDummyFE* >(&cfel[1]);

    elmat = 0.0;
    if (!xfe[0] || !xfe[1])
      return;

    if (xfe[0] != NULL)
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[0]->GetBaseFE()));
    else
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[1]->GetBaseFE()));

    
    // if (D==3)
    //   throw Exception(" D==3: len is not scaling correctly (h^2 instead of h)");

    int p = scafe->Order();

    const FlatXLocalGeometryInformation & xgeom(xfe[0]->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    int ndof_sca = scafe->GetNDof();
    int ndof_neg = xfe[0]!=NULL ? scafe->GetNDof() : 0;
    int ndof_pos = xfe[1]!=NULL ? scafe->GetNDof() : 0;
    int ndof_total = ndof_neg+ndof_pos;

    FlatVector<> shape_sca(ndof_sca,lh);

    FlatVector<> jump(ndof_total,lh);
    FlatVector<> dnshape(ndof_total,lh);

    FlatMatrixFixWidth<D> dshape_sca(ndof_sca,lh);

    FlatMatrix<> elmat_neg(ndof_neg,ndof_neg,lh);
    FlatMatrix<> elmat_pos(ndof_pos,ndof_pos,lh);

    FlatMatrix<> Nc(ndof_total,ndof_total,lh);
    FlatMatrix<> Ns(ndof_total,ndof_total,lh);

    Ns = 0.0;
    Nc = 0.0;

    /*
    FlatMatrix<> A(ndof_total,ndof_total,lh);

    FlatMatrix<> Lsys(ndof_total+2,ndof_total+2,lh);
    FlatMatrix<> L(ndof_total,ndof_total,lh); // lifting matrix

    ablockintegrator-> CalcElementMatrix (base_fel,
                                          eltrans, 
                                          A, lh);

    FlatMatrixFixWidth<2> constrb(ndof_total,lh);
    constrb = 0.0;

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
    */

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

    // const double b_t_neg = beta_neg->Evaluate(mipc);
    // const double b_t_pos = beta_pos->Evaluate(mipc);
    // const double a_t_neg = alpha_neg->Evaluate(mipc);
    // const double a_t_pos = alpha_pos->Evaluate(mipc);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;

    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HEAVISIDE:
      {
        if (xgeom.kappa[NEG] >= 0.5)
          { 
            kappa_neg = 1.0;
            kappa_pos = 0.0;
          }
        else
          { 
            kappa_neg = 0.0;
            kappa_pos = 1.0;
          }
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
        
      shape_sca = scafe->GetShape(mip.IP(), lh);

      jump.Range(0,ndof_neg) = (-b_neg) * shape_sca;
      jump.Range(ndof_neg,ndof_total) = (b_pos) * shape_sca;

      scafe->CalcMappedDShape (mip,dshape_sca);

      dnshape.Range(0,ndof_neg) = (a_neg*kappa_neg) * (dshape_sca * normal);
      dnshape.Range(ndof_neg,ndof_total) = (a_pos*kappa_pos) * (dshape_sca * normal);

      double ava = a_pos*0.5+a_neg*0.5;

      Nc -= weight * jump * Trans(dnshape);
      Ns += ava * weight * jump * Trans(jump);

    }

    if (minimal_stabilization)
    {
      throw Exception("not implemented");
      // Lsys = 0.0;

      // Lsys.Cols(0,ndof).Rows(0,ndof) = A;
      // Lsys.Cols(ndof,ndof+2).Rows(0,ndof) = constrb;
      // Lsys.Rows(ndof,ndof+2).Cols(0,ndof) = Trans(constrb);
      // LapackInverse (Lsys);

      // L = Lsys.Cols(0,ndof).Rows(0,ndof) * Trans(Nc);

      // elmat = Nc + Trans(Nc) + 1.0 * /*lam*(p+1)**/ p/h * Ns; 

      // elmat += 1.5 * Trans(L) * A * L;
    }
    else
      elmat = Nc + Trans(Nc) + lam*(p+1)*p/h * Ns; 
      // elmat = lam*(p+1)/p/h * Ns; 

  }

  static RegisterBilinearFormIntegrator<FictXNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initfictxnitsche2d_2b ("fictxnitsche", 2, 5);


}

/// coefficientfunction statt function-pointer
