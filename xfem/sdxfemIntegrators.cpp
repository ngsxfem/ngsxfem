#include "sdxfemIntegrators.hpp"

namespace ngfem
{


  template<int D>
  void SDXIntegrator<D> :: 
  AddPointContribution (const MappedIntegrationPoint<D,D> & mip,
                        FlatVector<> & shape,
                        FlatMatrixFixWidth<D> & dshape,
                        FlatMatrixFixWidth<D*D> & ddshape_ref_h1,
                        FlatVector<> & lapshape,
                        FlatVector<> & dudwshape,
                        FlatVector<> & diffopshape,
                        FlatMatrix<double> & elmat, 
                        int ndof_h1,
                        int ndof_x,
                        int order,
                        DOMAIN_TYPE dt,
                        double convmax,
                        const ScalarFiniteElement<D>* scafe,
                        const XFiniteElement * xfe) const
  {
    const int p = order;
    const int ndof = ndof_h1 + ndof_x;

    FlatVector<> shape_h1(ndof_h1,&shape(0)); //flat overlay
    FlatVector<> shape_x(ndof_x,&shape(ndof_h1)); //flat overlay

    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,&dshape(0,0)); //flat overlay
    FlatMatrixFixWidth<D> dshape_x(ndof_x,&dshape(ndof_h1,0)); //flat overlay

    FlatVector<> lapshape_h1(ndof_h1,&lapshape(0)); //flat overlay
    FlatVector<> lapshape_x(ndof_x,&lapshape(ndof_h1)); //flat overlay

    FlatVector<> dudwshape_h1(ndof_h1,&dudwshape(0)); //flat overlay
    FlatVector<> dudwshape_x(ndof_x,&dudwshape(ndof_h1)); //flat overlay

    const double beta = dt == POS ? beta_pos->Evaluate(mip) : beta_neg->Evaluate(mip);
    const double alpha = dt == POS ? alpha_pos->Evaluate(mip) : alpha_neg->Evaluate(mip);
    const double alpha_av = xfe ? 0.5 * (alpha_pos->Evaluate(mip) + alpha_neg->Evaluate(mip)) : alpha;
    const double mass = dt == POS ? mass_pos->Evaluate(mip) : mass_neg->Evaluate(mip);
        
    Vec<D> conv;
    if (dt == POS)
      conv_pos->Evaluate(mip,conv);
    else
      conv_neg->Evaluate(mip,conv);

    const double vol = (D == 2? 0.5 : 1.0/6.0) * mip.GetJacobiDet();
    const double h = D == 2? sqrt(2.0*vol) : cbrt(6.0*vol);
    const double Pe = 0.5 * h/p * convmax / alpha_av;

    // if (Pe < 1) {skip = true; continue;}
    const double stabparam = Pe < 1? 0 : (1.0-1.0/Pe) * 0.5 * h /convmax;

    Mat<D> Finv_T = Trans(mip.GetJacobianInverse());
    Vector<> transf(D*D);
    transf(0) = Finv_T(0,0)*Finv_T(0,0)+Finv_T(1,0)*Finv_T(1,0);
    transf(1) = Finv_T(0,0)*Finv_T(0,1)+Finv_T(1,0)*Finv_T(1,1);
    transf(2) = Finv_T(0,0)*Finv_T(0,1)+Finv_T(1,0)*Finv_T(1,1);
    transf(3) = Finv_T(0,1)*Finv_T(0,1)+Finv_T(1,1)*Finv_T(1,1);
        
    scafe->CalcShape(mip.IP(), shape_h1);
    scafe->CalcMappedDShape (mip,dshape_h1);
    scafe->CalcDDShape (mip.IP(),ddshape_ref_h1);
          
    lapshape_h1 = ddshape_ref_h1 * transf;
    dudwshape_h1 = dshape * conv;

    if (xfe)
    {
      shape_x = shape_h1;
      // dshape_x = dshape_h1;
      lapshape_x = lapshape_h1;
      dudwshape_x = dudwshape_h1;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xfe->GetSignsOfDof()[l] != dt){
          shape_x(l) = 0.0;
          // dshape_x.Row(l) = 0.0;
          lapshape_x(l) = 0.0;
          dudwshape_x(l) = 0.0;
        }
      }
    }

    diffopshape = - alpha * lapshape + dudwshape + mass * shape;
 
    const double fac = mip.GetWeight();

    elmat += (beta*fac*stabparam) * dudwshape * Trans(diffopshape);
  }


  template<int D>
  void SDXIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("SDXIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarFiniteElement<D> * scafe = NULL;

    if (D==3) throw Exception ("DDshape is not correct for D==3");
    
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;

    FlatVector<> shape(ndof,lh);
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D*D> ddshape_ref_h1(ndof_h1,lh);
    FlatVector<> lapshape(ndof,lh);
    FlatVector<> dudwshape(ndof,lh);
    FlatVector<> diffopshape(ndof,lh);

    int p = scafe->Order();

    double convmax = 0;
    { // START estimate convmax;
      DOMAIN_TYPE dt = POS;
      for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
      {
        if(!xfe)
        { 
          if (dummfe->GetDomainType() != dt)
            continue;
          IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
          for (int i = 0 ; i < pir.GetNIP(); i++)
          {
            MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
            Vec<D> conv;
            if (dt == POS)
              conv_pos->Evaluate(mip,conv);
            else
              conv_neg->Evaluate(mip,conv);
            convmax = max( convmax, max(abs(conv(0)),abs(conv(1))) );
          }
        }
        else
        {
          const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
          for (int i = 0; i < fquad.Size(); ++i)
          {

            IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
            MappedIntegrationPoint<D,D> mip(ip, eltrans);
            Vec<D> conv;
            if (dt == POS)
              conv_pos->Evaluate(mip,conv);
            else
              conv_neg->Evaluate(mip,conv);
            convmax = max( convmax, max(abs(conv(0)),abs(conv(1))) );
          }
        }
      }
    } // END estimate convmax;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          AddPointContribution(mip,
                               shape, dshape, ddshape_ref_h1,
                               lapshape, dudwshape, diffopshape,
                               elmat, 
                               ndof_h1, ndof_x,
                               p, dt, convmax,
                               scafe, xfe);
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {

          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          AddPointContribution(mip,
                               shape, dshape, ddshape_ref_h1,
                               lapshape, dudwshape, diffopshape,
                               elmat, 
                               ndof_h1, ndof_x,
                               p, dt, convmax,
                               scafe, xfe);
        }
      }
    }
  }

/*

  template<int D>
  void XNitscheConvScaledIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement<D> * xfe = NULL;
    const XDummyFE<D> * dummfe = NULL;
    const ScalarFiniteElement<D> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement<D>* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE<D>* >(&cfel[i]);
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

    int p = scafe->Order();

    double convmax = 0;
    { // START estimate convmax;
      bool sign=false;
      for (int k=0; k<2; sign=!sign, k++)
      {
        IntegrationRule pir; //partial integration rule
        if(!xfe)
        { 
          if (dummfe->GetSign() != sign)
            continue;
          pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
        }
        else{
          xfe->GetMasterElement().FillVolumeIntegrationRule(2*p,sign,pir);
        }
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          Vec<D> conv;
          if (sign)
            conv_pos->Evaluate(mip,conv);
          else
            conv_neg->Evaluate(mip,conv);
          convmax = max( convmax, max(abs(conv(0)),abs(conv(1))) );
        }
      }
    } // END estimate convmax;



    int ndof_x = xfe->GetNDof();
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> jump(ndof,lh);
    FlatMatrixFixWidth<2> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<2> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);
    const Array<int>& xsign = xfe->GetSignsOfDof();
    const MasterElement & masterel = xfe->GetMasterElement();
    
    // bool sign=false;
    IntegrationRule iir; //partial integration rule
    Array<Vec<2> > normals;

    masterel.FillInterfaceIntegrationRule(2*p,iir,normals);
    
    Vec<2> kappa = masterel.CalcKappa();
    // Vec<2> kappa(0.5,0.5);
    double kappa_neg = kappa(0);
    double kappa_pos = kappa(1);

    for (int i = 0 ; i < iir.GetNIP(); i++)
    {
      MappedIntegrationPoint<D,D> mip(iir[i], eltrans);

      Mat<D> inv_jac = mip.GetJacobianInverse();
      const double det = mip.GetJacobiDet();

      Vec<D> normal = det * Trans (inv_jac) * normals[i];       
      const double len = L2Norm (normal);
      normal /= len; // showing from positive to negative!

      const double weight = iir[i].Weight() * len;
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);
      const double b_neg = beta_neg->Evaluate(mip);
      const double b_pos = beta_pos->Evaluate(mip);
      const double lam = lambda->Evaluate(mip);

      jump.Range(0,ndof_h1) = (b_pos-b_neg) * scafe->GetShape(mip.IP(), lh);
      jump.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);

      scafe->CalcMappedDShape (mip,dshape_h1);
      xfe->CalcMappedDShape (mip,dshape_x);

      dshape.Range(0,ndof_h1) = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape_h1 * normal);
      dshape.Range(ndof_h1,ndof) = dshape_x * normal;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xsign[l]){
          jump(ndof_h1+l) *= -b_neg;
          dshape(ndof_h1+l) *= kappa_neg * a_neg;
        }
        else{
          jump(ndof_h1+l) *= b_pos;
          dshape(ndof_h1+l) *= kappa_pos * a_pos;
        }
      }

      const double ava = a_pos*0.5+a_neg*0.5;

      elmat -= weight * jump * Trans(dshape);
      elmat -= weight * dshape * Trans(jump);

      const double stabparam_diff = (p+1)*p/len * ava;
      const double stabparam_conv = convmax;
      const double stabparam = max(stabparam_diff,stabparam_conv);
      elmat += lam * stabparam * weight * jump * Trans(jump);

    }

  }
*/


  template<int D>
  void SDXSourceIntegrator<D> :: 
  AddPointContribution (const MappedIntegrationPoint<D,D> & mip,
                        FlatMatrixFixWidth<D> & dshape,
                        FlatVector<> & dudwshape,
                        FlatVector<double> & elvec, 
                        int ndof_h1,
                        int ndof_x,
                        int order,
                        DOMAIN_TYPE dt,
                        double convmax,
                        const ScalarFiniteElement<D>* scafe,
                        const XFiniteElement * xfe) const
  {
    const int p = order;
    const int ndof = ndof_h1 + ndof_x;

    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,&dshape(0,0)); //flat overlay

    FlatVector<> dudwshape_h1(ndof_h1,&dudwshape(0)); //flat overlay
    FlatVector<> dudwshape_x(ndof_x,&dudwshape(ndof_h1)); //flat overlay

    const double beta = dt == POS ? beta_pos->Evaluate(mip) : beta_neg->Evaluate(mip);
    const double alpha = dt == POS ? alpha_pos->Evaluate(mip) : alpha_neg->Evaluate(mip);
    const double alpha_av = xfe ? 0.5 * (alpha_pos->Evaluate(mip) + alpha_neg->Evaluate(mip)) : alpha;
    const double mass = dt == POS ? mass_pos->Evaluate(mip) : mass_neg->Evaluate(mip);

    const double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
        
    Vec<D> conv;
    if (dt == POS)
      conv_pos->Evaluate(mip,conv);
    else
      conv_neg->Evaluate(mip,conv);

    const double vol = (D == 2? 0.5 : 1.0/6.0) * mip.GetJacobiDet();
    const double h = D == 2? sqrt(2.0*vol) : cbrt(6.0*vol);
    const double Pe = 0.5 * h/p * convmax / alpha_av;

    // if (Pe < 1) {skip = true; continue;}
    const double stabparam = Pe < 1? 0 : (1.0-1.0/Pe) * 0.5 * h /convmax;

    scafe->CalcMappedDShape (mip,dshape_h1);
          
    dudwshape_h1 = dshape_h1 * conv;

    if (xfe)
    {
      dudwshape_x = dudwshape_h1;

      for (int l = 0; l < ndof_x; ++l)
      {
        if (xfe->GetSignsOfDof()[l] != dt){
          dudwshape_x(l) = 0.0;
        }
      }
    }

    const double fac = mip.GetWeight();

    elvec += (beta*fac*stabparam*coef) * dudwshape;
  }


  template<int D>
  void SDXSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    static Timer timer ("SDXSourceIntegrator::CalcElementVector");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarFiniteElement<D> * scafe = NULL;

    if (D==3) throw Exception ("DDshape is not correct for D==3");
    
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;

    FlatVector<> shape(ndof,lh);
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D*D> ddshape_ref_h1(ndof_h1,lh);
    FlatVector<> lapshape(ndof,lh);
    FlatVector<> dudwshape(ndof,lh);
    FlatVector<> diffopshape(ndof,lh);

    int p = scafe->Order();

    double convmax = 0;
    { // START estimate convmax;
      DOMAIN_TYPE dt = POS;
      for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
      {
        if(!xfe)
        { 
          if (dummfe->GetDomainType() != dt)
            continue;
          IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
          for (int i = 0 ; i < pir.GetNIP(); i++)
          {
            MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
            Vec<D> conv;
            if (dt == POS)
              conv_pos->Evaluate(mip,conv);
            else
              conv_neg->Evaluate(mip,conv);
            convmax = max( convmax, L2Norm(conv)); //max(abs(conv(0)),abs(conv(1))) );
          }
        }
        else
        {
          const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
          for (int i = 0; i < fquad.Size(); ++i)
          {

            IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
            MappedIntegrationPoint<D,D> mip(ip, eltrans);
            Vec<D> conv;
            if (dt == POS)
              conv_pos->Evaluate(mip,conv);
            else
              conv_neg->Evaluate(mip,conv);
            convmax = max( convmax, L2Norm(conv)); //max(abs(conv(0)),abs(conv(1))) );
          }
        }
      }
    } // END estimate convmax;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          AddPointContribution(mip,
                               dshape, dudwshape, 
                               elvec, 
                               ndof_h1, ndof_x,
                               p, dt, convmax,
                               scafe, xfe);
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {

          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          AddPointContribution(mip,
                               dshape, dudwshape, 
                               elvec, 
                               ndof_h1, ndof_x,
                               p, dt, convmax,
                               scafe, xfe);
        }
      }
    }
  }



  template class SDXIntegrator<2>;
  template class SDXIntegrator<3>;

  template class SDXSourceIntegrator<2>;
  template class SDXSourceIntegrator<3>;

  static RegisterBilinearFormIntegrator<SDXIntegrator<2> > initsdxcut2d ("sdstab", 2, 8);
  static RegisterBilinearFormIntegrator<SDXIntegrator<3> > initsdxcut3d ("sdstab", 3, 8);

  // static RegisterBilinearFormIntegrator<XNitscheConvScaledIntegrator<2> > initxnitscheconv2d ("xnitscheconv", 2, 7);

  static RegisterLinearFormIntegrator<SDXSourceIntegrator<2> > initsdxsource2d ("sdxsource", 2, 10);
  static RegisterLinearFormIntegrator<SDXSourceIntegrator<3> > initsdxsource3d ("sdxsource", 3, 10);

}

/// coefficientfunction statt function-pointer
