#include "xfemIntegrators.hpp"

namespace ngfem
{


  template<int D>
  void XMassIntegrator<D> ::
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof,&shape_total(ndof));

    int p = scafe->Order();
    
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
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          shape.Range(0,ndof) = scafe->GetShape(mip.IP(), lh);
          double fac = mip.GetWeight();
          elmat += (fac*coef) * shape * Trans(shape);
        }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D> * compr (lset_eval_p->GetCompositeRule<D>());
        const QuadratureRule<D> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ip;
          for (int d = 0; d < D; ++d)
            ip(d) = quad.points[i](d);
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * quad.weights[i];
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template<int D>
  void XSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
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
    
    elvec = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof,&shape_total(ndof));

    int p = scafe->Order();
    
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
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          shape = scafe->GetShape(mip.IP(), lh);
          double fac = mip.GetWeight();
          elvec += (fac*coef) * shape;
        }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D> * compr (lset_eval_p->GetCompositeRule<D>());
        const QuadratureRule<D> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ip;
          for (int d = 0; d < D; ++d)
            ip(d) = quad.points[i](d);
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * quad.weights[i];
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D>
  void SpaceTimeXMassIntegrator<D> ::
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof,&shape_total(ndof));

    int ps = scafe->OrderSpace();
    int pt = scafe->OrderTime();
    
    const double t1 = coef_tnew->EvaluateConst(); 
    const double t0 = coef_told->EvaluateConst();
    const double tau = t1 - t0;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule irs = SelectIntegrationRule (eltrans.GetElementType(), 2*ps);
        IntegrationRule irt = SelectIntegrationRule (ET_SEGM, 2*pt);
        for (int k = 0; k < irt.GetNIP(); ++k)
          for (int l = 0 ; l < irs.GetNIP(); l++)
          {
            MappedIntegrationPoint<D,D> mip(irs[l], eltrans);
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
              : coef_neg->Evaluate(mip);
            scafe->CalcShapeSpaceTime(mip.IP(), irt[k](0), shape, lh);
            double fac = mip.GetWeight() * irt[k].Weight() * tau;
            elmat += (fac*coef) * shape * Trans(shape);
          }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D+1> * compr (lset_eval_p->GetCompositeRule<D+1>());
        const QuadratureRule<D+1> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = quad.points[i](d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, quad.points[i](D), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * quad.weights[i];
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template<int D>
  void SpaceTimeXSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
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
    
    elvec = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof,&shape_total(ndof));

    int ps = scafe->OrderSpace();
    int pt = scafe->OrderTime();
    
    const double t1 = coef_tnew->EvaluateConst(); 
    const double t0 = coef_told->EvaluateConst();
    const double tau = t1 - t0;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule irs = SelectIntegrationRule (eltrans.GetElementType(), 2*ps);
        IntegrationRule irt = SelectIntegrationRule (ET_SEGM, 2*pt);
        for (int k = 0; k < irt.GetNIP(); ++k)
          for (int l = 0 ; l < irs.GetNIP(); l++)
          {
            MappedIntegrationPoint<D,D> mip(irs[l], eltrans);
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
              : coef_neg->Evaluate(mip);
            scafe->CalcShapeSpaceTime(mip.IP(), irt[k](0), shape, lh);
            double fac = mip.GetWeight() * irt[k].Weight() * tau;
            elvec += (fac*coef) * shape;
          }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_p = xfe->GetLocalGeometry();
        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D+1> * compr (lset_eval_p->GetCompositeRule<D+1>());
        const QuadratureRule<D+1> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = quad.points[i](d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, quad.points[i](D), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * quad.weights[i];
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }



/*  
  template<int D>
  void XLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
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
    
    if(!xfe && !dummfe)
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,&dshape(0,0)); //flat overlay
    FlatMatrixFixWidth<D> dshape_x(ndof_x,&dshape(ndof_h1,0)); //flat overlay
    int p = scafe->Order();
    
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
        double coef = sign? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

        scafe->CalcMappedDShape (mip,dshape_h1);
        if (xfe)
        {
          xfe->CalcMappedDShape (mip,dshape_x);
          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] == sign)
              dshape.Row(ndof_h1+l) = 0.0;
          }
        }

        double fac = mip.GetWeight();
        elmat += (fac*coef) * dshape * Trans(dshape);
      }
    }
  }

  template<int D>
  void XConvectionIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
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
    
    if(!xfe && !dummfe)
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D> dshape_h1(ndof_h1,&dshape(0,0)); //flat overlay
    FlatMatrixFixWidth<D> dshape_x(ndof_x,&dshape(ndof_h1,0)); //flat overlay
    FlatVector<> shape(ndof,lh);
    FlatVector<> dudwshape(ndof,lh);

    int p = scafe->Order();
    
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
          coef_pos->Evaluate(mip,conv);
        else
          coef_neg->Evaluate(mip,conv);

        scafe->CalcMappedDShape (mip,dshape_h1);
        shape.Range(0,ndof_h1) = scafe->GetShape(mip.IP(), lh);
        if (xfe)
        {
          shape.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);
          xfe->CalcMappedDShape (mip,dshape_x);
          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] == sign){
              dshape.Row(ndof_h1+l) = 0.0;
              shape(ndof_h1+l) = 0.0;
            }
          }
        }
        dudwshape = dshape * conv;
        
        double fac = mip.GetWeight();
        elmat += fac * shape * Trans(dudwshape);
      }
    }
  }


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  void XNitscheIntegrator<D, kappa_choice> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
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

    int ndof_x = xfe->GetNDof();
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> jump(ndof,lh);
    FlatMatrixFixWidth<2> dshape_h1(ndof_h1,lh);
    FlatMatrixFixWidth<2> dshape_x(ndof_x,lh);
    FlatVector<> dshape(ndof,lh);
    const Array<int>& xsign = xfe->GetSignsOfDof();
    int p = scafe->Order();
    const MasterElement & masterel = xfe->GetMasterElement();
    
    // bool sign=false;
    IntegrationRule iir; //partial integration rule
    Array<Vec<2> > normals;

    masterel.FillInterfaceIntegrationRule(2*p,iir,normals);
    

    MappedIntegrationPoint<D,D> mipt(iir[0], eltrans);
    const double b_t_neg = beta_neg->Evaluate(mipt);
    const double b_t_pos = beta_pos->Evaluate(mipt);
    const double a_t_neg = alpha_neg->Evaluate(mipt);
    const double a_t_pos = alpha_pos->Evaluate(mipt);

    double kappa_neg = 0.5;
    double kappa_pos = 0.5;
    switch (kappa_choice){
    case NITSCHE_VARIANTS::HALFHALF:
      break;
    case NITSCHE_VARIANTS::HANSBOBETA:
      {
        Vec<2> kappa = masterel.CalcKappa();
        double sum = kappa(0) * b_t_neg + kappa(1) * b_t_pos;
        kappa_neg = kappa(0) * b_t_neg / sum;
        kappa_pos = kappa(1) * b_t_pos / sum;
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
        Vec<2> kappa = masterel.CalcKappa();
        kappa_neg = kappa(0);
        kappa_pos = kappa(1);
        break;	      
      }
    }


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
      elmat += lam * (p+1)*p/len * ava * weight * jump * Trans(jump);

    }

  }



  template<int D>
  void XRobinIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement<D-1> * xfe = NULL;
    const XDummyFE<D-1> * dummfe = NULL;
    const ScalarFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement<D-1>* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE<D-1>* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarFiniteElement<D-1>* >(&cfel[i]);
    }
    
    elmat = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> shape(ndof,lh);
    int p = scafe->Order();
    
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
        xfe->GetMasterEdge().FillSurfaceIntegrationRule(2*p,sign,pir);
      }
      
      for (int i = 0 ; i < pir.GetNIP(); i++)
      {
        MappedIntegrationPoint<D-1,D> mip(pir[i], eltrans);
        double coef = sign? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
        
        shape.Range(0,ndof_h1) = scafe->GetShape(mip.IP(), lh);
        if (xfe)
        {
          shape.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);
          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] == sign)
              shape(ndof_h1+l) = 0.0;
          }
        }

        double fac = mip.GetWeight();
        elmat += (fac*coef) * shape * Trans(shape);
      }
    }
  }



  template<int D>
  void XSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
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
    
    elvec=0.0;
    
    if(!xfe && !dummfe)
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> shape(ndof,lh);
    int p = scafe->Order();
    
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
        double coef = sign? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
        
        shape.Range(0,ndof_h1) = scafe->GetShape(mip.IP(), lh);
        if (xfe)
        {
          shape.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);
          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] == sign)
              shape(ndof_h1+l) = 0.0;
          }
        }

        double fac = mip.GetWeight();
        elvec += (fac*coef) * shape;
      }
    }
  }

  template<int D>
  void XNeumannIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement<D-1> * xfe = NULL;
    const XDummyFE<D-1> * dummfe = NULL;
    const ScalarFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement<D-1>* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE<D-1>* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarFiniteElement<D-1>* >(&cfel[i]);
    }
    
    elvec = 0.0;
    
    if(!xfe && !dummfe)
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof_h1 = scafe->GetNDof();
    int ndof = ndof_h1+ndof_x;
    FlatVector<> shape(ndof,lh);
    int p = scafe->Order();
    
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
        xfe->GetMasterEdge().FillSurfaceIntegrationRule(2*p,sign,pir);
      }
      
      for (int i = 0 ; i < pir.GetNIP(); i++)
      {
        MappedIntegrationPoint<D-1,D> mip(pir[i], eltrans);
        double coef = sign? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
        
        shape.Range(0,ndof_h1) = scafe->GetShape(mip.IP(), lh);
        if (xfe)
        {
          shape.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);
          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] == sign)
              shape(ndof_h1+l) = 0.0;
          }
        }

        double fac = mip.GetWeight();
        elvec += (fac*coef) * shape;
      }
    }
  }
*/

  template class XMassIntegrator<2>;
  template class XMassIntegrator<3>;

  template class XSourceIntegrator<2>;
  template class XSourceIntegrator<3>;

  static RegisterBilinearFormIntegrator<XMassIntegrator<2> > initxh1cut2d ("xmass", 2, 2);
  static RegisterBilinearFormIntegrator<XMassIntegrator<3> > initxh1cut3d ("xmass", 3, 2);

  static RegisterLinearFormIntegrator<XSourceIntegrator<2> > initxh1source2d ("xsource", 2, 2);
  static RegisterLinearFormIntegrator<XSourceIntegrator<3> > initxh1source3d ("xsource", 3, 2);

  template class SpaceTimeXMassIntegrator<2>;
  template class SpaceTimeXMassIntegrator<3>;

  template class SpaceTimeXSourceIntegrator<2>;
  template class SpaceTimeXSourceIntegrator<3>;

  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<2> > initstxh1cut2d ("stxmass", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<3> > initstxh1cut3d ("stxmass", 3, 4);

  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<2> > initstxh1source2d ("stxsource", 2, 4);
  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<3> > initstxh1source3d ("stxsource", 3, 4);

/*
  static RegisterBilinearFormIntegrator<XLaplaceIntegrator<2> > initxlap2d ("xlaplace", 2, 2);
  static RegisterBilinearFormIntegrator<XLaplaceIntegrator<3> > initxlap3d ("xlaplace", 3, 2);

  static RegisterBilinearFormIntegrator<XConvectionIntegrator<2> > initxconv2d ("xconvection", 2, 2);
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<3> > initxconv3d ("xconvection", 3, 2);

  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HALFHALF> > initxnitsche2d_1 ("xnitsche_halfhalf", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_2 ("xnitsche_hansbo", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBO> > initxnitsche2d_2b ("xnitsche", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::HANSBOBETA> > initxnitsche2d_3 ("xnitsche_hansbobeta", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::BETA> > initxnitsche2d_4 ("xnitsche_beta", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHA> > initxnitsche2d_5 ("xnitsche_alpha", 2, 5);
  static RegisterBilinearFormIntegrator<XNitscheIntegrator<2,NITSCHE_VARIANTS::ALPHABETA> > initxnitsche2d_6 ("xnitsche_alphabeta", 2, 5);
  // static RegisterBilinearFormIntegrator<XNitscheIntegrator<3> > initxnitsche3d ("xnitsche", 3, 5);

  static RegisterBilinearFormIntegrator<XRobinIntegrator<2> > initxrob2d ("xrobin", 2, 2);
  static RegisterBilinearFormIntegrator<XRobinIntegrator<3> > initxrob3d ("xrobin", 3, 2);
*/
/*
  static RegisterLinearFormIntegrator<XNeumannIntegrator<2> > initxh1neumann2d ("xneumann", 2, 2);
  static RegisterLinearFormIntegrator<XNeumannIntegrator<3> > initxh1neumann3d ("xneumann", 3, 2);
*/

}

/// coefficientfunction statt function-pointer
