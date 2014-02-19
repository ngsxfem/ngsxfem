#include "stxfemIntegrators.hpp"

namespace ngfem
{

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

          double fac = mip.GetMeasure() * quad.weights[i] * tau;
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D>
  void SpaceTimeXLaplaceIntegrator<D> ::
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
    FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
    FlatMatrixFixWidth<D> dshapex(ndof,&dshape_total(ndof,0));

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

            scafe->CalcMappedDxShapeSpaceTime(mip, irt[k](0), dshape, lh);

            double fac = mip.GetWeight() * irt[k].Weight() * tau;
            elmat += (fac*coef) * dshape * Trans(dshape);
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

          scafe->CalcMappedDxShapeSpaceTime(mip, quad.points[i](D), dshape, lh);
          dshapex = dshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              dshapex.Row(l) = 0.0;
          }

          double fac = mip.GetMeasure() * quad.weights[i] * tau;
          elmat += (fac*coef) * dshape_total * Trans(dshape_total);
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

          double fac = mip.GetMeasure() * quad.weights[i] * tau;
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

  template class SpaceTimeXMassIntegrator<2>;
  template class SpaceTimeXMassIntegrator<3>;

  template class SpaceTimeXSourceIntegrator<2>;
  template class SpaceTimeXSourceIntegrator<3>;

  template class SpaceTimeXLaplaceIntegrator<2>;
  template class SpaceTimeXLaplaceIntegrator<3>;

  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<2> > initstxh1cut2d ("stxmass", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<3> > initstxh1cut3d ("stxmass", 3, 4);

  static RegisterBilinearFormIntegrator<SpaceTimeXLaplaceIntegrator<2> > initstxh1cut2dlap ("stxlaplace", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXLaplaceIntegrator<3> > initstxh1cut3dlap ("stxlaplace", 3, 4);

  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<2> > initstxh1source2d ("stxsource", 2, 4);
  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<3> > initstxh1source3d ("stxsource", 3, 4);

/*
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<2> > initxconv2d ("xconvection", 2, 2);
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<3> > initxconv3d ("xconvection", 3, 2);

  static RegisterBilinearFormIntegrator<XRobinIntegrator<2> > initxrob2d ("xrobin", 2, 2);
  static RegisterBilinearFormIntegrator<XRobinIntegrator<3> > initxrob3d ("xrobin", 3, 2);
*/
/*
  static RegisterLinearFormIntegrator<XNeumannIntegrator<2> > initxh1neumann2d ("xneumann", 2, 2);
  static RegisterLinearFormIntegrator<XNeumannIntegrator<3> > initxh1neumann3d ("xneumann", 3, 2);
*/

}
