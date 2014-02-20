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
  void SpaceTimeXTimeDerivativeIntegrator<D> ::
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

    FlatVector<> dtshape_total(ndof_total,lh);
    FlatVector<> dtshape(ndof,&dtshape_total(0));
    FlatVector<> dtshapex(ndof,&dtshape_total(ndof));

    int ps = scafe->OrderSpace();
    int pt = scafe->OrderTime();
    
    // const double t1 = coef_tnew->EvaluateConst(); 
    // const double t0 = coef_told->EvaluateConst();
    // const double tau = t1 - t0;

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
            scafe->CalcDtShapeSpaceTime(mip.IP(), irt[k](0), dtshape, lh);
            double fac = mip.GetWeight() * irt[k].Weight();
            elmat += (fac*coef) * shape * Trans(dtshape);
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
          scafe->CalcDtShapeSpaceTime(ips, quad.points[i](D), dtshape, lh);

          shapex = shape;
          dtshapex = dtshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
            {
              shapex(l) = 0.0;
              dtshapex(l) = 0.0;
            }
          }
          double fac = mip.GetMeasure() * quad.weights[i];
          
          elmat += (fac*coef) * shape_total * Trans(dtshape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D>
  void SpaceTimeXConvectionIntegrator<D> ::
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

    FlatVector<> bgradshape_total(ndof_total,lh);
    FlatVector<> bgradshape(ndof,&bgradshape_total(0));
    FlatVector<> bgradshapex(ndof,&bgradshape_total(ndof));

    FlatMatrixFixWidth<D> gradshape(ndof,lh);

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
            Vec<D> conv;
            if (dt == POS)
                coef_pos->Evaluate(mip,conv);
            else
                coef_neg->Evaluate(mip,conv);

            scafe->CalcShapeSpaceTime(mip.IP(), irt[k](0), shape, lh);
            scafe->CalcMappedDxShapeSpaceTime(mip, irt[k](0), gradshape, lh);
            bgradshape = gradshape * conv;

            double fac = mip.GetWeight() * irt[k].Weight() * tau;
            elmat += fac * shape * Trans(bgradshape);
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

          Vec<D> conv;
          if (dt == POS)
              coef_pos->Evaluate(mip,conv);
          else
              coef_neg->Evaluate(mip,conv);

          scafe->CalcShapeSpaceTime(ips, quad.points[i](D), shape, lh);
          scafe->CalcMappedDxShapeSpaceTime(mip, quad.points[i](D), gradshape, lh);
          bgradshape = gradshape * conv;

          shapex = shape;
          bgradshapex = bgradshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
            {
              shapex(l) = 0.0;
              bgradshapex(l) = 0.0;
            }
          }

          double fac = mip.GetMeasure() * quad.weights[i] * tau;
          elmat += fac * shape_total * Trans(bgradshape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D, TIME t>
  void SpaceTimeXTraceMassIntegrator<D,t> ::
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
    
    const double tracetime = t == PAST ? 0.0 : 1.0;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule irs = SelectIntegrationRule (eltrans.GetElementType(), 2*ps);
        for (int l = 0 ; l < irs.GetNIP(); l++)
        {
            MappedIntegrationPoint<D,D> mip(irs[l], eltrans);
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
                : coef_neg->Evaluate(mip);
            scafe->CalcShapeSpaceTime(mip.IP(), tracetime, shape, lh);
            double fac = mip.GetWeight();
            elmat += (fac*coef) * shape * Trans(shape);
        }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_st_p = xfe->GetLocalGeometry();

        const XLocalGeometryInformation * lset_eval_p = t == PAST ? 
            xfe->GetLocalGeometry()->GetPastTrace() : xfe->GetLocalGeometry()->GetFutureTrace();

        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D> * compr (lset_eval_p->GetCompositeRule<D>());
        const QuadratureRule<D> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ips(quad.points[i]);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, tracetime, shape, lh);
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

          double fac = mip.GetMeasure() * quad.weights[i] * tau;
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D, TIME t>
  void SpaceTimeXTraceSourceIntegrator<D,t> ::
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
    
    const double tracetime = t == PAST ? 0.0 : 1.0;

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule irs = SelectIntegrationRule (eltrans.GetElementType(), 2*ps);
        for (int l = 0 ; l < irs.GetNIP(); l++)
        {
            MappedIntegrationPoint<D,D> mip(irs[l], eltrans);
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
                : coef_neg->Evaluate(mip);
            scafe->CalcShapeSpaceTime(mip.IP(), tracetime, shape, lh);
            double fac = mip.GetWeight();
            elvec += (fac*coef) * shape;
        }
      }
      else
      {
        const XLocalGeometryInformation * lset_eval_p = t == PAST ? 
            xfe->GetLocalGeometry()->GetPastTrace() : xfe->GetLocalGeometry()->GetFutureTrace();
        if (lset_eval_p == NULL)
          throw Exception(" no local geometry");
        const CompositeQuadratureRule<D> * compr (lset_eval_p->GetCompositeRule<D>());
        const QuadratureRule<D> & quad(compr->GetRule(dt));
        for (int i = 0; i < quad.Size(); ++i)
        {
          IntegrationPoint ips(quad.points[i]);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, tracetime, shape, lh);
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

  template class SpaceTimeXLaplaceIntegrator<2>;
  template class SpaceTimeXLaplaceIntegrator<3>;

  template class SpaceTimeXConvectionIntegrator<2>;
  template class SpaceTimeXConvectionIntegrator<3>;

  template class SpaceTimeXTimeDerivativeIntegrator<2>;
  template class SpaceTimeXTimeDerivativeIntegrator<3>;

  template class SpaceTimeXTraceMassIntegrator<2,PAST>;
  template class SpaceTimeXTraceMassIntegrator<3,PAST>;
  template class SpaceTimeXTraceMassIntegrator<2,FUTURE>;
  template class SpaceTimeXTraceMassIntegrator<3,FUTURE>;

  template class SpaceTimeXSourceIntegrator<2>;
  template class SpaceTimeXSourceIntegrator<3>;

  template class SpaceTimeXTraceSourceIntegrator<2,PAST>;
  template class SpaceTimeXTraceSourceIntegrator<3,PAST>;
  template class SpaceTimeXTraceSourceIntegrator<2,FUTURE>;
  template class SpaceTimeXTraceSourceIntegrator<3,FUTURE>;

  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<2> > initstxh1cut2d ("stx_mass", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXMassIntegrator<3> > initstxh1cut3d ("stx_mass", 3, 4);

  static RegisterBilinearFormIntegrator<SpaceTimeXLaplaceIntegrator<2> > initstxh1cut2dlap ("stx_laplace", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXLaplaceIntegrator<3> > initstxh1cut3dlap ("stx_laplace", 3, 4);

  static RegisterBilinearFormIntegrator<SpaceTimeXConvectionIntegrator<2> > initstxh1cut2dconv ("stx_convection", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXConvectionIntegrator<3> > initstxh1cut3dconv ("stx_convection", 3, 4);

  static RegisterBilinearFormIntegrator<SpaceTimeXTimeDerivativeIntegrator<2> > initstxh1cut2dtimeder ("stx_timeder", 2, 2);
  static RegisterBilinearFormIntegrator<SpaceTimeXTimeDerivativeIntegrator<3> > initstxh1cut3dtimeder ("stx_timeder", 3, 2);

  static RegisterBilinearFormIntegrator<SpaceTimeXTraceMassIntegrator<2,PAST> > initstxh1cut2dtracmp ("stx_tracemass_past", 2, 2);
  static RegisterBilinearFormIntegrator<SpaceTimeXTraceMassIntegrator<3,PAST> > initstxh1cut3dtracmp ("stx_tracemass_past", 3, 2);
  static RegisterBilinearFormIntegrator<SpaceTimeXTraceMassIntegrator<2,FUTURE> > initstxh1cut2dtracmf ("stx_tracemass_future", 2, 2);
  static RegisterBilinearFormIntegrator<SpaceTimeXTraceMassIntegrator<3,FUTURE> > initstxh1cut3dtracmf ("stx_tracemass_future", 3, 2);

  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<2> > initstxh1source2d ("stx_source", 2, 4);
  static RegisterLinearFormIntegrator<SpaceTimeXSourceIntegrator<3> > initstxh1source3d ("stx_source", 3, 4);

  static RegisterLinearFormIntegrator<SpaceTimeXTraceSourceIntegrator<2,PAST> > initstxh1cut2dtracsp ("stx_tracesource_past", 2, 2);
  static RegisterLinearFormIntegrator<SpaceTimeXTraceSourceIntegrator<3,PAST> > initstxh1cut3dtracsp ("stx_tracesource_past", 3, 2);
  static RegisterLinearFormIntegrator<SpaceTimeXTraceSourceIntegrator<2,FUTURE> > initstxh1cut2dtracsf ("stx_tracesource_future", 2, 2);
  static RegisterLinearFormIntegrator<SpaceTimeXTraceSourceIntegrator<3,FUTURE> > initstxh1cut3dtracsf ("stx_tracesource_future", 3, 2);

/*
  static RegisterBilinearFormIntegrator<XRobinIntegrator<2> > initxrob2d ("xrobin", 2, 2);
  static RegisterBilinearFormIntegrator<XRobinIntegrator<3> > initxrob3d ("xrobin", 3, 2);

  static RegisterLinearFormIntegrator<XNeumannIntegrator<2> > initxh1neumann2d ("xneumann", 2, 2);
  static RegisterLinearFormIntegrator<XNeumannIntegrator<3> > initxh1neumann3d ("xneumann", 3, 2);
*/

}
