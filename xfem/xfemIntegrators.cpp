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
    static Timer timer ("XMassIntegrator::CalcElementMatrix");
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetWeight();
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D>
  void XLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XLaplaceIntegrator::CalcElementMatrix");
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
    FlatMatrixFixWidth<D> dshapex(ndof,&dshape_total(ndof,0));

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
          double coef = dt == POS ? alpha_pos : alpha_neg;
          scafe->CalcMappedDShape(mip, dshape);
          double fac = mip.GetWeight();
          elmat += (fac*coef) * dshape * Trans(dshape);
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
          double coef = dt == POS ? alpha_pos : alpha_neg;
          
          scafe->CalcMappedDShape(mip, dshape);
          dshapex = dshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              dshapex.Row(l) = 0.0;
          }

          double fac = mip.GetWeight();
          elmat += (fac*coef) * dshape_total * Trans(dshape_total);
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
    static Timer timer ("XSourceIntegrator::CalcElementMatrix");
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetWeight();
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }



  template<int D>
  void XRobinIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XRobinIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarFiniteElement<D-1>* >(&cfel[i]);
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
          MappedIntegrationPoint<D-1,D> mip(pir[i], eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          shape.Range(0,ndof) = scafe->GetShape(mip.IP(), lh);
          double fac = mip.GetWeight();
          elmat += (fac*coef) * shape * Trans(shape);
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D-1> & fcompr(xgeom.GetCompositeRule<D-1>());
        const FlatQuadratureRule<D-1> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D-1,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetWeight();
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template<int D>
  void XNeumannIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("XNeumannIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarFiniteElement<D-1>* >(&cfel[i]);
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
          MappedIntegrationPoint<D-1,D> mip(pir[i], eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          shape = scafe->GetShape(mip.IP(), lh);
          double fac = mip.GetWeight();
          elvec += (fac*coef) * shape;
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D-1> & fcompr(xgeom.GetCompositeRule<D-1>());
        const FlatQuadratureRule<D-1> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D-1,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetWeight();
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template class XMassIntegrator<2>;
  template class XMassIntegrator<3>;

  template class XLaplaceIntegrator<2>;
  template class XLaplaceIntegrator<3>;

  template class XSourceIntegrator<2>;
  template class XSourceIntegrator<3>;

  template class XRobinIntegrator<2>;
  template class XRobinIntegrator<3>;

  template class XNeumannIntegrator<2>;
  template class XNeumannIntegrator<3>;


  static RegisterBilinearFormIntegrator<XMassIntegrator<2> > initxh1cut2d ("xmass", 2, 2);
  static RegisterBilinearFormIntegrator<XMassIntegrator<3> > initxh1cut3d ("xmass", 3, 2);

  static RegisterBilinearFormIntegrator<XLaplaceIntegrator<2> > initxh1cut2dlap ("xlaplace", 2, 2);
  static RegisterBilinearFormIntegrator<XLaplaceIntegrator<3> > initxh1cut3dlap ("xlaplace", 3, 2);

  static RegisterLinearFormIntegrator<XSourceIntegrator<2> > initxh1source2d ("xsource", 2, 2);
  static RegisterLinearFormIntegrator<XSourceIntegrator<3> > initxh1source3d ("xsource", 3, 2);

  static RegisterBilinearFormIntegrator<XRobinIntegrator<2> > initxrob2d ("xrobin", 2, 2);
  static RegisterBilinearFormIntegrator<XRobinIntegrator<3> > initxrob3d ("xrobin", 3, 2);

  static RegisterLinearFormIntegrator<XNeumannIntegrator<2> > initxh1neumann2d ("xneumann", 2, 2);
  static RegisterLinearFormIntegrator<XNeumannIntegrator<3> > initxh1neumann3d ("xneumann", 3, 2);

/*
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<2> > initxconv2d ("xconvection", 2, 2);
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<3> > initxconv3d ("xconvection", 3, 2);
*/

  template<int D>
  void NoXLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("NoXLaplaceIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> & scafe = 
      dynamic_cast<const ScalarFiniteElement<D>&> (base_fel);
    
    int order=scafe.Order();
    elmat = 0.0;
    
    ELEMENT_TYPE eltype = eltrans.GetElementType();

    ScalarFieldEvaluator & lset_eval = * ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,lh);

    CompositeQuadratureRule<D> * cquad = new CompositeQuadratureRule<D>() ;
    XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT, lset_eval, 
                                                                          *cquad, lh, 
                                                                          2*order, 2, 
                                                                          0,0);
    DOMAIN_TYPE dt_el = xgeom->MakeQuadRule();
    FlatXLocalGeometryInformation fxgeom(*xgeom,lh);
      
    int ndof = scafe.GetNDof();
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      if(dt_el == NEG || dt_el == POS)
      { 
        if (dt_el != dt)
          continue;
        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*order);
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          double coef = dt == POS ? alpha_pos : alpha_neg;
          scafe.CalcMappedDShape(mip, dshape);
          double fac = mip.GetWeight();
          elmat += (fac*coef) * dshape * Trans(dshape);
        }
      }
      else
      {
        QuadratureRule<D> dtquad(cquad->GetRule(dt));
        for (int i = 0; i < dtquad.Size(); ++i)
        {
          IntegrationPoint ip(&dtquad.points[i](0),dtquad.weights[i]);
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? alpha_pos : alpha_neg;
          scafe.CalcMappedDShape(mip, dshape);
          double fac = mip.GetWeight();
          elmat += (fac*coef) * dshape * Trans(dshape);
        } // quad rule
      } // if xfe
    } // loop over domain types
    delete xgeom;

  }


  template class NoXLaplaceIntegrator<2>;
  template class NoXLaplaceIntegrator<3>;

  static RegisterBilinearFormIntegrator<NoXLaplaceIntegrator<2> > initxh1cut2dnolap ("noxlaplace", 2, 3);
  static RegisterBilinearFormIntegrator<NoXLaplaceIntegrator<3> > initxh1cut3dnolap ("noxlaplace", 3, 3);


}

/// coefficientfunction statt function-pointer
