#include "xfemIntegrators.hpp"

namespace ngfem
{

  template<int D>
  void CastXScalarFiniteElements (const FiniteElement & base_fel,
                                  const ScalarFiniteElement<D> * & scafe,
                                  const XFiniteElement * & xfe,
                                  const XDummyFE * & dummfe)               
  {
    auto cfel = dynamic_cast<const CompoundFiniteElement&> (base_fel);
    xfe = nullptr;
    dummfe = nullptr;
    scafe = nullptr;
    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==nullptr)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==nullptr)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==nullptr)
        scafe = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel[i]);
    }
    if (!scafe) 
      throw Exception(" no scalar fe found!");
    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");
  }
  
  
  template<int D>
  void XMassIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("XMassInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;
    
    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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
          elmat += mip.GetWeight() * coef * shape * Trans(shape);
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
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          elmat += mip.GetWeight() * coef * shape_total * Trans(shape_total);
        }
      }
    }
  }


  template<int D>
  void XConvectionIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XConvectionInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatMatrixFixWidth<D> dshape(ndof,lh);

    FlatVector<> dudwshape_total(ndof_total,lh);
    FlatVector<> dudwshape(ndof,&dudwshape_total(0));
    FlatVector<> dudwshapex(ndof_x,&dudwshape_total(ndof));

    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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

          scafe->CalcMappedDShape(mip, dshape);
          scafe->CalcShape(mip.IP(), shape);

          dudwshape = dshape * conv;

          double fac = mip.GetWeight();
          elmat += fac * shape * Trans(dudwshape);
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
          
          scafe->CalcMappedDShape(mip, dshape);
          scafe->CalcShape(mip.IP(), shape);

          shapex = shape;
          dudwshape = dshape * conv;
          dudwshapex = dudwshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
            {
              shapex(l) = 0.0;
              dudwshapex(l) = 0.0;
            }
          }

          double fac = mip.GetWeight();
          elmat += fac * shape_total * Trans(dudwshape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }


  template<int D>
  void XLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XLaplaceInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
    FlatMatrixFixWidth<D> dshapex(ndof_x,&dshape_total(ndof,0));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("XSourceInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elvec = 0.0;

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XRobinInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D-1> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("XNeumannInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D-1> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elvec = 0.0;

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
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

  template class XConvectionIntegrator<2>;
  template class XConvectionIntegrator<3>;

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

  static RegisterBilinearFormIntegrator<XConvectionIntegrator<2> > initxconv2d ("xconvection", 2, 2);
  static RegisterBilinearFormIntegrator<XConvectionIntegrator<3> > initxconv3d ("xconvection", 3, 2);

}

/// coefficientfunction statt function-pointer
