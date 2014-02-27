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
    static Timer timer ("SpaceTimeXMassIntegrator::CalcElementMatrix");
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
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
    static Timer timer ("SpaceTimeXLaplaceIntegrator::CalcElementMatrix");
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcMappedDxShapeSpaceTime(mip, fquad.points(i,D), dshape, lh);
          dshapex = dshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              dshapex.Row(l) = 0.0;
          }

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
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
    static Timer timer ("SpaceTimeXTimeDerivativeIntegrator::CalcElementMatrix");
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

        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D), shape, lh);
          scafe->CalcDtShapeSpaceTime(ips, fquad.points(i,D), dtshape, lh);

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
          double fac = mip.GetMeasure() * fquad.weights(i);
          
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
    static Timer timer ("SpaceTimeXConvectionIntegrator::CalcElementMatrix");
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);

          Vec<D> conv;
          if (dt == POS)
              coef_pos->Evaluate(mip,conv);
          else
              coef_neg->Evaluate(mip,conv);

          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D), shape, lh);
          scafe->CalcMappedDxShapeSpaceTime(mip, fquad.points(i,D), gradshape, lh);
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

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
          elmat += fac * shape_total * Trans(bgradshape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

/*
  template<int D>
  void SpaceTimeXTransportIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXLaplaceIntegrator::CalcElementMatrix");
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

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;

    FlatMatrixFixWidth<D+2> bmat_total(ndof_total,lh);
    FlatMatrixFixWidth<D+2> bmat(ndof,&bmat_total(0,0));
    FlatMatrixFixWidth<D+2> bmatx(ndof,&bmat_total(ndof,0));

    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof,&shape_total(ndof));

    FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
    FlatMatrixFixWidth<D> dshapex(ndof,&dshape_total(ndof,0));

    FlatVector<> dtshape_total(ndof_total,lh);
    FlatVector<> dtshape(ndof,&dtshape_total(0));
    FlatVector<> dtshapex(ndof,&dtshape_total(ndof));

    int ps = scafe->OrderSpace();
    int pt = scafe->OrderTime();
    
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
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcMappedDxShapeSpaceTime(mip, fquad.points(i,D), dshape, lh);
          dshapex = dshape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              dshapex.Row(l) = 0.0;
          }

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
          elmat += (fac*coef) * dshape_total * Trans(dshape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

*/

  template<int D, TIME t>
  void SpaceTimeXTraceMassIntegrator<D,t> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXTraceMassIntegrator::CalcElementMatrix");
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

        const FlatXLocalGeometryInformation & xgeom( t == PAST ?
                                                     xfe->GetFlatLocalGeometryDownTrace() :
                                                     xfe->GetFlatLocalGeometryUpTrace() );
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, tracetime, shape, lh);
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
  void SpaceTimeXSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXSourceIntegrator::CalcElementVector");
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
            MappedIntegrationPoint<D,D> mips(irs[l], eltrans);
            DimMappedIntegrationPoint<D+1> mip(irs[l],eltrans);
            mip.Point().Range(0,D) = mips.GetPoint();
            mip.Point()[D] = t0 + irt[k](0) * tau;
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
              : coef_neg->Evaluate(mip);

            scafe->CalcShapeSpaceTime(mips.IP(), irt[k](0), shape, lh);
            double fac = mips.GetWeight() * irt[k].Weight() * tau;
            elvec += (fac*coef) * shape;
          }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D+1> & fcompr(xgeom.GetCompositeRule<D+1>());
        const FlatQuadratureRule<D+1> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips;
          for (int d = 0; d < D; ++d)
            ips(d) = fquad.points(i,d);

          MappedIntegrationPoint<D,D> mips(ips, eltrans);
          DimMappedIntegrationPoint<D+1> mip(ips,eltrans);
          mip.Point().Range(0,D) = mips.GetPoint();
          mip.Point()[D] = t0 + fquad.points(i,D) * tau;
          double coef = dt == POS ? coef_pos->Evaluate(mip) 
            : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mips.GetMeasure() * fquad.weights(i) * tau;
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
    static Timer timer ("SpaceTimeXTraceSourceIntegrator::CalcElementVector");
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
            double coef = dt == POS ? coef_pos->Evaluate(mip) * scale_pos
                : coef_neg->Evaluate(mip) * scale_neg;
            scafe->CalcShapeSpaceTime(mip.IP(), tracetime, shape, lh);
            double fac = mip.GetWeight();
            elvec += (fac*coef) * shape;
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom( t == PAST ?
                                                     xfe->GetFlatLocalGeometryDownTrace() :
                                                     xfe->GetFlatLocalGeometryUpTrace() );
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) * scale_pos : coef_neg->Evaluate(mip) * scale_neg;
          scafe->CalcShapeSpaceTime(ips, tracetime, shape, lh);
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
  void SpaceTimeXRobinIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXRobinIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarSpaceTimeFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarSpaceTimeFiniteElement<D-1>* >(&cfel[i]);
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
            MappedIntegrationPoint<D-1,D> mip(irs[l], eltrans);
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
              : coef_neg->Evaluate(mip);
            scafe->CalcShapeSpaceTime(mip.IP(), irt[k](0), shape, lh);
            double fac = mip.GetWeight() * irt[k].Weight() * tau;
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
          IntegrationPoint ips;
          for (int d = 0; d < D-1; ++d)
            ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D-1,D> mip(ips, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D-1), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
          elmat += (fac*coef) * shape_total * Trans(shape_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }



  template<int D>
  void SpaceTimeXNeumannIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXNeumannIntegrator::CalcElementVector");
    RegionTimer reg (timer);
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarSpaceTimeFiniteElement<D-1> * scafe = NULL;

    for (int i = 0; i < cfel.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfel[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfel[i]);
      if (scafe==NULL)
        scafe = dynamic_cast<const ScalarSpaceTimeFiniteElement<D-1>* >(&cfel[i]);
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
            MappedIntegrationPoint<D-1,D> mips(irs[l], eltrans);
            DimMappedIntegrationPoint<D+1> mip(irs[l],eltrans);
            mip.Point().Range(0,D) = mips.GetPoint();
            mip.Point()[D] = t0 + irt[k](0) * tau;
            double coef = dt == POS ? coef_pos->Evaluate(mip) 
              : coef_neg->Evaluate(mip);

            scafe->CalcShapeSpaceTime(mips.IP(), irt[k](0), shape, lh);
            double fac = mips.GetWeight() * irt[k].Weight() * tau;
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
          IntegrationPoint ips;
          for (int d = 0; d < D-1; ++d)
            ips(d) = fquad.points(i,d);

          MappedIntegrationPoint<D-1,D> mips(ips, eltrans);
          DimMappedIntegrationPoint<D+1> mip(ips,eltrans);
          mip.Point().Range(0,D) = mips.GetPoint();
          mip.Point()[D] = t0 + fquad.points(i,D-1) * tau;
          double coef = dt == POS ? coef_pos->Evaluate(mip) 
            : coef_neg->Evaluate(mip);

          scafe->CalcShapeSpaceTime(ips, fquad.points(i,D-1), shape, lh);
          shapex = shape;

          for (int l = 0; l < ndof_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          }

          double fac = mips.GetMeasure() * fquad.weights(i) * tau;
          elvec += (fac*coef) * shape_total;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }




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

  template class SpaceTimeXRobinIntegrator<2>;
  template class SpaceTimeXRobinIntegrator<3>;

  template class SpaceTimeXNeumannIntegrator<2>;
  template class SpaceTimeXNeumannIntegrator<3>;

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

  static RegisterBilinearFormIntegrator<SpaceTimeXRobinIntegrator<2> > initstxh1cut_rob_2d ("stx_robin", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXRobinIntegrator<3> > initstxh1cut_rob_3d ("stx_robin", 3, 4);

  static RegisterLinearFormIntegrator<SpaceTimeXNeumannIntegrator<2> > initstxh1neumann2d ("stx_neumann", 2, 4);
  static RegisterLinearFormIntegrator<SpaceTimeXNeumannIntegrator<3> > initstxh1neumann3d ("stx_neumann", 3, 4);

}
