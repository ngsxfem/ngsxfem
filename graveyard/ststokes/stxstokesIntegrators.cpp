#include "stxstokesIntegrators.hpp"

namespace ngfem
{

  template<int D>
  void SpaceTimeXStokesIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("SpaceTimeXStokesIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const ScalarSpaceTimeFiniteElement<D> & scafeuv =
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D>&> (cfel[0]);

    const CompoundFiniteElement & cfelp = 
      dynamic_cast<const CompoundFiniteElement&> (cfel[2]);

    cout << " here " << endl; getchar();

    const XFiniteElement * xfe = NULL;
    const XDummyFE * dummfe = NULL;
    const ScalarSpaceTimeFiniteElement<D> * scafep = NULL;

    for (int i = 0; i < cfelp.GetNComponents(); ++i)
    {
      if (xfe==NULL)
        xfe = dynamic_cast<const XFiniteElement* >(&cfelp[i]);
      if (dummfe==NULL)
        dummfe = dynamic_cast<const XDummyFE* >(&cfelp[i]);
      if (scafep==NULL)
        scafep = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfelp[i]);
    }
    
    elmat = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    int ndofuv = scafeuv.GetNDof();

    int ndofp_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndofp = scafep->GetNDof();
    int ndofp_total = ndofp+ndofp_x;

    FlatVector<> shapep_total(ndofp_total,lh);
    FlatVector<> shapep(ndofp,&shapep_total(0));
    FlatVector<> shapepx(ndofp,&shapep_total(ndofp));
    
    FlatMatrixFixWidth<D> dshapeu(ndofuv,lh);

    int orderuv_s = scafeuv.OrderSpace();
    int orderuv_t = scafeuv.OrderTime();
    int orderp_s = scafep->OrderSpace();
    int orderp_t = scafep->OrderTime();

    int ps = max(orderuv_s-1,orderp_s)*2;
    int pt = max(orderuv_t,orderp_t)*2;

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

            scafep->CalcShapeSpaceTime(mip.IP(), irt[k](0), shapep, lh);
            double fac = mip.GetWeight() * irt[k].Weight() * tau;

            // elmat += (fac*coef) * shapep * Trans(shapep);
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

          scafep->CalcShapeSpaceTime(ips, fquad.points(i,D), shapep, lh);
          shapepx = shapep;

          for (int l = 0; l < ndofp_x; ++l)
          {
            if (xfe->GetSignsOfDof()[l] != dt)
              shapepx(l) = 0.0;
          }

          double fac = mip.GetMeasure() * fquad.weights(i) * tau;
          // elmat += (fac*coef) * shapep_total * Trans(shapep_total);
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template class SpaceTimeXStokesIntegrator<2>;
  template class SpaceTimeXStokesIntegrator<3>;

  static RegisterBilinearFormIntegrator<SpaceTimeXStokesIntegrator<2> > initstxh1cut2d ("stx_stokes", 2, 4);
  static RegisterBilinearFormIntegrator<SpaceTimeXStokesIntegrator<3> > initstxh1cut3d ("stx_stokes", 3, 4);

}
