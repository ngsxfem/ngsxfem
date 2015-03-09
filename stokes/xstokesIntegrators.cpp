#include "xstokesIntegrators.hpp"

namespace ngfem
{

  template<int D>
  void XStokesIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("XStokesIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);


    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const CompoundFiniteElement & feuv_comp =
      dynamic_cast<const CompoundFiniteElement&> (cfel[0]);

    const CompoundFiniteElement & fep_comp =
      dynamic_cast<const CompoundFiniteElement&> (cfel[D]);

    const ScalarFiniteElement<D> & feuv =
      dynamic_cast<const ScalarFiniteElement<D>&> (feuv_comp[0]);

    const XFiniteElement * feuvx =
      dynamic_cast<const XFiniteElement *> (&feuv_comp[1]);
    const XDummyFE * dummy_feuvx =
      dynamic_cast<const XDummyFE *> (&feuv_comp[1]);

    const ScalarFiniteElement<D> & fep = 
      dynamic_cast<const ScalarFiniteElement<D>&> (fep_comp[0]);

    const XFiniteElement * fepx =
      dynamic_cast<const XFiniteElement *> (&fep_comp[1]);
    const XDummyFE * dummy_fepx =
      dynamic_cast<const XDummyFE *> (&fep_comp[1]);

    elmat = 0.0;

    if (!feuvx && !dummy_feuvx) 
      throw Exception(" not containing X-elements (velocity)?");

    if (!fepx && !dummy_fepx) 
      throw Exception(" not containing X-elements (pressure)?");

    int ndofuv_x = feuvx!=NULL ? feuvx->GetNDof() : 0;
    int ndofuv = feuv.GetNDof();
    int ndofuv_total = ndofuv+ndofuv_x;

    shared_ptr<IntRange> dofrangeuv[D];
    shared_ptr<IntRange> dofrangeuv_x[D];
    int cnt = 0;
    for (int d = 0; d < D; ++d)
    {
      dofrangeuv[d] = make_shared<IntRange>(cnt,cnt+ndofuv);
      cnt += ndofuv;
      dofrangeuv_x[d] = make_shared<IntRange>(cnt,cnt+ndofuv_x);
      cnt += ndofuv_x;
    }

    int ndofp = fep.GetNDof();
    int ndofp_x = fepx!=NULL ? fepx->GetNDof() : 0;
    int ndofp_total = ndofp+ndofp_x;

    IntRange dofrangep(cnt, cnt+ndofp);
    cnt+= ndofp;
    IntRange dofrangep_x(cnt, cnt+ndofp_x);

    int ndof_total = ndofp_total + D*ndofuv_total;

    FlatVector<> shapep_total(ndofp_total,lh);
    FlatVector<> shapep(ndofp,&shapep_total(0));
    FlatVector<> shapepx(ndofp,&shapep_total(ndofp));

    FlatVector<> shapeuv_total(ndofuv_total,lh);
    FlatVector<> shapeuv(ndofuv,&shapeuv_total(0));
    FlatVector<> shapeuvx(ndofuv,&shapeuv_total(ndofuv));
    
    int orderuv = feuv.Order();
    int orderp = fep.Order();
    int ps = max(orderuv-1,orderp);

    Mat<D*D+1> mat = 0;

    FlatMatrixFixHeight<D*D+1> bmat(ndof_total,lh);

    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      // (nu  grud u, grad v)
      double val = (dt==POS) ? coef_pos -> EvaluateConst() : coef_neg -> EvaluateConst();
      for (int i = 0; i < D*D; i++)
        mat(i, i) = val;

      // (du1/dx+du2/dy, p)
      for (int d = 0; d < D; ++d)
        mat(D*D,(D+1)*d) = mat((D+1)*d,D*D) = -1;


      HeapReset hr(lh);

      if(!fepx)
      { 
        if (dummy_fepx->GetDomainType() != dt)
          continue;
        IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*ps);
        for (int l = 0 ; l < ir.GetNIP(); l++)
        {
          MappedIntegrationPoint<D,D> mip(ir[l], eltrans);

          FlatMatrixFixWidth<D> gradu(ndofuv, lh);
          feuv.CalcMappedDShape (mip, gradu);

          // the shape functions of the pressure
          FlatVector<> vecp(ndofp, lh);
          fep.CalcShape (mip.IP(), vecp);

          bmat = 0;

          // \nabla u
          // the first nd_u shape functions belong to u_x, the next nd_u belong to u_y:
          for (int d = 0; d < D; ++d)
            bmat.Rows(D*d,D*(d+1)).Cols(*dofrangeuv[d]) = Trans (gradu);

          // \nabla u^T
          // the first nd_u shape functions belong to u_x, the next nd_u belong to u_y:
          for (int d = 0; d < D; ++d)
            for (int k = 0; k < D; ++k)
              bmat.Row(d*D+k).Range(*dofrangeuv[k]) = gradu.Col(d);

          // ... and finally nd_p shape functions for the pressure:
          bmat.Row(D*D).Range(dofrangep) = vecp;

          elmat += mip.GetWeight() * Trans(bmat) * mat * bmat;
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(feuvx->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ips(&fquad.points(i,0),fquad.weights(i));

          // for (int d = 0; d < D; ++d)
          //   ips(d) = fquad.points(i,d);
          MappedIntegrationPoint<D,D> mip(ips, eltrans);

          FlatMatrixFixWidth<D> gradu(ndofuv, lh);
          feuv.CalcMappedDShape (mip, gradu);

          // the shape functions of the pressure
          FlatVector<> vecp(ndofp, lh);
          fep.CalcShape (mip.IP(), vecp);

          FlatMatrixFixWidth<D> gradux(ndofuv, lh);
          FlatVector<> vecpx(ndofp, lh);

          vecpx = vecp;
          gradux = gradu;

          for (int l = 0; l < ndofp_x; ++l)
          {
            if (fepx->GetSignsOfDof()[l] != dt)
              vecpx(l) = 0.0;
          }

          for (int l = 0; l < ndofuv_x; ++l)
          {
            if (feuvx->GetSignsOfDof()[l] != dt)
              gradux(l) = 0.0;
          }

          bmat = 0.0;
          // the first nd_u shape functions belong to u_x, the next nd_u belong to u_y:
          for (int d = 0; d < D; ++d)
            bmat.Rows(D*d,D*(d+1)).Cols(*dofrangeuv[d]) = Trans (gradu);

          // \nabla u^T
          for (int d = 0; d < D; ++d)
            for (int k = 0; k < D; ++k)
              bmat.Row(d*D+k).Range(*dofrangeuv[k]) = gradu.Col(d); // du_k / dx_d


          for (int d = 0; d < D; ++d)
            bmat.Rows(D*d,D*(d+1)).Cols(*dofrangeuv_x[d]) = Trans (gradux);

          // \nabla u_x^T
          for (int d = 0; d < D; ++d)
            for (int k = 0; k < D; ++k)
              bmat.Row(d*D+k).Range(*dofrangeuv_x[k]) = gradux.Col(d);

          // ... and finally nd_p shape functions for the pressure:
          bmat.Row(D*D).Range(dofrangep) = vecp;
          bmat.Row(D*D).Range(dofrangep_x) = vecpx;

          elmat += mip.GetWeight() * Trans(bmat) * mat * bmat;
        } // quad rule
      } // if xfe
    } // loop over domain types
  }

  template class XStokesIntegrator<2>;
  template class XStokesIntegrator<3>;

  static RegisterBilinearFormIntegrator<XStokesIntegrator<2> > initstxh1cut2d ("xstokes", 2, 2);
  static RegisterBilinearFormIntegrator<XStokesIntegrator<3> > initstxh1cut3d ("xstokes", 3, 2);

}
