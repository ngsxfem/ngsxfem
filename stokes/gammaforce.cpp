#include "gammaforce.hpp"

namespace ngfem
{

  template<int D>
  void GammaForceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    static Timer timer ("GammaForceIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);


    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const CompoundFiniteElement & feuv_comp =
      dynamic_cast<const CompoundFiniteElement&> (cfel[0]);

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D>&> (feuv_comp[0]);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&feuv_comp[1]);
    const XDummyFE * dummy_fe =
      dynamic_cast<const XDummyFE *> (&feuv_comp[1]);

    // cout << " here a " << endl; getchar();

    elvec = 0.0;

    if (!xfe && !dummy_fe) 
      throw Exception(" not containing X-elements (velocity)?");

    if(!xfe)
      { 
        if (dummy_fe->GetDomainType() != IF)
          return;
        else
          throw Exception("dummy on the interface ?!?!?!");
      }

    int ndofuv_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndofuv = scafe.GetNDof();
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

    int ndof_total = base_fel.GetNDof();

    FlatVector<> shapeuv_total(ndofuv_total,lh);
    FlatVector<> shapeuv(ndofuv,&shapeuv_total(0));
    FlatVector<> shapeuvx(ndofuv_x,&shapeuv_total(ndofuv));
    
    // int orderuv = scafe.Order();

    // int ps = orderuv+1;

    // int p = scafe.Order();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    FlatMatrixFixWidth<D*D> bmat(ndof_total,lh);
    FlatMatrixFixWidth<D*D> bmat_P(ndof_total,lh);
    FlatMatrixFixWidth<D> vshape(ndofuv_total,lh);
    
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

      const double coef_val = coef->Evaluate(mip);

      scafe.CalcShape (ip,shapeuv_total);
      for (int i = 0; i < D; ++i)
	{
	  vshape.Col(i) = shapeuv_total;
	}      

      elvec += (coef_val * weight) * vshape * normal;

    }
  }

  template class GammaForceIntegrator<2>;
  template class GammaForceIntegrator<3>;

  static RegisterLinearFormIntegrator<GammaForceIntegrator<2>> initgammaforce2d ("xGammaForce", 2, 1);
  static RegisterLinearFormIntegrator<GammaForceIntegrator<3>> initgammaforce3d ("xGammaForce", 3, 1);

}
