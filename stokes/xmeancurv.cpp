#include "xmeancurv.hpp"

namespace ngfem
{

  template<int D, bool improved>
  void XLBMeanCurvIntegrator<D,improved> ::
  CalcElementVector (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    static Timer timer ("XmodLBMeanCurvIntegrator::CalcElementMatrix");
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
    FlatVector<> shapeuvx(ndofuv,&shapeuv_total(ndofuv));
    
    // int orderuv = scafe.Order();

    // int ps = orderuv+1;

    // int p = scafe.Order();

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    FlatMatrixFixWidth<D*D> bmat(ndof_total,lh);
    FlatMatrixFixWidth<D*D> bmat_P(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape(ndofuv_total,lh);
    
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


      Vec<D> improved_normal = normal;
      // coef_lset num. diff
      if (improved)
      {
        Vec<D> dlset_ref = 0.0;
        double eps = 1e-7;
        for (int d = 0; d < D; ++d)
        {
          IntegrationPoint ipl(ip); 
          ipl(d) -= eps;
          MappedIntegrationPoint<D,D> mipl(ipl, eltrans);

          IntegrationPoint ipr(ip); ipr(d) += eps;
          MappedIntegrationPoint<D,D> mipr(ipr, eltrans);
          dlset_ref(d) = 0.5/eps * (coef_lset->Evaluate(mipr) - coef_lset->Evaluate(mipl));
        }
        improved_normal = Trans(Finv) * dlset_ref;
        improved_normal /= L2Norm(improved_normal);
      }

      // const double coef_val = coef->Evaluate(mip);

      scafe.CalcMappedDShape (mip,dshape);

      bmat = 0;

      // the first nd_u shape functions belong to u_x, the next nd_u belong to u_y:
      for (int d = 0; d < D; ++d)
        bmat.Rows(*dofrangeuv[d]).Cols(D*d,D*(d+1)) = dshape;
      if (!xfe->Empty())
      {
        throw Exception("not yet - how to deal with disc. vel?");
        // gradux = gradu;

        // for (int l = 0; l < ndofuv_x; ++l)
        //   {
        //     if (xfe->GetSignsOfDof()[l] != dt)
        //       gradux(l) = 0.0;
        //   }

        // for (int d = 0; d < D; ++d)
        //   bmat.Rows(D*d,D*(d+1)).Cols(*dofrangeuv_x[d]) = Trans (gradu);
      }
      
      // std::cout << " dshape = " << dshape << std::endl;
      // std::cout << " bmat = " << bmat << std::endl;

      Mat<D> P;
      P = Id<D>() - normal * Trans(normal);

      Mat<D> improved_P;
      improved_P = Id<D>() - improved_normal * Trans(improved_normal);
      
      Mat<D*D> project;
      project = 0.0;
      for (int d = 0; d < D; ++d)
        project.Rows(d*D,(d+1)*D).Cols(d*D,(d+1)*D) = P;
      
      FlatVector<> P_as_vec(D*D, &improved_P(0,0)) ;
      
      bmat_P = bmat * project;

      elvec -= (coef_val * weight) * (bmat_P * P_as_vec);

    }
  }

  template class XLBMeanCurvIntegrator<2,false>;
  template class XLBMeanCurvIntegrator<3,false>;

  static RegisterLinearFormIntegrator<XLBMeanCurvIntegrator<2,false> > initstxh1cut2d ("xLBmeancurv", 2, 1);
  static RegisterLinearFormIntegrator<XLBMeanCurvIntegrator<3,false> > initstxh1cut3d ("xLBmeancurv", 3, 1);

  template class XLBMeanCurvIntegrator<2,true>;
  template class XLBMeanCurvIntegrator<3,true>;

  static RegisterLinearFormIntegrator<XLBMeanCurvIntegrator<2,true> > initxmodLBmeancurv2d ("xmodLBmeancurv", 2, 2);
  static RegisterLinearFormIntegrator<XLBMeanCurvIntegrator<3,true> > initxmodLBmeancurv3d ("xmodLBmeancurv", 3, 2);



  /*
  //  Evaluate (u_x, u_y) 
  class DiffOpIdU : public DiffOp<DiffOpIdU>
  {
    // 2 components:
    // u1 u2

  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2 };
    enum { DIFFORDER = 0 };
  
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (bfel);
      const ScalarFiniteElement<2> & fel_u = 
        dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
    
      int nd_u = fel_u.GetNDof();

      FlatVector<> vecu(nd_u, lh);
      fel_u.CalcShape (mip.IP(), vecu);

      mat = 0;
      mat.Row(0).Range(cfel.GetRange(0)) = vecu;
      mat.Row(1).Range(cfel.GetRange(1)) = vecu;
    }
  };


  class StokesUIntegrator 
    : public T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
  {
  public:
    ///
    StokesUIntegrator (const Array<shared_ptr<CoefficientFunction>> & )
      :  T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
         (DiagDMat<2> (make_shared<ConstantCoefficientFunction>(1)))
    { ; }

    ///
    virtual string Name () const { return "Stokes IdU"; }
  };

  */



}
