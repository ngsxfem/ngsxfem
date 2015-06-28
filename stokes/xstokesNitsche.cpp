#include "xstokesNitsche.hpp"

namespace ngfem
{

  template <int D>
  void XStokesNitscheIntegrator<D>::
  CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    static Timer timer ("XStokesNitscheIntegrator::CalcElementMatrix");
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

    if(!feuvx)
      return;

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
    int ps = orderuv;

    Mat<D> mat = 0;

    FlatMatrixFixWidth<D> bmat(ndof_total,lh);
    FlatMatrixFixWidth<D> bmatjump(ndof_total,lh);

    bmat = 0.0;
    bmatjump = 0.0;

    FlatMatrix<> Nc(ndof_total,ndof_total,lh);
    FlatMatrix<> Ns(ndof_total,ndof_total,lh);
    
    Ns = 0.0;
    Nc = 0.0;

    const FlatArray<DOMAIN_TYPE>& xusign = feuvx->GetSignsOfDof();
    const FlatArray<DOMAIN_TYPE>& xpsign = fepx->GetSignsOfDof();
    
    //int p = scafe->Order();

    const FlatXLocalGeometryInformation & xgeom(feuvx->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());
    
    double kappa_neg;
    double kappa_pos;

    if (xgeom.kappa[NEG] >= 0.5)
      kappa_neg = 1.0;
    else
      kappa_neg = 0.0;
    kappa_pos = 1.0 - kappa_neg;

    const double lam = lambda->EvaluateConst();

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
      
      const double a_neg = alpha_neg->Evaluate(mip);
      const double a_pos = alpha_pos->Evaluate(mip);

      const double ava = a_pos*kappa_pos+a_neg*kappa_neg;
      
      shapep = fep.GetShape(mip.IP(), lh);
      shapep*=-1.0;
      
      FlatMatrixFixWidth<D> gradu(ndofuv, lh);
      feuv.CalcMappedDShape (mip, gradu);
      FlatVector<> dudn(ndofuv,lh);
      dudn = gradu * normal;

      bmat = 0.0;
      for (int d = 0; d<D; ++d)
      {
        bmat.Rows(*dofrangeuv[d]).Col(d)=(ava)*dudn;
        bmat.Rows(*dofrangeuv_x[d]).Col(d)=dudn;
        bmat.Rows(dofrangep).Col(d)=normal[d]*shapep;
        bmat.Rows(dofrangep_x).Col(d)=normal[d]*shapep;
        for (int k = 0; k < D; ++k)
        {
          bmat.Rows(*dofrangeuv[k]).Col(d) += (ava) * normal(k) * gradu.Col(d);
          bmat.Rows(*dofrangeuv_x[k]).Col(d) += normal(k) * gradu.Col(d);
        }
      }    
      
      shapeuv = feuv.GetShape(mip.IP(),lh);
      
      for (int d = 0; d<D; ++d)
      {
	  
        bmatjump.Rows(*dofrangeuv[d]).Col(d)=0.0;
        bmatjump.Rows(*dofrangeuv_x[d]).Col(d)=shapeuv;
        bmatjump.Rows(dofrangep).Col(d)=0.0;
        bmatjump.Rows(dofrangep_x).Col(d)=0.0;
      }  
      
      for (int d=0; d<D; ++d)
      {
        for (int l = 0; l < ndofuv_x; ++l)
        {
          if (xusign[l] == NEG){
            bmatjump.Row(dofrangeuv[d]->Next()+l)[d] *= -1.0;
            bmat.Row(dofrangeuv[d]->Next()+l) *= kappa_neg * a_neg; 
          }
          else{
            bmatjump.Row(dofrangeuv[d]->Next()+l)[d] *= 1.0;
            bmat.Row(dofrangeuv[d]->Next()+l) *= kappa_pos * a_pos; 
          }
        }
      }
      
      for (int l = 0; l < ndofp_x; ++l)
      {
        if (xpsign[l] == NEG){
          //bmatjump.Row(dofrangep.Next()+l) *= -1.0;
          bmat.Row(dofrangep.Next()+l) *= kappa_neg; 
        }
        else{
          //bmatjump.Row(dofrangep.Next()+l) *= 1.0;
          bmat.Row(dofrangep.Next()+l) *= kappa_pos; 
        }
      }
      
      /*
        cout<<"bmat = "<<bmat<<endl;
        cout<<"bmatjump = "<<bmatjump<<endl;
        getchar();
      */

      Nc = -weight * bmatjump * Trans(bmat);
      Ns = weight * bmatjump * Trans(bmatjump);

      elmat+= Nc + Trans(Nc) + ava*lam*(ps+1)*ps/h * Ns; 
      //elmat+= Nc - Trans(Nc) + lam*(ps+1)*ps/h * Ns; 
      //elmat+= lam*(ps+1)*ps/h * Ns; 
      // FlatMatrix<> lapmat(ndof_total,ndof_total,lh);
      // laplace->CalcElementMatrix(

      // cout<<"elmat = "<<elmat<<endl;
      // getchar();

      // Vector<double> lami(elmat.Height());
      // Matrix<double> evecs(elmat.Height());
                        
      // CalcEigenSystem (elmat, lami, evecs);
      // cout << "lami = " << endl << lami << endl;
      // getchar();
      
    }
    
  }

  template class XStokesNitscheIntegrator<2>;
  template class XStokesNitscheIntegrator<3>;

  static RegisterBilinearFormIntegrator<XStokesNitscheIntegrator<2>> inisnitschestokes2d("xstokesnitsche",2,3);
  static RegisterBilinearFormIntegrator<XStokesNitscheIntegrator<3>> inisnitschestokes3d("xstokesnitsche",3,3);

  template <int D>
  void XStokesNitscheRhsIntegrator<D>::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("XStokesNitscheModLBIntegrator::CalcElementMatrix");
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

    elvec = 0.0;

    if(!feuvx)
      return;

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

    IntRange dofrangep_x(cnt, cnt+ndofp_x);

    int ndof_total = ndofp_total + D*ndofuv_total;

    FlatVector<> shapeuv_total(ndofuv_total,lh);
    FlatVector<> shapeuv(ndofuv,&shapeuv_total(0));
    FlatVector<> shapeuvx(ndofuv,&shapeuv_total(ndofuv));
    
    int orderuv = feuv.Order();
    int orderp = fep.Order();
    int ps = orderuv;

    Mat<D> mat = 0;

    FlatVector<> bmatjump(ndof_total,lh);

    bmatjump = 0.0;

    const FlatArray<DOMAIN_TYPE>& xusign = feuvx->GetSignsOfDof();

    const FlatXLocalGeometryInformation & xgeom(feuvx->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());
    
    double kappa_neg;
    double kappa_pos;

    if (xgeom.kappa[NEG] >= 0.5)
      kappa_neg = 1.0;
    else
      kappa_neg = 0.0;
    kappa_pos = 1.0 - kappa_neg;

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

      shapeuv = feuv.GetShape(mip.IP(),lh);
      
      bmatjump = 0.0;
      for (int d = 0; d<D; ++d)
      {
        bmatjump.Range(*dofrangeuv[d])= normal[d] * shapeuv;
        bmatjump.Range(*dofrangeuv_x[d])= normal[d] * shapeuv;
      }  

      // cout << "before: bmatjump = \n" << bmatjump << endl; getchar();
      
      for (int d=0; d<D; ++d)
      {
        for (int l = 0; l < ndofuv_x; ++l)
        {
          if (xusign[l] == NEG){
            bmatjump(dofrangeuv[d]->Next()+l) *= kappa_pos;
          }
          else{
            bmatjump(dofrangeuv[d]->Next()+l) *= kappa_neg;
          }
        }
      }

      // cout << "after: bmatjump = \n" << bmatjump << endl; getchar();
      
      elvec+= weight * coef_val * bmatjump; 
    }
    
  }

  template class XStokesNitscheRhsIntegrator<2>;
  template class XStokesNitscheRhsIntegrator<3>;

  static RegisterLinearFormIntegrator<XStokesNitscheRhsIntegrator<2>> inisnitschestokesrhs2d("xstokesnitscherhs",2,1);
  static RegisterLinearFormIntegrator<XStokesNitscheRhsIntegrator<3>> inisnitschestokesrhs3d("xstokesnitscherhs",3,1);

  template <int D>
  void XStokesNitscheModLBIntegrator<D>::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("XStokesNitscheRhsIntegrator::CalcElementMatrix");
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

    elvec = 0.0;

    if(!feuvx)
      return;

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

    IntRange dofrangep_x(cnt, cnt+ndofp_x);

    int ndof_total = ndofp_total + D*ndofuv_total;

    FlatMatrixFixWidth<D> dshapeuv_total(ndofuv_total,lh);
    FlatMatrixFixWidth<D> dshapeuv(ndofuv,&dshapeuv_total(0,0));
    FlatMatrixFixWidth<D> dshapeuvx(ndofuv,&dshapeuv_total(ndofuv,0));
    FlatMatrixFixWidth<D*D> bmat(D*ndofuv_total,lh);
    FlatMatrixFixWidth<D*D> bmat_P(D*ndofuv_total,lh);
    FlatVector<> P2_bmat_P(D*ndofuv_total,lh);
    
    int orderuv = feuv.Order();
    int orderp = fep.Order();
    int ps = orderuv;

    Mat<D> mat = 0;

    const FlatArray<DOMAIN_TYPE>& xusign = feuvx->GetSignsOfDof();

    const FlatXLocalGeometryInformation & xgeom(feuvx->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    IntegrationPoint ipc(0.0,0.0,0.0);
    MappedIntegrationPoint<D,D> mipc(ipc, eltrans);

    const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());
    
    double kappa_neg;
    double kappa_pos;

    if (xgeom.kappa[NEG] >= 0.5)
      kappa_neg = 1.0;
    else
      kappa_neg = 0.0;
    kappa_pos = 1.0 - kappa_neg;

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
      
      const double coef_val = coef_sigma->Evaluate(mip);

      Vec<D> improved_normal = normal;
      // coef_lset num. diff
      if (improved) //improved)
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
        improved_normal /= - L2Norm(improved_normal);
      }

      // cout << " normal = " << normal << endl;
      // cout << " improved_normal = " << improved_normal << endl;
      // getchar();
      
      
      feuv.CalcMappedDShape (mip,dshapeuv);
      bmat = 0.0;
      for (int d = 0; d < D; ++d)
      {
        bmat.Rows(*dofrangeuv[d]).Cols(D*d,D*(d+1)) = dshapeuv;
        bmat.Rows(*dofrangeuv_x[d]).Cols(D*d,D*(d+1)) = dshapeuv;
      }

      // cout << " bmat = " << bmat << endl;

      Vec<D*D> P_vec;
      FlatMatrix<> P(D,D,&P_vec(0));
      P = Id<D>() - normal * Trans(normal);

      Mat<D*D> project;
      project = 0.0;
      for (int d = 0; d < D; ++d)
        project.Rows(d*D,(d+1)*D).Cols(d*D,(d+1)*D) = P;
      
      Vec<D*D> P_impr_vec;
      FlatMatrix<> P_impr(D,D,&P_impr_vec(0));
      P_impr = Id<D>() - improved_normal * Trans(improved_normal);
      
      bmat_P = bmat * project; 

      P2_bmat_P = bmat * P_impr_vec; 

      for (int d=0; d<D; ++d)
      {
        for (int l = 0; l < ndofuv_x; ++l)
        {
          if (xusign[l] == NEG){
            P2_bmat_P(dofrangeuv[d]->Next()+l) *= kappa_pos;
          }
          else{
            P2_bmat_P(dofrangeuv[d]->Next()+l) *= kappa_neg;
          }
        }
      }

      elvec.Range(0,D*ndofuv_total) -= weight * coef_val * P2_bmat_P;
    }
  }

  template class XStokesNitscheModLBIntegrator<2>;
  template class XStokesNitscheModLBIntegrator<3>;

  static RegisterLinearFormIntegrator<XStokesNitscheModLBIntegrator<2>> inisnitschestokesmodlb2d("xstokesnitschemodlb",2,2);
  static RegisterLinearFormIntegrator<XStokesNitscheModLBIntegrator<3>> inisnitschestokesmodlb3d("xstokesnitschemodlb",3,2);

}




