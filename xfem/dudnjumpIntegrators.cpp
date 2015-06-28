/*********************************************************************/
/* File:   CDR_Integrators.cpp                                       */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   09. Jan. 2011                                             */
/*********************************************************************/

/*
  

 */


#include <fem.hpp>


namespace ngfem
{

  void MySetRefBaryCenter(ELEMENT_TYPE eltype, FlatVector<> & point)
  {
    switch (eltype)
    {
    case ET_POINT: point = 0.0; break;
    case ET_SEGM: point = 0.5; break;
    case ET_TRIG: point = 1.0/3.0; break;
    case ET_QUAD: point = 0.5; break;
    case ET_TET: point = 1.0/4.0; break;
    case ET_PYRAMID: point(0) = 1.0/2.0; point(1) = 1.0/2.0; point(2) = 1.0/4.0; break;
    case ET_PRISM: point(0) = 1.0/3.0; point(1) = 1.0/3.0; point(2) = 1.0/2.0; break;
    case ET_HEX: point = 1.0/2.0; break;
    };
  };


  template <int D>
  class dudnJumpIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef;
  public:
    dudnJumpIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetBilinearFormIntegrator(coeffs),coef(coeffs[0])
    { 
    }

    virtual ~dudnJumpIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual bool IsSymmetric () const 
    { return true; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new dudnJumpIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("dudnJumpIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                  const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                  const FiniteElement & volumefel2, int LocalFacetNr2,
                                  const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("dudnJumpIntegratorIntegrator");

      if (LocalFacetNr2==-1) throw Exception("dudnJumpIntegrator: LocalFacetNr2==-1");

      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel1);
      ELEMENT_TYPE eltype1 = volumefel1.ElementType();
      int nd1 = fel1_l2->GetNDof();

      const ScalarFiniteElement<D> * fel2_l2 = NULL;
      ELEMENT_TYPE eltype2 = eltype1;

      fel2_l2 = dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel2);
      eltype2 = volumefel2.ElementType();
      int nd2 = fel2_l2->GetNDof();
      double maxorder = max(fel1_l2->Order(),fel2_l2->Order());

      FlatVector<double> eltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype1,eltype1_ref_barcenter);
      IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans1);
      Vec<D> pointv = barycenterv.GetPoint();
      FlatVector<double> eltype2_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype2,eltype2_ref_barcenter);
      IntegrationPoint refbarycenterv2(eltype2_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv2 (refbarycenterv2, eltrans2);
      Vec<D> pointv2 = barycenterv2.GetPoint();
      Vec<D> pointsdiff = pointv - pointv2;
      
      // *testout << " pointv  = " << pointv << endl;
      // *testout << " pointv2 = " << pointv2 << endl;


      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      FlatVector<> mat2_shape(nd2, lh);
      FlatVector<> mat2_dudn(nd2, lh);
      
      FlatMatrixFixHeight<1> bmat(nd1+nd2, lh);
      FlatMatrixFixHeight<1> dbmat(nd1+nd2, lh);
      Mat<1> dmat;

      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1,ElVertices1); 
      Facet2ElementTrafo transform2(eltype2,ElVertices2); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
        normal_ref1(i) = normals1[LocalFacetNr1][i];
        normal_ref2(i) = normals2[LocalFacetNr2][i];
      }
      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*maxorder);
	
      if (maxorder==0) maxorder=1;
   
      bmat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
      {
        IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
	  
        MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);
        // double lam = coef->Evaluate(sip1);

        // Mat<D> jac1 = sip1.GetJacobian();
        Mat<D> inv_jac1 = sip1.GetJacobianInverse();
        double det1 = sip1.GetJacobiDet();

        Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
        double len1 = L2Norm (normal1);
        normal1 /= len1;

        fel1_l2->CalcShape(sip1.IP(), mat1_shape);
        Vec<D> invjac_normal1 = inv_jac1 * normal1;
        mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
        IntegrationPoint ip2 = (LocalFacetNr2!=-1) ? transform2(LocalFacetNr2, ir_facet[l]) : ip1;
        MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);
        // double lam2 = coef_lam->Evaluate(sip2);
        // Mat<D> jac2 = sip2.GetJacobian();
        Mat<D> inv_jac2 = sip2.GetJacobianInverse();
        double det2 = sip2.GetJacobiDet();
	  
        Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;       
        double len2 = L2Norm (normal2); 
        if(abs(len1-len2)>1e-6){
          std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
          throw Exception ("dudnJumpIntegrator: len1!=len2");
        }
        normal2 /= len2;
        Vec<D> invjac_normal2;;
        fel2_l2->CalcShape(sip2.IP(), mat2_shape);
        invjac_normal2 = inv_jac2 * normal2;
        mat2_dudn = fel2_l2->GetDShape (sip2.IP(), lh) * invjac_normal2;
        bmat.Row(0).Range (0   , nd1)   = mat1_dudn;	    
        bmat.Row(0).Range (nd1   , nd1+nd2)   = mat2_dudn;

        dmat(0,0) = 1.0;

        const double orthdist = abs(InnerProduct(pointsdiff,normal1));

        *testout << " orthdist h = " << orthdist << endl;
        // *testout << "       len1 = " << len1 << endl;

        dmat *= coef->Evaluate(sip1) * len1 / orthdist * ir_facet[l].Weight();
        dbmat = dmat * bmat;
        elmat += Trans (bmat) * dbmat;
        // elmat = 1.0;
      }
      if (LocalFacetNr2==-1) elmat=0.0;
    }
  };

  // \int_F du/dn dv/dn on boundary facets... 
  template <int D>
  class dudnboundBLIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef;
  public:
    dudnboundBLIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetBilinearFormIntegrator(coeffs),coef(coeffs[0])
    { 
    }

    virtual ~dudnboundBLIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 1; }

    virtual bool IsSymmetric () const 
    { return true; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new dudnboundBLIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("dudnboundBLIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                  const ElementTransformation & seltrans,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
    {
      
      static int timer = NgProfiler::CreateTimer ("dudnboundBLIntegratorIntegrator");

      if (D==3) throw Exception("dudnboundBLIntegratorIntegrator - scaling might be wrong");

      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel);
      ELEMENT_TYPE eltype1 = volumefel.ElementType();
      ELEMENT_TYPE seltype1 = seltrans.GetElementType();
      int nd1 = fel1_l2->GetNDof();


      FlatVector<double> eltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype1,eltype1_ref_barcenter);
      IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans);
      Vec<D> pointv = barycenterv.GetPoint();

      FlatVector<double> seltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(seltype1,seltype1_ref_barcenter);
      IntegrationPoint refbarycenterf(seltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D-1,D> barycenterf (refbarycenterf, seltrans);
      Vec<D> pointf = barycenterf.GetPoint();
      Vec<D> pointsdiff = pointv - pointf;

      double maxorder = fel1_l2->Order();
      
      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      
      FlatMatrixFixHeight<1> bmat(nd1, lh);
      FlatMatrixFixHeight<1> dbmat(nd1, lh);
      Mat<1> dmat;

      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1,ElVertices); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

      Vec<D> normal_ref1;
      for (int i=0; i<D; i++){
        normal_ref1(i) = normals1[LocalFacetNr][i];
      }
      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*maxorder);
	
      if (maxorder==0) maxorder=1;
   
      bmat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
      {
        IntegrationPoint ip1 = transform1(LocalFacetNr, ir_facet[l]);
	  
        MappedIntegrationPoint<D,D> sip1 (ip1, eltrans);
        // double lam = coef->Evaluate(sip1);

        // Mat<D> jac1 = sip1.GetJacobian();
        Mat<D> inv_jac1 = sip1.GetJacobianInverse();
        double det1 = sip1.GetJacobiDet();

        Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
        double len1 = L2Norm (normal1);
        normal1 /= len1;

        fel1_l2->CalcShape(sip1.IP(), mat1_shape);
        Vec<D> invjac_normal1 = inv_jac1 * normal1;
        mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
        bmat.Row(0).Range (0   , nd1)   = mat1_dudn;	    

        dmat(0,0) = 1.0;

        const double orthdist = 1.0 * abs(InnerProduct(pointsdiff,normal1));

        *testout << " orthdist h = " << orthdist << endl;
        // *testout << "       len1 = " << len1 << endl;

        dmat *= coef->Evaluate(sip1) * len1 / orthdist* ir_facet[l].Weight();
        dbmat = dmat * bmat;
        elmat+= Trans (bmat) * dbmat;
      }
    }
  };


  // \int_F g dv/dn on boundary facets... 
  template <int D>
  class dudnboundLFIntegrator : public FacetLinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef;
  public:
    dudnboundLFIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetLinearFormIntegrator(coeffs), coef(coeffs[0])
    { 
    }

    virtual ~dudnboundLFIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 1; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new dudnboundLFIntegrator (coeffs);
    }
    
    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elvec,
                                    LocalHeap & lh) const
    {
      throw Exception("dudnboundLFIntegrator::CalcElementVector - not implemented!");
    }
    
    virtual void CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
                                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                  const ElementTransformation & seltrans,
                                  FlatVector<double> & elvec,
                                  LocalHeap & lh) const
    {
      
      static int timer = NgProfiler::CreateTimer ("dudnboundLFIntegratorIntegrator");

      if (D==3) throw Exception("dudnboundLFIntegratorIntegrator - scaling might be wrong");

      NgProfiler::RegionTimer reg (timer);
      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel);
      ELEMENT_TYPE eltype1 = volumefel.ElementType();
      ELEMENT_TYPE seltype1 = seltrans.GetElementType();
      int nd1 = fel1_l2->GetNDof();


      FlatVector<double> eltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype1,eltype1_ref_barcenter);
      IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans);
      Vec<D> pointv = barycenterv.GetPoint();

      FlatVector<double> seltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(seltype1,seltype1_ref_barcenter);
      IntegrationPoint refbarycenterf(seltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D-1,D> barycenterf (refbarycenterf, seltrans);
      Vec<D> pointf = barycenterf.GetPoint();
      Vec<D> pointsdiff = pointv - pointf;
      

      double maxorder = fel1_l2->Order();
      
      elvec= 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      
      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1,ElVertices); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

      Vec<D> normal_ref1;
      for (int i=0; i<D; i++){
        normal_ref1(i) = normals1[LocalFacetNr][i];
      }
      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*maxorder);
	
      if (maxorder==0) maxorder=1;
   
      for (int l = 0; l < ir_facet.GetNIP(); l++)
      {
        IntegrationPoint ip1 = transform1(LocalFacetNr, ir_facet[l]);
	  
        MappedIntegrationPoint<D,D> sip1 (ip1, eltrans);
        // double lam = coef->Evaluate(sip1);

        // Mat<D> jac1 = sip1.GetJacobian();
        Mat<D> inv_jac1 = sip1.GetJacobianInverse();
        double det1 = sip1.GetJacobiDet();

        Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
        double len1 = L2Norm (normal1);
        normal1 /= len1;
          
        const double orthdist = 1.0 * abs(InnerProduct(pointsdiff,normal1));

        fel1_l2->CalcShape(sip1.IP(), mat1_shape);
        Vec<D> invjac_normal1 = inv_jac1 * normal1;
        mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
	  
        *testout << " orthdist h = " << orthdist << endl;
        // *testout << "       len1 = " << len1 << endl;

        const double fac = coef->Evaluate(sip1) * len1 / orthdist * ir_facet[l].Weight();
        elvec += fac * mat1_dudn;
      }
    }
  };

  static RegisterBilinearFormIntegrator<dudnJumpIntegrator<2> > init_gp_2d_1 ("dudnjumpdvdnjump", 2, 1);
  static RegisterBilinearFormIntegrator<dudnboundBLIntegrator<2> > init_gp_2d_2 ("dudndvdn", 2, 1);
  static RegisterLinearFormIntegrator<dudnboundLFIntegrator<2> > init_gp_2d_3 ("dudndvdn", 2, 1);
}
