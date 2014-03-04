#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"
// #include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{

  void SetRefBaryCenter(ELEMENT_TYPE eltype, FlatVector<> & point)
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
  class LowOrderGhostPenaltyIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_lam_neg;
    CoefficientFunction *coef_lam_pos;
  public:
    LowOrderGhostPenaltyIntegrator (Array<CoefficientFunction*> & coeffs) 
      : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_lam_neg  = coeffs[0];
      coef_lam_pos  = coeffs[1];
    }

    virtual ~LowOrderGhostPenaltyIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 0; }
    
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new LowOrderGhostPenaltyIntegrator (coeffs);
    }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("LowOrderGhostPenaltyIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                  const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                  const FiniteElement & volumefel2, int LocalFacetNr2,
                                  const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh// ,
                                  // BitArray* twice
      ) const
    {
      static int timer = NgProfiler::CreateTimer ("LowOrderGhostPenaltyIntegrator");

      if (LocalFacetNr2==-1) throw Exception("LowOrderGhostPenaltyIntegrator: LocalFacetNr2==-1");
      if (D==3) throw Exception("scaling is not correct!");
      NgProfiler::RegionTimer reg (timer);

      const CompoundFiniteElement * dcfel1 = 
        dynamic_cast<const CompoundFiniteElement*> (&volumefel1);
      const CompoundFiniteElement * dcfel2 = 
        dynamic_cast<const CompoundFiniteElement*> (&volumefel2);
      if (dcfel1==NULL || dcfel2==NULL){
        cout << "not compound!" << endl;
        cout << " leaving " << endl;
        return;
      }
      
      const CompoundFiniteElement & cfel1 = *dcfel1;       
      const CompoundFiniteElement & cfel2 = *dcfel2;       

      const XFiniteElement * xfe1 = NULL;
      const ScalarFiniteElement<D> * scafe1 = NULL;
      const XFiniteElement * xfe2 = NULL;
      const ScalarFiniteElement<D> * scafe2 = NULL;
      
      for (int i = 0; i < cfel1.GetNComponents(); ++i)
      {
        if (xfe1==NULL)
          xfe1 = dynamic_cast<const XFiniteElement* >(&cfel1[i]);
        if (scafe1==NULL)
          scafe1 = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel1[i]);
      }

      for (int i = 0; i < cfel2.GetNComponents(); ++i)
      {
        if (xfe2==NULL)
          xfe2 = dynamic_cast<const XFiniteElement* >(&cfel2[i]);
        if (scafe2==NULL)
          scafe2 = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel2[i]);
      }
    
      elmat = 0.0;

      if (!xfe1 || !xfe2){ // not a ghost edge
        return;
      }

      ELEMENT_TYPE eltype1 = volumefel1.ElementType();
      int nd1 = volumefel1.GetNDof();
      int ndof_x1 = xfe1->GetNDof();
      int ndof_sca1 = scafe1->GetNDof();
      FlatVector<> mat1_dudn(nd1, lh);
      FlatVector<> mat1_dudn_sca(ndof_sca1, &mat1_dudn(0));
      FlatVector<> mat1_dudn_x(ndof_x1, &mat1_dudn(ndof_sca1));

      ELEMENT_TYPE eltype2 = volumefel2.ElementType();
      int nd2 = volumefel2.GetNDof();
      int ndof_x2 = xfe2->GetNDof();
      int ndof_sca2 = scafe2->GetNDof();
      FlatVector<> mat2_dudn(nd2, lh);
      FlatVector<> mat2_dudn_sca(ndof_sca2, &mat2_dudn(0));
      FlatVector<> mat2_dudn_x(ndof_x2, &mat2_dudn(ndof_sca2));

      int maxorder = max(scafe1->Order(),scafe2->Order());
      if (maxorder==0) maxorder=1;

      FlatMatrixFixWidth<1> bmat(nd1+nd2, lh);

      // const MasterElement & masterel1 = xfe1->GetMasterElement();
      // const MasterElement & masterel2 = xfe2->GetMasterElement();

      // Edge * commonedge = NULL;
      // for (int i = 0; i < masterel1.NBaseEdges(); ++i)
      //   for (int j = 0; j < masterel2.NBaseEdges(); ++j)
      //     if (masterel1.GetBaseEdge(i) == masterel2.GetBaseEdge(j)){
      //       commonedge = masterel1.GetBaseEdge(i);
      //     }

      Facet2ElementTrafo transform1(eltype1,ElVertices1); 
      Facet2ElementTrafo transform2(eltype2,ElVertices2); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      HeapReset hr(lh);

      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

      FlatVector<double> eltype1_ref_barcenter(3,lh);
      SetRefBaryCenter(eltype1,eltype1_ref_barcenter);
      IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans1);
      Vec<D> pointv = barycenterv.GetPoint();
      FlatVector<double> eltype2_ref_barcenter(3,lh);
      SetRefBaryCenter(eltype2,eltype2_ref_barcenter);
      IntegrationPoint refbarycenterv2(eltype2_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv2 (refbarycenterv2, eltrans2);
      Vec<D> pointv2 = barycenterv2.GetPoint();
      Vec<D> pointsdiff = pointv - pointv2;

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
        normal_ref1(i) = normals1[LocalFacetNr1][i];
        normal_ref2(i) = normals2[LocalFacetNr2][i];
      }

      const int p = maxorder;

      const FlatArray<DOMAIN_TYPE>& xsign1 = xfe1->GetSignsOfDof();
      const FlatArray<DOMAIN_TYPE>& xsign2 = xfe2->GetSignsOfDof();
      
      DOMAIN_TYPE dt = POS;
      for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
      {

        const IntegrationRule & ir_facet =
          SelectIntegrationRule (etfacet, 2*p);

        // IntegrationRule ir_facet;
        // commonedge->FillSurfaceIntegrationRule(2*p,sign,ir_facet);
        bmat = 0.0;

        for (int l = 0; l < ir_facet.GetNIP(); l++)
        {
          IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
          MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);

          double lam = 0;
          if (dt == POS)
            lam = coef_lam_pos->Evaluate(sip1);
          else
            lam = coef_lam_neg->Evaluate(sip1);
          Mat<D> inv_jac1 = sip1.GetJacobianInverse();
          double det1 = sip1.GetJacobiDet();

          Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
          double len1 = L2Norm (normal1);
          normal1 /= len1;
          Vec<D> invjac_normal1 = inv_jac1 * normal1;
          mat1_dudn_sca = scafe1->GetDShape (sip1.IP(), lh) * invjac_normal1;
          mat1_dudn_x = mat1_dudn_sca;
          for (int i = 0; i < ndof_sca1; ++i)
            if (xsign1[i]!=dt)
              mat1_dudn_x(i) = 0;              
          // mat1_dudn_x = xfe1->GetDShape (sip1.IP(), lh) * invjac_normal1;

          const double orthdist = abs(InnerProduct(pointsdiff,normal1));

          IntegrationPoint ip2 = transform2(LocalFacetNr2, ir_facet[l]);
          MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);
          Mat<D> inv_jac2 = sip2.GetJacobianInverse();
          double det2 = sip2.GetJacobiDet();
          Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;       
          double len2 = L2Norm (normal2); 
          if(abs(len1-len2)>1e-6){
            std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
            throw Exception ("LowOrderGhostPenaltyIntegrator: len1!=len2");
          }
          normal2 /= len2;
          Vec<D> invjac_normal2;;
          invjac_normal2 = inv_jac2 * normal2;
          mat2_dudn_sca = scafe2->GetDShape (sip2.IP(), lh) * invjac_normal2;
          mat2_dudn_x = mat2_dudn_sca;
          for (int i = 0; i < ndof_sca2; ++i)
            if (xsign2[i]!=dt)
              mat2_dudn_x(i) = 0;              
          // mat2_dudn_x = xfe2->GetDShape (sip2.IP(), lh) * invjac_normal2;

          // for (int k = 0; k < ndof_x1; ++k)
          //   if (xsign1[k] != sign)
          //     mat1_dudn_x(k) = 0;
          // for (int k = 0; k < ndof_x2; ++k)
          //   if (xsign2[k] != sign)
          //     mat2_dudn_x(k) = 0;

          bmat = 0.0;
          bmat.Col(0).Range(  0,    nd1) = mat1_dudn;
          bmat.Col(0).Range(nd1,nd1+nd2) = mat2_dudn;

          elmat += lam * orthdist * len1 * ir_facet[l].Weight() * bmat * Trans (bmat);

        }
      }
      
    }

  };

  template class LowOrderGhostPenaltyIntegrator<2>;
  template class LowOrderGhostPenaltyIntegrator<3>;

  static RegisterBilinearFormIntegrator<LowOrderGhostPenaltyIntegrator<2> > init_gp_2d ("lo_ghostpenalty", 2, 2);
  static RegisterBilinearFormIntegrator<LowOrderGhostPenaltyIntegrator<3> > init_gp_3d ("lo_ghostpenalty", 3, 2);

}
