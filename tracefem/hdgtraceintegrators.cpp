#define FILE_HDGTRACEINTEGRATORS_CPP
#include "traceintegrators.hpp"
#include "hdgtraceintegrators.hpp"

namespace ngfem
{

  template<int D>
  void HDGTraceLaplaceBeltramiIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & ccfel,
                     const ElementTransformation & eltrans,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("HDGTraceLBInt::CalcElMat");
    RegionTimer reg (timer);

    static Timer timer_facet_prepare1 ("HDGTraceLBInt::facetprepare 1");
    static Timer timer_facet_prepare2 ("HDGTraceLBInt::facetprepare 2");
    static Timer timer_facet_prepare3 ("HDGTraceLBInt::facetprepare 3");
    static Timer timer_facet_prepare4 ("HDGTraceLBInt::facetprepare 4");
    
    auto & cfel = static_cast<const CompoundFiniteElement&> (ccfel);
    
    const XFiniteElement * tracel2fe =
      dynamic_cast<const XFiniteElement *> (&cfel[0]);

    const XFiniteElement * tracefacetfe =
      dynamic_cast<const XFiniteElement *> (&cfel[1]);
    
    elmat = 0.0;
    if (!tracel2fe) return;
    if (!tracefacetfe) return;

    static Timer timer2 ("HDGTraceLBInt::CalcElementMatrix (cut el.)");
    RegionTimer reg2 (timer2);
    
    auto & fel_l2 = dynamic_cast<const ScalarFiniteElement<D> & >(tracel2fe->GetBaseFE());
    auto & fel_facet = dynamic_cast<const FacetVolumeFiniteElement<D> &> (tracefacetfe->GetBaseFE());

    int order = fel_l2.Order();
    
    ELEMENT_TYPE eltype = cfel.ElementType();

    if (eltype != ET_TET)
      throw Exception(" oh damn,... I really expected a tetra here..");
    

    Vec<3> interface_normal_ref;
    { // determine interface_normal_ref
      static Timer timera ("HDGTraceLBInt::determine_interface_normal_ref");
      RegionTimer reg (timera);
      
      const POINT3D * verts = ElementTopology::GetVertices(eltype);
      const int nv =  ElementTopology::GetNVertices(eltype);
      FlatVector<> lset_vvals(nv,lh);
      for (int v = 0; v < nv; ++v)
      {
        IntegrationPoint vip(&verts[v][0],0);
        lset_vvals(v) = tracel2fe->GetFlatLocalGeometry().EvaluateLsetAtPoint<D,D>(vip);
      }
      auto tet_p1_fe = new (lh) ScalarFE<ET_TET,1>;
      FlatMatrixFixWidth<3> dshape_ref(nv,lh);
      IntegrationPoint ip(0.0,0.0,0.0,0.0);
      tet_p1_fe->CalcDShape(ip,dshape_ref);
      interface_normal_ref = Trans(dshape_ref)*lset_vvals;
    }

    auto trig_p1_fe = new (lh) ScalarFE<ET_TRIG,1>;

    IntRange l2_dofs = cfel.GetRange(0);
    IntRange facet_dofs = cfel.GetRange(1);

    int nd_l2 = tracel2fe->GetNDof();
    int nd_facet = tracefacetfe->GetNDof();
    int nd = cfel.GetNDof();  

    int base_l2 = 0;
    int base_facet = base_l2+nd_l2;

    
    FlatVector<> mat_l2(nd_l2, lh);
    FlatVector<> mat_dudn(nd_l2, lh);
    FlatVector<> mat_facet(nd_facet, lh);
    FlatMatrix<> mat_gradgrad (nd_l2, lh);

    FlatMatrix<> mat_coupling(nd_l2, nd, lh);
    FlatVector<> jump(nd, lh);
    mat_coupling = 0.0;

 

    // The facet contribution
    {
      int nfacet = ElementTopology::GetNFacets(eltype);
      
      Facet2ElementTrafo transform(eltype); 
      FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

      for (int k = 0; k < nfacet; k++)
      {
        timer_facet_prepare1.Start();
        HeapReset hr(lh);
        ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        Vec<D> facet_normal_ref = normals[k];

        if (etfacet != ET_TRIG)
          throw Exception(" oh damn,... I really expected a triangle here..");

        const POINT3D * verts = ElementTopology::GetVertices(etfacet);
        const int nv =  ElementTopology::GetNVertices(etfacet);
        FlatVector<> lset_vvals(nv,lh);
        timer_facet_prepare1.Stop();
        
        timer_facet_prepare2.Start();
        for (int v = 0; v < nv; ++v)
        {
          IntegrationPoint fip(&verts[v][0],0);
          IntegrationPoint vip = transform(k, fip);
          lset_vvals(v) = tracel2fe->GetFlatLocalGeometry().EvaluateLsetAtPoint<D,D>(vip);
        }
        timer_facet_prepare2.Stop();
        
        timer_facet_prepare3.Start();

        auto lset_eval = ScalarFieldEvaluator::Create(D-1,*trig_p1_fe,lset_vvals,lh);
        shared_ptr<XLocalGeometryInformation> xgeom = nullptr;
        CompositeQuadratureRule<D-1> cquad2d;
        xgeom = XLocalGeometryInformation::Create(etfacet, ET_POINT,
                                                  *lset_eval, cquad2d, lh,
                                                  2*order, 0, 0, 0); // at least order - 1 is needed here...
        DOMAIN_TYPE facet_domain;

        {
          static Timer timerb2 ("HDGTraceLBInt::make cutgeom on facet");
          RegionTimer reg (timerb2);
          facet_domain = xgeom->MakeQuadRule();
        }

        timer_facet_prepare3.Stop();
        timer_facet_prepare4.Start();

        Array<int> comp_facetdofs;
        comp_facetdofs += l2_dofs;
        comp_facetdofs += fel_facet.GetFacetDofs(k) + base_facet;

        timer_facet_prepare4.Stop();

        if (facet_domain == IF)
        {
          static Timer timerb ("HDGTraceLBInt::edge part");
          RegionTimer reg (timerb);
          HeapReset hr(lh);
          
          const QuadratureRuleCoDim1<D-1> & edge_quad(cquad2d.GetInterfaceRule());
	    
          FlatMatrix<> comp_jumps(edge_quad.Size(), comp_facetdofs.Size(), lh);
          FlatMatrix<> comp_facjumps(edge_quad.Size(), comp_facetdofs.Size(), lh);
          FlatMatrix<> facdudn(edge_quad.Size(), nd_l2, lh);
          
          for (int l = 0; l < edge_quad.Size(); ++l)
          {
            IntegrationPoint fip(&edge_quad.points[l](0),edge_quad.weights[l]);

            Vec<D-1> tang_ref;
            tang_ref[0]=   edge_quad.normals[l][1];
            tang_ref[1]= - edge_quad.normals[l][0];
            tang_ref /= L2Norm(tang_ref);

            FlatMatrix<> facettrafoJ = transform.GetJacobian(k,lh);
            Vec<D> tang_vol_ref = facettrafoJ * tang_ref;
         
            IntegrationPoint vip = transform(k, fip);

            MappedIntegrationPoint<3,3> mip(vip, eltrans);

            double alpha = coef_alpha -> Evaluate (mip);
            
            Vec<3> tang = mip.GetJacobian() * tang_vol_ref;
            double measure_ratio = L2Norm(tang);
            tang /= measure_ratio;

            double weight = edge_quad.weights[l] * measure_ratio;

            Mat<D> inv_jac = mip.GetJacobianInverse();
            double h = cbrt(mip.GetMeasure());

            Vec<D> facet_normal = Trans (inv_jac) * facet_normal_ref;       
            facet_normal /= L2Norm (facet_normal);

            Vec<3> interface_normal = Trans(inv_jac) * interface_normal_ref ;
            interface_normal /= L2Norm(interface_normal);

            // normal to edge and normal to interface...
            Vec<3> conormal = Cross(tang,interface_normal);

            mat_facet = 0.0;
            fel_facet.Facet(k).CalcShape (vip,mat_facet.Range(fel_facet.GetFacetDofs(k)));
          
            fel_l2.CalcShape(vip, mat_l2);

            Vec<D> invjac_conormal = inv_jac * conormal;
            mat_dudn = fel_l2.GetDShape (mip.IP(), lh) * invjac_conormal;

            jump.Range(l2_dofs) = mat_l2;
            jump.Range(facet_dofs) = - alpha * mat_facet;
          
            comp_jumps.Row(l) = jump(comp_facetdofs);
            double fac = param_IP_edge*weight*alpha;
            comp_facjumps.Row(l) = (fac * (order+1) * (order+1) * (order+1) * (order+1)/h) * jump(comp_facetdofs);
            facdudn.Row(l) = (-weight) * mat_dudn;

          }

          FlatMatrix<> comp_elmat(comp_facetdofs.Size(), comp_facetdofs.Size(), lh);

          // penalty term:
          comp_elmat = Trans (comp_facjumps) * comp_jumps | Lapack;
          elmat.Rows(comp_facetdofs).Cols(comp_facetdofs) += comp_elmat;

          FlatMatrix<> comp_elmat2(nd_l2, comp_facetdofs.Size(), lh);

          comp_elmat2 = Trans(facdudn) * comp_jumps | Lapack;
          mat_coupling.Cols(comp_facetdofs) += comp_elmat2;
          
        }

        {
          static Timer timerc ("HDGTraceLBInt::facet part (normal diffusion)");
          RegionTimer reg (timerc);
          HeapReset hr(lh);
          
          IntegrationRule ir_facet(etfacet, 2*order);
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
	    
          FlatMatrix<> comp_jumps(ir_facet.GetNIP(), comp_facetdofs.Size(), lh);
          FlatMatrix<> comp_facjumps(ir_facet.GetNIP(), comp_facetdofs.Size(), lh);
	    
          MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);
          
	    for (int l = 0; l < ir_facet.GetNIP(); l++)
	      {
	        MappedIntegrationPoint<D,D> & mip = mir[l];

                double alpha = coef_alpha -> Evaluate (mip);
                
	        Mat<D> inv_jac = mip.GetJacobianInverse();

                double param_normaldiff = coef_param_normaldiff -> Evaluate (mip);

                double det = mip.GetMeasure();

	        Vec<D> facet_normal = det * Trans (inv_jac) * facet_normal_ref;       
	        double len = L2Norm (facet_normal);
	        facet_normal /= len;

	        mat_facet = 0.0;
	        fel_facet.Facet(k).CalcShape (ir_facet_vol[l], 
	        			      mat_facet.Range(fel_facet.GetFacetDofs(k)));
	        fel_l2.CalcShape(ir_facet_vol[l], mat_l2);

	        jump.Range(l2_dofs) = mat_l2;
	        jump.Range(facet_dofs) = -mat_facet;

	        comp_jumps.Row(l) = jump(comp_facetdofs);
	        double h = cbrt(mip.GetMeasure());
	        double weight = len * ir_facet[l].Weight() * 1.0/h;
	        comp_facjumps.Row(l) = (param_IP_facet*param_normaldiff*alpha
                                        *weight*(order+1)*(order+1)/h)
                  *jump(comp_facetdofs);
	      }
          // facet penalty term:
          FlatMatrix<> comp_elmat(comp_facetdofs.Size(), comp_facetdofs.Size(), lh);
          comp_elmat = Trans (comp_facjumps) * comp_jumps | Lapack;
          elmat.Rows(comp_facetdofs).Cols(comp_facetdofs) += comp_elmat;

        }
      }
      
      elmat.Rows(l2_dofs) += mat_coupling;
      elmat.Cols(l2_dofs) += Trans(mat_coupling);
    }


    {
      static Timer timerd ("HDGTraceLBInt::volume part (normal diffusion)");
      RegionTimer reg (timerd);
      
      FlatMatrixFixWidth<D> dshape(nd_l2,lh);
      FlatMatrixFixWidth<D> dshape_ref(nd_l2,lh);
      FlatVector<> dshapen(nd_l2,lh);
    
      IntegrationRule ir(eltype, 2*order-2);

      // loop over integration points
      for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        // calculate Jacobi matrix in the integration point
        MappedIntegrationPoint<D,D> mip(ir[i], eltrans);

        double param_normaldiff = coef_param_normaldiff -> Evaluate (mip);
        // alpha(x) (diffusion coefficient)
        double alpha = coef_alpha -> Evaluate (mip);

        Mat<D> inv_jac = mip.GetJacobianInverse();
        Vec<D> interface_normal = Trans(inv_jac) * interface_normal_ref ;
        interface_normal /= L2Norm(interface_normal);
        
        fel_l2.CalcMappedDShape (mip, dshape);
        dshapen = dshape * interface_normal;
        
        // integration weight and Jacobi determinant
        double fac = mip.IP().Weight() * mip.GetMeasure();

        const double h = D==2 ? sqrt(mip.GetMeasure()) : cbrt(mip.GetMeasure());
        
        elmat.Rows(l2_dofs).Cols(l2_dofs) += (fac/h*alpha*param_normaldiff) * dshapen * Trans(dshapen);
        
      }
    }

    {
      static Timer timere ("HDGTraceLBInt::interface part (Lapl.Beltr.)");
      RegionTimer reg (timere);

      FlatMatrixFixWidth<D> dshape(nd_l2,lh);
      FlatMatrixFixWidth<D> proj(nd_l2,lh);
      const FlatXLocalGeometryInformation & xgeom(tracel2fe->GetFlatLocalGeometry());
      const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
      const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

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
        const double alpha = coef_alpha->Evaluate(mip);

        fel_l2.CalcMappedDShape(mip, dshape);
        proj = dshape * (Id<D>() - normal * Trans(normal)) ;
        elmat.Rows(l2_dofs).Cols(l2_dofs) += (alpha * weight) * proj * Trans(proj);
      }

    }
    
  }

  template class HDGTraceLaplaceBeltramiIntegrator<2>;
  template class HDGTraceLaplaceBeltramiIntegrator<3>;

  static RegisterBilinearFormIntegrator<HDGTraceLaplaceBeltramiIntegrator<2> > inithdgtracelaplacebeltrami2d ("hdgtracelaplacebeltrami", 2, 1);
  static RegisterBilinearFormIntegrator<HDGTraceLaplaceBeltramiIntegrator<3> > inithdgtracelaplacebeltrami3d ("hdgtracelaplacebeltrami", 3, 1);

};
