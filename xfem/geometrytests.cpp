/*********************************************************************/
/* File:   geometrytests.cpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   5. Jul. 2015                                              */
/*********************************************************************/

///HACKED.... INCLUDES are done outside

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;
#include <typeinfo>  //for 'typeid' to work  

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
namespace ngcomp
{ 
  template<int D>
  void CalcGradientOfCoeff(shared_ptr<CoefficientFunction> coef, const MappedIntegrationPoint<D,D>& mip,
                           Vec<D>& der, LocalHeap& lh)
  {
    HeapReset hr(lh);
    // bmatu = 0;
    // evaluate dshape by numerical diff
    //fel_u, eltrans, sip, returnval, lh

    const IntegrationPoint& ip = mip.IP();//volume_ir[i];
    const ElementTransformation & eltrans = mip.GetTransformation();
    
    Vec<D> der_ref;
  
    double eps = 1e-7;
    for (int j = 0; j < D; j++)   // d / dxj
    {
      IntegrationPoint ipl(ip);
      ipl(j) -= eps;
      MappedIntegrationPoint<D,D> mipl(ipl, eltrans);

      IntegrationPoint ipr(ip);
      ipr(j) += eps;
      MappedIntegrationPoint<D,D> mipr(ipr, eltrans);

      const double valright = coef->Evaluate(mipr);
      const double valleft = coef->Evaluate(mipl);
      
      der_ref[j] = (1.0/(2*eps)) * (valright-valleft);
    }
    der = Trans(mip.GetJacobianInverse()) * der_ref;
  }

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcGeometryTest : public NumProc
  {
  protected:
    shared_ptr<CoefficientFunction> lset;
    shared_ptr<GridFunction> gf_lset_p1;
    shared_ptr<GridFunction> gf_lset_p2;
    shared_ptr<GridFunction> deform;
  public:


    NumProcGeometryTest (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset", ""), true);
      gf_lset_p1  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p1", ""), true);
      gf_lset_p2  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p2", ""), true);
      deform  = apde->GetGridFunction (flags.GetStringFlag ("deformation", ""));
    }

    virtual ~NumProcGeometryTest()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcGeometryTest";
    }


    virtual void Do (LocalHeap & clh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;

      int ne=ma->GetNE();
      int nedges=ma->GetNEdges();
      int nf=ma->GetNFaces();
      int nv=ma->GetNV();
      int nse=ma->GetNSE();

      shared_ptr<FESpace> fes_deform = deform->GetFESpace();

      gf_lset_p2->GetVector() = 0.0;
      
      if (gf_lset_p1 != nullptr)
      {
        LocalHeap lh(clh.Split());
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(lh);
          Ngs_Element ngel = ma->GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();
          if (eltype!= ET_TRIG)
            throw Exception("only trigs for now..");

          ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
          
          IntegrationPoint ip1(0.0,0.0);
          MappedIntegrationPoint<D,D> mip1(ip1,eltrans);
          double lset1 = lset->Evaluate(mip1);
          
          IntegrationPoint ip2(1.0,0.0);
          MappedIntegrationPoint<D,D> mip2(ip2,eltrans);
          double lset2 = lset->Evaluate(mip2);

          IntegrationPoint ip3(0.0,1.0);
          MappedIntegrationPoint<D,D> mip3(ip3,eltrans);
          double lset3 = lset->Evaluate(mip3);

          Array<IntegrationPoint> ips(0);
          ips.Append(ip2);
          ips.Append(ip3);
          ips.Append(ip1);
          
          Array<double> lsetvals(0);
          lsetvals.Append(lset2);
          lsetvals.Append(lset3);
          lsetvals.Append(lset1);
          
          Array<int> vnums;
          Array<int> dof;
          ma->GetElPNums(elnr, vnums);
          for (int i = 0; i < 3; ++i)
          {
            gf_lset_p1->GetFESpace()->GetVertexDofNrs(vnums[i],dof);
            FlatVector<> val(1,&lsetvals[i]);
            gf_lset_p1->GetVector().SetIndirect(dof,val);
              // deform->GetVector().SetIndirect(dnums,values);
            gf_lset_p2->GetFESpace()->GetVertexDofNrs(vnums[i],dof);
            gf_lset_p2->GetVector().SetIndirect(dof,val);
            
          }


          Array<int> facets;
          Array<int> verts;
          Array<int> facetverts;
          ma->GetElEdges(elnr,facets);
          ma->GetElVertices(elnr,verts);
          
          Array<int> dnums;
          for (int f = 0; f < 3; ++f)
          {
            int facet = facets[f];
            ma->GetFacetPNums(facet,facetverts);
            
            int v1 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[0])
                v1 = i;
            int v2 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[1])
                v2 = i;

            // cout << " v1 = " << v1 << endl;
            // cout << " v2 = " << v2 << endl;
            
            IntegrationPoint curr_ip(0.0,0.0);
            curr_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            curr_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);

            MappedIntegrationPoint<D,D> curr_mip(curr_ip,eltrans);
            double lset_edge = lset->Evaluate(curr_mip);

            double edge_value = lset_edge - 0.5 * lsetvals[v1] - 0.5 * lsetvals[v2];
            edge_value *= -8.0;
            
            gf_lset_p2->GetFESpace()->GetEdgeDofNrs(facet,dof);
            FlatVector<> val(1,&edge_value);
            gf_lset_p2->GetVector().SetIndirect(dof,val);
              // deform->GetVector().SetIndirect(dnums,values);
            
          }
          
        }
      }
      
      // cout << " typeid(lset) = " << typeid(lset).name() << endl;
      
      // auto gf_lset = dynamic_pointer_cast<GridFunction>(lset);

      // cout << " gf_lset = " << gf_lset << endl;
      
      if (gf_lset_p2 != nullptr)
      {
        auto gf_fes_p2 = gf_lset_p2->GetFESpace();
        LocalHeap lh(clh.Split());
        int totalits = 0;
        int deformpoints = 0;
        int corrected_points = 0;
        int accepted_points = 0;
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(lh);
          Ngs_Element ngel = ma->GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();
          if (eltype!= ET_TRIG)
            throw Exception("only trigs for now..");

          const FiniteElement & fe = gf_fes_p2->GetFE(elnr,lh);
          const ScalarFiniteElement<D>& sca_fe = dynamic_cast<const ScalarFiniteElement<D>&>(fe);
          Array<int> sca_dnums;
          gf_fes_p2->GetDofNrs(elnr,sca_dnums);
          FlatVector<> sca_values(sca_dnums.Size(),lh);
          FlatVector<> shape(sca_dnums.Size(),lh);
          FlatMatrixFixWidth<D> dshape(sca_dnums.Size(),lh);
          gf_lset_p2->GetVector().GetIndirect(sca_dnums,sca_values);
          
          
          ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
        
          Array<IntegrationPoint> ips(0);
          IntegrationPoint ip1(0.0,0.0);
          MappedIntegrationPoint<D,D> mip1(ip1,eltrans);
          double lset1 = lset->Evaluate(mip1);
          
          IntegrationPoint ip2(1.0,0.0);
          MappedIntegrationPoint<D,D> mip2(ip2,eltrans);
          double lset2 = lset->Evaluate(mip2);

          IntegrationPoint ip3(0.0,1.0);
          MappedIntegrationPoint<D,D> mip3(ip3,eltrans);
          double lset3 = lset->Evaluate(mip3);

          // cout << " lset1 = " << lset1 << endl;
          // cout << " lset2 = " << lset2 << endl;
          // cout << " lset3 = " << lset3 << endl;

          const double h = pow(mip1.GetJacobiDet(),1.0/D);

          if ( (abs(lset1/h) > 1.0) && (abs(lset2/h) > 1.0) && (abs(lset3/h) > 1.0))
            continue;
          
          
          ips.Append(ip2);
          ips.Append(ip3);
          ips.Append(ip1);

          // IntegrationPoint ip4(0.5,0.0);
          // ips.Append(ip4);
          // MappedIntegrationPoint<D,D> mip4(ip4,eltrans);
          // double lset4 = lset->Evaluate(mip4);
          
          auto eval_linear = [lset1,lset2,lset3] (const IntegrationPoint & ip)
            { return (1-ip(0)-ip(1))*lset1+ip(0)*lset2+ip(1)*lset3; };

          // IntegrationPoints
          
          Array<int> facets;
          Array<int> verts;
          Array<int> facetverts;
          ma->GetElEdges(elnr,facets);
          ma->GetElVertices(elnr,verts);

          Array<int> dnums;
          for (int f = 0; f < 3; ++f)
          {
            int facet = facets[f];
            ma->GetFacetPNums(facet,facetverts);
            
            int v1 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[0])
                v1 = i;
            int v2 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[1])
                v2 = i;

            // cout << " v1 = " << v1 << endl;
            // cout << " v2 = " << v2 << endl;
            
            IntegrationPoint curr_ip(0.0,0.0);
            curr_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            curr_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);

            IntegrationPoint old_ip(0.0,0.0);
            old_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            old_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);
            Vec<D> old_ref_point;
            old_ref_point(0) = old_ip(0);
            old_ref_point(1) = old_ip(1);
            // MappedIntegrationPoint<D,D> first_point(ips[v1],eltrans);
            // MappedIntegrationPoint<D,D> second_point(ips[v2],eltrans);
            // Vec<D> tang = second_point.GetPoint() - first_point.GetPoint();
            // Vec<D> normal; normal(0) = -tang(1); normal(1) = tang(0);
            // normal /= L2Norm(normal);
            
            Vec<D> grad;
            grad(0) = lset2 - lset1;
            grad(1) = lset3 - lset1;

            // cout << " grad = " << grad << endl;

            MappedIntegrationPoint<D,D> mip_old (curr_ip, eltrans);
            
            // Vec<D> normal = grad;
            // double len = L2Norm(normal);

            // normal /= len;

            // cout << " len = " << len << endl;
            
            double lset_lin = eval_linear(curr_ip);

            if (abs(lset_lin) > h ) break;
            // cout << " curr_ip = " << curr_ip << endl;

            deformpoints++;
            
            Vec<D> old_coord = mip_old.GetPoint();

            Vec<D> dist;
            dist = 0.0;

            for (int it = 0; it < 100; ++it)
            {
              // ElementTransformation & new_eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
              sca_fe.CalcShape(curr_ip, shape);
              double lset_prec = InnerProduct(shape, sca_values); // LATER: consider deform
            
              // MappedIntegrationPoint<D,D> curr_mip(curr_ip,new_eltrans);
              // cout << " curr_mip = " << curr_mip << endl;


              // cout << " curr_mip.GetPoint() = " << curr_mip.GetPoint() << endl;
            
              // cout << " lset_prec = " << lset_prec << endl;
              // cout << " lset_lin = " << lset_lin << endl;
              // if (abs(lset_prec) > 0.75*h)
              //   break;
              
              double f0 = lset_lin - lset_prec;

              // cout << " f0 = " << f0 << endl;

              if (abs(f0) < 1e-12)
                break;
              
              totalits++;
            
              sca_fe.CalcDShape(curr_ip, dshape);

              // cout << " sca_values = " << sca_values << endl;
              
              Vec<D> grad = Trans(dshape) * sca_values; 


              Vec<D> normal = grad;
              double len = L2Norm(normal);

              normal /= len;
              
              // cout << " grad = " << grad << endl;
              
              const double dphidn = InnerProduct(grad,normal);

              // cout << " dphidn = " << dphidn << endl;
              // cout << " normal = " << normal << endl;
            
              Vec<D> update;
              // FlatVector<> values(D,&update(0));
              // // cout << " values = " << values << endl;
              // fes_deform->GetEdgeDofNrs(facet, dnums);
              // // cout << " dnums = " << dnums << endl;
              // deform->GetVector().GetIndirect(dnums,values);
              // cout << " values = " << values << endl;
            
              // cout << " update = " << update << endl;

              update = f0 / dphidn * normal;

              
              // cout << " update = " << update << endl;

              curr_ip(0) += update(0);
              curr_ip(1) += update(1);

              dist += update;
              // cout << " curr_ip = " << curr_ip << endl;
              
              // deform->GetVector().SetIndirect(dnums,values);
              // getchar();
            }

            double distnorm = L2Norm(dist);
            if (distnorm > 0.1)
            {
              dist *= 0.1 / distnorm; 
              curr_ip(0) = old_ref_point(0) + dist(0);
              curr_ip(1) = old_ref_point(1) + dist(1);
              corrected_points++;
            }
            else
            {
              accepted_points++;
            }
            
            MappedIntegrationPoint<D,D> mip_new (curr_ip, eltrans);
            Vec<D> new_coord = mip_new.GetPoint();
            // Vec<D> deform_vec = new_coord - old_coord;
            Vec<D> deform_vec = new_coord - old_coord;
            deform_vec *= -8.0;

            // cout << " old_coord = " << old_coord << endl;
            // cout << " new_coord = " << new_coord << endl;
            
            FlatVector<> values(D,&deform_vec(0));
            fes_deform->GetEdgeDofNrs(facet, dnums);
            deform->GetVector().SetIndirect(dnums,values);
            
            // std::cout << " -- " << std::endl;

            // getchar();
            // break;
          }
          // break;
        }
        cout << " deformpoints = " << deformpoints << endl;
        cout << " totalits = " << totalits << endl;
        cout << " totalits/deformpoints = " << totalits/deformpoints << endl;
        cout << " accepted_points = " << accepted_points << endl;
        cout << " corrected_points = " << corrected_points << endl;
      }
      else
      {
        ma->SetDeformation(deform);
        LocalHeap lh(clh.Split());
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(lh);
          Ngs_Element ngel = ma->GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();
          if (eltype!= ET_TRIG)
            throw Exception("only trigs for now..");
          
          ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
        
          Array<IntegrationPoint> ips(0);
          IntegrationPoint ip1(0.0,0.0);
          MappedIntegrationPoint<D,D> mip1(ip1,eltrans);
          double lset1 = lset->Evaluate(mip1);
          
          IntegrationPoint ip2(1.0,0.0);
          MappedIntegrationPoint<D,D> mip2(ip2,eltrans);
          double lset2 = lset->Evaluate(mip2);

          IntegrationPoint ip3(0.0,1.0);
          MappedIntegrationPoint<D,D> mip3(ip3,eltrans);
          double lset3 = lset->Evaluate(mip3);

          ips.Append(ip2);
          ips.Append(ip3);
          ips.Append(ip1);

          // IntegrationPoint ip4(0.5,0.0);
          // ips.Append(ip4);
          // MappedIntegrationPoint<D,D> mip4(ip4,eltrans);
          // double lset4 = lset->Evaluate(mip4);
          
          auto eval_linear = [lset1,lset2,lset3] (const IntegrationPoint & ip)
            { return (1-ip(0)-ip(1))*lset1+ip(0)*lset2+ip(1)*lset3; };

          // IntegrationPoints
          
          Array<int> facets;
          Array<int> verts;
          Array<int> facetverts;
          ma->GetElEdges(elnr,facets);
          ma->GetElVertices(elnr,verts);

          const double h = pow(mip1.GetJacobiDet(),1.0/D);
          
          Array<int> dnums;
          for (int f = 0; f < 3; ++f)
          {
            int facet = facets[f];
            ma->GetFacetPNums(facet,facetverts);
            
            int v1 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[0])
                v1 = i;
            int v2 = -1;
            for (int i = 0; i < verts.Size(); i++)
              if (verts[i] == facetverts[1])
                v2 = i;

            // cout << " v1 = " << v1 << endl;
            // cout << " v2 = " << v2 << endl;
            
            IntegrationPoint curr_ip(0.0,0.0);
            curr_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            curr_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);

            // MappedIntegrationPoint<D,D> first_point(ips[v1],eltrans);
            // MappedIntegrationPoint<D,D> second_point(ips[v2],eltrans);
            // Vec<D> tang = second_point.GetPoint() - first_point.GetPoint();
            // Vec<D> normal; normal(0) = -tang(1); normal(1) = tang(0);
            // normal /= L2Norm(normal);

            for (int it = 0; it < 100; ++it)
            {
              ElementTransformation & new_eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);

            
              MappedIntegrationPoint<D,D> curr_mip(curr_ip,new_eltrans);
              // cout << " curr_mip = " << curr_mip << endl;

              double lset_prec = lset->Evaluate(curr_mip); // LATER: consider deform
              double lset_lin = eval_linear(curr_ip);

              // cout << " curr_mip.GetPoint() = " << curr_mip.GetPoint() << endl;
            
              // cout << " lset_prec = " << lset_prec << endl;
              // cout << " lset_lin = " << lset_lin << endl;
              if (abs(lset_prec) > 0.75*h)
                break;
              
              double f0 = lset_prec - lset_lin;

              cout << " f0 = " << f0 << endl;

              if (abs(f0) < 1e-6)
                break;
            
              Vec<D> grad;
              CalcGradientOfCoeff<D>(lset, curr_mip, grad, lh);

              double len = L2Norm(grad);
              Vec<D> normal = grad;

              normal /= len;

              const double dphidn = InnerProduct(grad,normal);
              // cout << " normal = " << normal << endl;
            
              Vec<D> update;
              FlatVector<> values(D,&update(0));
              // cout << " values = " << values << endl;
              fes_deform->GetEdgeDofNrs(facet, dnums);
              // cout << " dnums = " << dnums << endl;
              deform->GetVector().GetIndirect(dnums,values);
              // cout << " values = " << values << endl;
            
              // cout << " update = " << update << endl;

              update += f0 / dphidn * normal;
            
              // cout << " update = " << update << endl;


              deform->GetVector().SetIndirect(dnums,values);
            }
            std::cout << " -- " << std::endl;

            // getchar();
            // break;
          }
          // break;
        }

      }

      ma->SetDeformation(nullptr);
    }    
    

  };
}
static RegisterNumProc<NumProcGeometryTest<2> > npxgeomtest("xgeomtest");
