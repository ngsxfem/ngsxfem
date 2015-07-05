/*********************************************************************/
/* File:   geometrytests.cpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   5. Jul. 2015                                              */
/*********************************************************************/

///HACKED.... INCLUDES are done outside

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

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
    shared_ptr<GridFunction> deform;
  public:


    NumProcGeometryTest (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset", ""));
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
      

// #pragma omp parallel
      for (int it = 0; it < 10; ++it)
      {
        LocalHeap lh(clh.Split());
// #pragma omp for schedule(static)
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

          Array<int> dnums;
          for (int f = 1; f < 3; ++f)
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

            cout << " v1 = " << v1 << endl;
            cout << " v2 = " << v2 << endl;
            
            IntegrationPoint curr_ip(0.0,0.0);
            curr_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            curr_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);
            
            MappedIntegrationPoint<D,D> curr_mip(curr_ip,eltrans);
            cout << " curr_mip = " << curr_mip << endl;

            double lset_prec = lset->Evaluate(curr_mip); // LATER: consider deform
            double lset_lin = eval_linear(curr_ip);

            cout << " curr_mip.GetPoint() = " << curr_mip.GetPoint() << endl;
            
            cout << " lset_prec = " << lset_prec << endl;
            cout << " lset_lin = " << lset_lin << endl;

            double f0 = lset_prec - lset_lin;

            // cout << " f0 = " << f0 << endl;
            
            Vec<D> grad;
            CalcGradientOfCoeff<D>(lset, curr_mip, grad, lh);

            double len = L2Norm(grad);
            Vec<D> normal = grad;
            normal /= len;

            // cout << " normal = " << normal << endl;
            
            Vec<D> update;
            FlatVector<> values(D,&update(0));
            // cout << " values = " << values << endl;
            fes_deform->GetEdgeDofNrs(facet, dnums);
            // cout << " dnums = " << dnums << endl;
            deform->GetVector().GetIndirect(dnums,values);
            // cout << " values = " << values << endl;
            
            // cout << " update = " << update << endl;

            update += 10.0 * f0 / len * normal;
            
            // cout << " update = " << update << endl;


            deform->GetVector().SetIndirect(dnums,values);

            // getchar();
            // break;
          }
          // break;
        }
        ma->SetDeformation(deform);

      }
      
    }    
    

  };
}
static RegisterNumProc<NumProcGeometryTest<2> > npxgeomtest("xgeomtest");
