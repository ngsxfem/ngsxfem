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
    shared_ptr<GridFunction> gf_lset_ho;
    shared_ptr<GridFunction> deform;
    bool no_cut_off;
    double threshold=0.1;
    double volume_ctrl=-1;

    bool dynamic_search_dir=false; //take normal of accurate level set 

    // statistics
    double * n_maxits;
    double * n_totalits;
    double * n_deformed_points;
    double * n_accepted_points;
    double * n_corrected_points;
  public:


    NumProcGeometryTest (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset", ""), false);
      gf_lset_p1  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p1", ""), false);
      gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_ho", ""), true);
      if (!gf_lset_ho)
        gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p2", ""), false);
      deform  = apde->GetGridFunction (flags.GetStringFlag ("deformation", ""));
      no_cut_off = flags.GetDefineFlag("nocutoff");
      
      dynamic_search_dir = flags.GetDefineFlag("dynamic_search_dir");
      
      threshold = flags.GetNumFlag("threshold",0.1);
      
      volume_ctrl = flags.GetNumFlag("volume",-1);

      //statistic variables:
      apde->AddVariable ("npgeomtest.maxits", 0.0, 6);
      n_maxits = & apde->GetVariable ("npgeomtest.maxits");
      apde->AddVariable ("npgeomtest.totalits", 0.0, 6);
      n_totalits = & apde->GetVariable ("npgeomtest.totalits");
      apde->AddVariable ("npgeomtest.deformed_points", 0.0, 6);
      n_deformed_points = & apde->GetVariable ("npgeomtest.deformed_points");
      apde->AddVariable ("npgeomtest.accepted_points", 0.0, 6);
      n_accepted_points = & apde->GetVariable ("npgeomtest.accepted_points");
      apde->AddVariable ("npgeomtest.corrected_points", 0.0, 6);
      n_corrected_points = & apde->GetVariable ("npgeomtest.corrected_points");
    }

    virtual ~NumProcGeometryTest()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcGeometryTest";
    }

    void SearchCorrespondingPoint (
      const ScalarFiniteElement<D> & sca_fe, FlatVector<> sca_values,
      const Vec<D> & init_point, double goal_val, const Vec<D> & init_search_dir,
      Vec<D> & final_point, LocalHeap & lh)
    {
      HeapReset hr(lh);
      
      IntegrationPoint curr_ip;
      for (int d = 0; d < D; ++d) curr_ip(d) = init_point(d);

      Vec<D> search_dir = init_search_dir;

      FlatVector<> shape(sca_fe.GetNDof(),lh);
      FlatMatrixFixWidth<D> dshape(sca_fe.GetNDof(),lh);

      int it = 0;
      for (it = 0; it < 100; ++it)
      {
        sca_fe.CalcShape(curr_ip, shape);
        sca_fe.CalcDShape(curr_ip, dshape);

        const double curr_val = InnerProduct(shape, sca_values);
        const Vec<D> curr_grad = Trans(dshape) * sca_values;

        const double curr_defect = goal_val - curr_val;
        if (abs(curr_defect) < 1e-12)
          break;

        if (dynamic_search_dir)
        {
          search_dir = curr_grad;
          search_dir /= L2Norm(search_dir);
        }
        
        const double dphidn = InnerProduct(curr_grad,search_dir);

        for (int d = 0; d < D; ++d)
          curr_ip(d) += curr_defect / dphidn * search_dir(d);
      }

#pragma omp critical (totalits)
      *n_totalits += it;
#pragma omp critical (maxits)
      *n_maxits = max((double)it,*n_maxits);

      if (it == 100)
        final_point = init_point;
      else
        for (int d = 0; d < D; ++d) final_point(d) = curr_ip(d);
    }


    

    void Test (LocalHeap & lh)
    {
      double volume = 0.0;
      
      int ne=ma->GetNE();

      for (int elnr = 0; elnr < ne; ++elnr)
      {
        HeapReset hr(lh);
        Ngs_Element ngel = ma->GetElement(elnr);
        ELEMENT_TYPE eltype = ngel.GetType();
        Array<int> dofs;
        gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dofs);
        FlatVector<> vals(dofs.Size(),lh);
        gf_lset_p1->GetVector().GetIndirect(dofs,vals);

        bool cut_els = false;

        if (vals[0] == 0.0) cut_els = true;
        
        DOMAIN_TYPE first_dt = vals[0] > 0 ? POS : NEG;
        for (int i = 1; i < 3; ++i)
          if (first_dt == POS)
          {
            if (vals[i] <= 0.0)
              cut_els = true;
          }
          else
          {
            if (vals[i] >= 0.0)
              cut_els = true;
          }

        if (!cut_els)
        {
          if (first_dt==NEG)
          {
            ma->SetDeformation(deform);
            ElementTransformation & eltrans_curved = ma->GetTrafo (ElementId(VOL,elnr), lh);
            
            IntegrationRule pir = SelectIntegrationRule (eltrans_curved.GetElementType(), 4);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans_curved);
              volume += mip.GetWeight();
            }

          }
          continue;
        }
        ma->SetDeformation(nullptr);
        ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);

        ScalarFieldEvaluator * lset_eval_p
          = ScalarFieldEvaluator::Create(D,*gf_lset_ho,eltrans,lh);

        auto cquad = new CompositeQuadratureRule<D>() ;

        ELEMENT_TYPE et_time = ET_POINT;

        auto xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                       *cquad, lh, 
                                                       4, 0, 
                                                       0, 0);

        xgeom->MakeQuadRule();
        FlatXLocalGeometryInformation fxgeom(*xgeom,lh);
        const FlatCompositeQuadratureRule<D> & fcompr(fxgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(NEG));

        ma->SetDeformation(deform);
        ElementTransformation & eltrans_curved = ma->GetTrafo (ElementId(VOL,elnr), lh);
        
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans_curved);
          volume += mip.GetWeight();
        }
        
        ma->SetDeformation(nullptr);

      }
      cout << " volume = " << volume << endl;
      if (volume_ctrl>0)
        cout << " volume error = " << abs(volume-volume_ctrl) << endl;
      ma->SetDeformation(deform);
    }
    
    //only setting P1 interpolant for now..
    void SetInterpolants(LocalHeap & clh)
    {
      int ne=ma->GetNE();
      gf_lset_p1->GetVector() = 0.0;
      // gf_lset_ho->GetVector() = 0.0;

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
            gf_lset_ho->GetFESpace()->GetVertexDofNrs(vnums[i],dof);
            // gf_lset_ho->GetVector().SetIndirect(dof,val);
            
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
            
            gf_lset_ho->GetFESpace()->GetEdgeDofNrs(facet,dof);
            FlatVector<> val(1,&edge_value);
            // gf_lset_ho->GetVector().SetIndirect(dof,val);
              // deform->GetVector().SetIndirect(dnums,values);
            
          }
          
        }
      }

    }
    
    
    virtual void Do (LocalHeap & clh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;

      shared_ptr<BaseVector> factor = deform->GetVector().CreateVector();
      *factor = 0.0;
      deform->GetVector() = 0.0;
      
      int ne=ma->GetNE();

      shared_ptr<FESpace> fes_deform = deform->GetFESpace();

      SetInterpolants(clh);
      
      {
        auto gf_fes_ho = gf_lset_ho->GetFESpace();
        auto gf_fes_p1 = gf_lset_p1->GetFESpace();
        LocalHeap lh(clh.Split());
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(lh);
          Ngs_Element ngel = ma->GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();
          if (eltype!= ET_TRIG)
            throw Exception("only trigs for now..");

          const FiniteElement & fe_ho = gf_fes_ho->GetFE(elnr,lh);
          const ScalarFiniteElement<D>& sca_fe_ho = dynamic_cast<const ScalarFiniteElement<D>&>(fe_ho);
          Array<int> dnums_lset_ho;
          gf_fes_ho->GetDofNrs(elnr,dnums_lset_ho);
          FlatVector<> lset_vals_ho(dnums_lset_ho.Size(),lh);
          gf_lset_ho->GetVector().GetIndirect(dnums_lset_ho,lset_vals_ho);
          
          const FiniteElement & fe_p1 = gf_fes_p1->GetFE(elnr,lh);
          const ScalarFiniteElement<D>& sca_fe_p1 = dynamic_cast<const ScalarFiniteElement<D>&>(fe_p1);
          Array<int> dnums_lset_p1;
          FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),lh);
          gf_fes_p1->GetDofNrs(elnr,dnums_lset_p1);
          gf_lset_p1->GetVector().GetIndirect(dnums_lset_p1,lset_vals_p1);
          
          FlatVector<> shape(dnums_lset_ho.Size(),lh);
          FlatMatrixFixWidth<D> dshape(dnums_lset_ho.Size(),lh);
          
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

          IntegrationPoint ip4(0.0,0.0,1.0);
          MappedIntegrationPoint<D,D> mip4(ip4,eltrans);
          double lset4 = lset->Evaluate(mip4);

          
          Vec<D> grad;
          grad(0) = lset2 - lset1;
          grad(1) = lset3 - lset1;
          if (D==3)
            grad(2) = lset4 - lset1;
          
          Vec<D> normal = grad;
          double len = L2Norm(normal);
          normal /= len;

          const double h = pow(mip1.GetJacobiDet(),1.0/D);

          if ( (abs(lset1/h) > 1.0) && (abs(lset2/h) > 1.0) && (abs(lset3/h) > 1.0))
            if (!no_cut_off)
              continue;
          
          
          ips.Append(ip2);
          ips.Append(ip3);
          if (D==3)
            ips.Append(ip4);
          ips.Append(ip1);

          
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
            
            IntegrationPoint curr_ip(0.0,0.0);
            curr_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            curr_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);

            IntegrationPoint old_ip(0.0,0.0);
            old_ip(0) = 0.5 * ips[v1](0) + 0.5 * ips[v2](0);
            old_ip(1) = 0.5 * ips[v1](1) + 0.5 * ips[v2](1);
            Vec<D> old_ref_point;
            old_ref_point(0) = old_ip(0);
            old_ref_point(1) = old_ip(1);

            MappedIntegrationPoint<D,D> mip_old (curr_ip, eltrans);
            

            // cout << " len = " << len << endl;
            
            double lset_lin = eval_linear(curr_ip);

            if (abs(lset_lin) > h && !no_cut_off ) break;
            // cout << " curr_ip = " << curr_ip << endl;

            *n_deformed_points+=1.0;
            
            Vec<D> old_coord = mip_old.GetPoint();

            Vec<D> orig_point;
            for (int d = 0; d < D; ++d) orig_point(d) = curr_ip(d);
            Vec<D> final_point;
            
            SearchCorrespondingPoint(sca_fe_ho, lset_vals_ho, orig_point,
                                     lset_lin, normal, final_point, lh);
            for (int d = 0; d < D; ++d) curr_ip(d) = final_point(d);
            
            Vec<D> dist = final_point - orig_point;

            double distnorm = L2Norm(dist);
            if (distnorm > threshold)
            {
              dist *= threshold / distnorm; 
              curr_ip(0) = old_ref_point(0) + dist(0);
              curr_ip(1) = old_ref_point(1) + dist(1);
              *n_corrected_points+=1.0;
            }
            else
            {
              *n_accepted_points+=1.0;
            }

            
            MappedIntegrationPoint<D,D> mip_new (curr_ip, eltrans);
            Vec<D> new_coord = mip_new.GetPoint();
            Vec<D> deform_vec;
            FlatVector<> values(D,&deform_vec(0));
            fes_deform->GetEdgeDofNrs(facet, dnums);
            deform->GetVector().GetIndirect(dnums,values);
            deform_vec += -8.0 * (new_coord - old_coord);
            // deform_vec *= -8.0;
            
            deform->GetVector().SetIndirect(dnums,values);

            factor->GetIndirect(dnums,values);
            values(0) += 1.0;
            factor->SetIndirect(dnums,values);
            
          }
        }

        Array<int> dnums(1);
        for (int i = 0; i < factor->Size(); ++i)
        {
          FlatVector<> val_fac(D,lh);
          FlatVector<> values(D,lh);
          dnums[0] = i;
          deform->GetVector().GetIndirect(dnums,values);
          factor->GetIndirect(dnums,val_fac);
          if (val_fac(0) > 0)
            values *= 1.0/val_fac(0);
          deform->GetVector().SetIndirect(dnums,values);
        }
        
        cout << " deformpoints = " << *n_deformed_points << endl;
        cout << " totalits = " << *n_totalits << endl;
        cout << " totalits/deformpoints = " << *n_totalits/ *n_deformed_points << endl;
        cout << " accepted_points = " << *n_accepted_points << endl;
        cout << " corrected_points = " << *n_corrected_points << endl;
        cout << " maxits = " << *n_maxits << endl;
      }

      ma->SetDeformation(deform);

      Test(clh);
    }    
    

  };




  class NumProcUnsetDeformation : public NumProc
  {
  public:
    shared_ptr<MeshAccess> ma;
    NumProcUnsetDeformation (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      ma = apde->GetMeshAccess();
    }

    virtual ~NumProcUnsetDeformation()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcUnsetDeformation";
    }


    virtual void Do (LocalHeap & clh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;
      ma->SetDeformation(nullptr);
      
    }
  };
    
}
static RegisterNumProc<NumProcUnsetDeformation > npxgeomtest("unsetdeformation");
static RegisterNumProc<NumProcGeometryTest<2> > npxgeomtestasdf("xgeomtest");
