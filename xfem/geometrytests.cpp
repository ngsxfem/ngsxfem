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
    double accept_threshold=0.1;
    double reject_threshold=0.4;
    double volume_ctrl=-1;

    bool dynamic_search_dir=false; //take normal of accurate level set 

    double lower_lset_bound=0.0; //domain of interest for deformation: lower bound
    double upper_lset_bound=0.0; //domain of interest for deformation: lower bound

    bool no_edges = false;
    bool no_faces = false;
    bool no_cells = false;
    
    // statistics
    double * n_maxits;
    double * n_totalits;

    double * n_deformed_points_edge;
    double * n_accepted_points_edge;
    double * n_rejected_points_edge;
    double * n_corrected_points_edge;

    double * n_deformed_points_face;
    double * n_accepted_points_face;
    double * n_rejected_points_face;
    double * n_corrected_points_face;

    double * n_deformed_points_cell;
    double * n_accepted_points_cell;
    double * n_rejected_points_cell;
    double * n_corrected_points_cell;
  public:


    NumProcGeometryTest (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset", ""), true);
      gf_lset_p1  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p1", ""), false);
      gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_ho", ""), true);
      if (!gf_lset_ho)
        gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p2", ""), false);
      deform  = apde->GetGridFunction (flags.GetStringFlag ("deformation", ""));

      no_edges = flags.GetDefineFlag("no_edges");
      no_faces = flags.GetDefineFlag("no_faces");
      no_cells = flags.GetDefineFlag("no_cells");
      
      lower_lset_bound = flags.GetNumFlag("lower_lset_bound",0.0);
      upper_lset_bound = flags.GetNumFlag("upper_lset_bound",0.0);
      
      dynamic_search_dir = flags.GetDefineFlag("dynamic_search_dir");
      
      accept_threshold = flags.GetNumFlag("accept_threshold",flags.GetNumFlag("threshold",0.1));
      reject_threshold = flags.GetNumFlag("reject_threshold",accept_threshold * 4.0);
      
      volume_ctrl = flags.GetNumFlag("volume",-1);

      //statistic variables:
      apde->AddVariable ("npgeomtest.maxits", 0.0, 6);
      n_maxits = & apde->GetVariable ("npgeomtest.maxits");
      apde->AddVariable ("npgeomtest.totalits", 0.0, 6);
      n_totalits = & apde->GetVariable ("npgeomtest.totalits");

      apde->AddVariable ("npgeomtest.deformed_points_edge", 0.0, 6);
      n_deformed_points_edge = & apde->GetVariable ("npgeomtest.deformed_points_edge");
      apde->AddVariable ("npgeomtest.accepted_points_edge", 0.0, 6);
      n_accepted_points_edge = & apde->GetVariable ("npgeomtest.accepted_points_edge");
      apde->AddVariable ("npgeomtest.corrected_points_edge", 0.0, 6);
      n_corrected_points_edge = & apde->GetVariable ("npgeomtest.corrected_points_edge");
      apde->AddVariable ("npgeomtest.rejected_points_edge", 0.0, 6);
      n_rejected_points_edge = & apde->GetVariable ("npgeomtest.rejected_points_edge");

      apde->AddVariable ("npgeomtest.deformed_points_face", 0.0, 6);
      n_deformed_points_face = & apde->GetVariable ("npgeomtest.deformed_points_face");
      apde->AddVariable ("npgeomtest.accepted_points_face", 0.0, 6);
      n_accepted_points_face = & apde->GetVariable ("npgeomtest.accepted_points_face");
      apde->AddVariable ("npgeomtest.corrected_points_face", 0.0, 6);
      n_corrected_points_face = & apde->GetVariable ("npgeomtest.corrected_points_face");
      apde->AddVariable ("npgeomtest.rejected_points_face", 0.0, 6);
      n_rejected_points_face = & apde->GetVariable ("npgeomtest.rejected_points_face");

      apde->AddVariable ("npgeomtest.deformed_points_cell", 0.0, 6);
      n_deformed_points_cell = & apde->GetVariable ("npgeomtest.deformed_points_cell");
      apde->AddVariable ("npgeomtest.accepted_points_cell", 0.0, 6);
      n_accepted_points_cell = & apde->GetVariable ("npgeomtest.accepted_points_cell");
      apde->AddVariable ("npgeomtest.corrected_points_cell", 0.0, 6);
      n_corrected_points_cell = & apde->GetVariable ("npgeomtest.corrected_points_cell");
      apde->AddVariable ("npgeomtest.rejected_points_cell", 0.0, 6);
      n_rejected_points_cell = & apde->GetVariable ("npgeomtest.rejected_points_cell");
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
      const Vec<D> & init_point, double goal_val,
      const Vec<D> & init_search_dir, Vec<D> & final_point, LocalHeap & lh)
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
      int nv=ma->GetNV();
      gf_lset_p1->GetVector() = 0.0;
      
      LocalHeap lh(clh.Split());
      for (int vnr = 0; vnr < nv; ++vnr)
      {
        Vec<D> point;
        ma->GetPoint<D>(vnr,point);
        Mat<1,D> pointmat;
        pointmat.Row(0) = point;
        IntegrationPoint ip(0.0);
        FE_ElementTransformation<0,D> eltrans(ET_POINT,pointmat);
        MappedIntegrationPoint<0,D> mip(ip,eltrans);

        double val_lset;
        if (lset)
          val_lset = lset->Evaluate(mip);
        else
        {
          Array<int> dof;
          gf_lset_ho->GetFESpace()->GetVertexDofNrs(vnr, dof);
          FlatVector<> fval(1,&val_lset);
          gf_lset_ho->GetVector().GetIndirect(dof,fval);
        }
        
        Array<int> dof;
        gf_lset_p1->GetFESpace()->GetVertexDofNrs(vnr,dof);
        FlatVector<> val(1,&val_lset);
        gf_lset_p1->GetVector().SetIndirect(dof,val);
      }
      
      if (false)
      {
        gf_lset_p1->GetVector() = 0.0;
        // gf_lset_ho->GetVector() = 0.0;
        int ne=ma->GetNE();
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
          for (int i = 0; i < D+1; ++i)
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
          for (int f = 0; f < D+1; ++f)
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
      cout << " This is the Do-call on refinement level " << refinements++ << std::endl;

      * n_maxits = 0.0;
      * n_totalits = 0.0;

      * n_deformed_points_edge = 0.0;
      * n_accepted_points_edge = 0.0;
      * n_rejected_points_edge = 0.0;
      * n_corrected_points_edge = 0.0;

      * n_deformed_points_face = 0.0;
      * n_accepted_points_face = 0.0;
      * n_rejected_points_face = 0.0;
      * n_corrected_points_face = 0.0;

      * n_deformed_points_cell = 0.0;
      * n_accepted_points_cell = 0.0;
      * n_rejected_points_cell = 0.0;
      * n_corrected_points_cell = 0.0;
      
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
          if (eltype!= ET_TRIG && eltype!=ET_TET)
            throw Exception("only simplices for now..");

          const FiniteElement & fe_ho = gf_fes_ho->GetFE(elnr,lh);
          const ScalarFiniteElement<D>& sca_fe_ho = dynamic_cast<const ScalarFiniteElement<D>&>(fe_ho);
          Array<int> dnums_lset_ho;
          gf_fes_ho->GetDofNrs(elnr,dnums_lset_ho);
          FlatVector<> lset_vals_ho(dnums_lset_ho.Size(),lh);
          gf_lset_ho->GetVector().GetIndirect(dnums_lset_ho,lset_vals_ho);
          
          const FiniteElement & fe_p1 = gf_fes_p1->GetFE(elnr,lh);
          const ScalarFiniteElement<D>& sca_fe_p1 = dynamic_cast<const ScalarFiniteElement<D>&>(fe_p1);
          Array<int> dnums_lset_p1;
          gf_fes_p1->GetDofNrs(elnr,dnums_lset_p1);
          FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),lh);
          gf_lset_p1->GetVector().GetIndirect(dnums_lset_p1,lset_vals_p1);

          // element only predicts its own deformation if (linear) intersected
          {
            bool has_pos = (lset_vals_p1[D] > lower_lset_bound);
            bool has_neg = (lset_vals_p1[D] < upper_lset_bound);
            for (int d = 0; d < D; ++d)
            {
              if (lset_vals_p1[d] > lower_lset_bound)
                has_pos = true;
              if (lset_vals_p1[d] < upper_lset_bound)
                has_neg = true;
            }
            if (! (has_pos && has_neg) ) 
              continue;
          } 
          
          FlatVector<> shape_p1(dnums_lset_p1.Size(),lh);
          FlatMatrixFixWidth<D> dshape_p1(dnums_lset_p1.Size(),lh);
          
          ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
        
          IntegrationPoint ip_center(1.0/(D+1),1.0/(D+1),1.0/(D+1));
          MappedIntegrationPoint<D,D> mip_center(ip_center,eltrans);
          
          sca_fe_p1.CalcShape(ip_center,shape_p1);
          const double lset_center_p1 = InnerProduct(shape_p1,lset_vals_p1);
          const double h = pow(mip_center.GetJacobiDet(),1.0/D);
          

          sca_fe_p1.CalcDShape(ip_center,dshape_p1);
          Vec<D> grad = Trans(dshape_p1) * lset_vals_p1;

          Vec<D> normal =  grad;
          double len = L2Norm(normal);
          normal /= len;
          
          Array<int> edges;
          Array<int> faces;
          Array<int> verts;

          ma->GetElVertices(elnr,verts);
          ma->GetElEdges(elnr,edges);
          ma->GetElFaces(elnr,faces);

          Array<int> edge_verts;
          Array<int> face_verts;

          if (!no_edges)
          for (int ref_edge_nr = 0; ref_edge_nr < D+1; ++ref_edge_nr)
          {
            HeapReset hr(lh);
            const int edge_order = deform->GetFESpace()->GetOrder();
            if (edge_order <= 1) continue;
            
            const int global_edge_nr = edges[ref_edge_nr];
            if (deform->GetFESpace()->IsDirichletEdge(global_edge_nr))
              continue;
            const int v1 = ElementTopology::GetEdges(eltype)[ref_edge_nr][0];
            const int v2 = ElementTopology::GetEdges(eltype)[ref_edge_nr][1];

            ma->GetEdgePNums(global_edge_nr, edge_verts);

            H1HighOrderFE<ET_SEGM> & edge_fe = *(new (lh) H1HighOrderFE<ET_SEGM>(edge_order));
            edge_fe.SetVertexNumbers(edge_verts);

            FlatVector<> shape_edge_ho(edge_fe.GetNDof(),lh);
            
            const int edge_dofs_offset = 2;
            int inner_edge_dofs = edge_fe.GetNDof() - edge_dofs_offset;
            FlatVector<> shape_only_edge_ho(inner_edge_dofs,&shape_edge_ho(edge_dofs_offset));
            FlatMatrix<> edge_mass_mat(inner_edge_dofs,inner_edge_dofs,lh);
            FlatMatrixFixWidth<D> edge_rhs_mat(inner_edge_dofs,lh);
            edge_mass_mat = 0.0;
            edge_rhs_mat = 0.0;
            
            const IntegrationRule & ir_edge = SelectIntegrationRule (ET_SEGM, 2*edge_order);

            IntegrationPoint curr_vol_ip;
            FlatMatrixFixWidth<D> deform_contribution(inner_edge_dofs,lh);
            
            for (int l = 0; l < ir_edge.GetNIP(); l++)
            {
              const IntegrationPoint & curr_ip_edge(ir_edge[l]); //one integration point case...
              
              for (int d = 0; d < D; ++d)
                curr_vol_ip(d) = curr_ip_edge(0) * ElementTopology::GetVertices(eltype)[v1][d]
                  + (1.0-curr_ip_edge(0)) * ElementTopology::GetVertices(eltype)[v2][d];

              IntegrationPoint old_ip(curr_vol_ip);
              Vec<D> old_ref_point = old_ip.Point();

              sca_fe_p1.CalcShape(curr_vol_ip,shape_p1);
              double lset_lin = InnerProduct(shape_p1,lset_vals_p1);
            
              Vec<D> orig_point;
              for (int d = 0; d < D; ++d) orig_point(d) = curr_vol_ip(d);
              Vec<D> final_point;
            
              //statistics:
              *n_deformed_points_edge+=1.0;
              SearchCorrespondingPoint(sca_fe_ho, lset_vals_ho, orig_point,
                                       lset_lin, normal, final_point, lh);
              for (int d = 0; d < D; ++d) curr_vol_ip(d) = final_point(d);
            
              Vec<D> dist = final_point - orig_point;

              double distnorm = L2Norm(dist);
              if (distnorm > accept_threshold)
              {
                if (distnorm < reject_threshold)
                {
                  dist *= accept_threshold / distnorm;
                  //statistics:
                  *n_corrected_points_edge+=1.0;
                  MappedIntegrationPoint<D,D> mip(old_ip,eltrans);
                  cout << " old_ip = " << old_ip << endl;
                  cout << " corrected mip.GetPoint() = " << mip.GetPoint() << endl;
                }
                else
                {
                  dist = 0.0;
                  //statistics:
                  *n_rejected_points_edge+=1.0;
                }
              }
              else
              {
                //statistics:
                *n_accepted_points_edge+=1.0;
              }


              edge_fe.CalcShape(curr_ip_edge,shape_edge_ho);
              edge_mass_mat += curr_ip_edge.Weight() * shape_only_edge_ho * Trans(shape_only_edge_ho);

              Vec<D> transf_dist = mip_center.GetJacobian() * dist;
              edge_rhs_mat += curr_ip_edge.Weight() * shape_only_edge_ho * Trans(transf_dist);
              
            }
            
            CalcInverse(edge_mass_mat);
            deform_contribution = edge_mass_mat * edge_rhs_mat;
            
            Array<int> dnums_deform;
            fes_deform->GetEdgeDofNrs(global_edge_nr, dnums_deform);

            FlatMatrixFixWidth<D> deform_vec(dnums_deform.Size(),lh);
            FlatVector<> deform_vec_as_vec(D*dnums_deform.Size(),&deform_vec(0,0));
            deform->GetVector().GetIndirect(dnums_deform,deform_vec_as_vec);

            deform_vec += deform_contribution;
            deform->GetVector().SetIndirect(dnums_deform,deform_vec_as_vec);

            // count the number of times that a deformation value
            // has been added for the d.o.f.
            FlatVector<> values(D*dnums_deform.Size(),lh);
            factor->GetIndirect(dnums_deform,values);
            for (int k = 0; k < values.Size(); ++k)
              values(k) += 1.0;
            factor->SetIndirect(dnums_deform,values);
            
          }



          int nfaces = ElementTopology::GetNFaces(eltype);

          if (!no_faces)
          for (int ref_face_nr = 0; ref_face_nr < nfaces; ++ref_face_nr)
          {
            HeapReset hr(lh);

            const int face_order = deform->GetFESpace()->GetOrder();
            if (face_order <= 2) continue;

            ELEMENT_TYPE etface = ElementTopology::GetFaceType (eltype, ref_face_nr);

            if (etface != ET_TRIG)
              throw Exception ("only trig for now - fe creation is hard coded");
            
            const int global_face_nr = faces[ref_face_nr];

            const int v1 = ElementTopology::GetFaces(eltype)[ref_face_nr][0];
            const int v2 = ElementTopology::GetFaces(eltype)[ref_face_nr][1];
            const int v3 = ElementTopology::GetFaces(eltype)[ref_face_nr][2];

            ma->GetFacePNums(global_face_nr, face_verts);
            if (deform->GetFESpace()->IsDirichletFace(global_face_nr))
              continue;

            H1HighOrderFE<ET_TRIG> & face_fe = *(new (lh) H1HighOrderFE<ET_TRIG>(face_order));

            face_fe.SetVertexNumbers(face_verts);

            FlatVector<> shape_face_ho(face_fe.GetNDof(),lh);
            const int face_dofs_offset = 3 * face_order;
            int inner_face_dofs = face_fe.GetNDof() - face_dofs_offset;
            FlatVector<> shape_only_face_ho(inner_face_dofs,&shape_face_ho(face_dofs_offset));
            FlatMatrix<> face_mass_mat(inner_face_dofs,inner_face_dofs,lh);
            FlatMatrixFixWidth<D> face_rhs_mat(inner_face_dofs,lh);
            face_mass_mat = 0.0;
            face_rhs_mat = 0.0;
            
            const IntegrationRule & ir_face = SelectIntegrationRule (etface, 2*face_order);

            IntegrationPoint curr_vol_ip;
            FlatMatrixFixWidth<D> deform_contribution(inner_face_dofs,lh);
            
            for (int l = 0; l < ir_face.GetNIP(); l++)
            {
              const IntegrationPoint & curr_ip_face(ir_face[l]); //one integration point case...
              
              for (int d = 0; d < D; ++d)
              {
                curr_vol_ip(d) = 0.0;
                curr_vol_ip(d) += curr_ip_face(0) * ElementTopology::GetVertices(eltype)[v1][d];
                curr_vol_ip(d) += curr_ip_face(1) * ElementTopology::GetVertices(eltype)[v2][d];
                curr_vol_ip(d) += (1-curr_ip_face(0)-curr_ip_face(1)) * ElementTopology::GetVertices(eltype)[v3][d];
              }

              IntegrationPoint old_ip(curr_vol_ip);
              Vec<D> old_ref_point = old_ip.Point();

              sca_fe_p1.CalcShape(curr_vol_ip,shape_p1);
              double lset_lin = InnerProduct(shape_p1,lset_vals_p1);
            
              Vec<D> orig_point;
              for (int d = 0; d < D; ++d) orig_point(d) = curr_vol_ip(d);
              Vec<D> final_point;
            
              //statistics:
              *n_deformed_points_face+=1.0;
              SearchCorrespondingPoint(sca_fe_ho, lset_vals_ho, orig_point,
                                       lset_lin, normal, final_point, lh);
              for (int d = 0; d < D; ++d) curr_vol_ip(d) = final_point(d);
            
              Vec<D> dist = final_point - orig_point;

              double distnorm = L2Norm(dist);
              if (distnorm > accept_threshold)
              {
                if (distnorm < reject_threshold)
                {
                  dist *= accept_threshold / distnorm;
                  //statistics:
                  *n_corrected_points_face+=1.0;
                }
                else
                {
                  dist = 0.0;
                  //statistics:
                  *n_rejected_points_face+=1.0;
                }
              }
              else
              {
                //statistics:
                *n_accepted_points_face+=1.0;
              }


              face_fe.CalcShape(curr_ip_face,shape_face_ho);
              
              face_mass_mat += curr_ip_face.Weight() * shape_only_face_ho * Trans(shape_only_face_ho);

              Vec<D> transf_dist = mip_center.GetJacobian() * dist;
              // cout << " mip_center.GetJacobian() = " << mip_center.GetJacobian() << endl;
              // cout << " normal = " << normal << endl;
              // cout << " dist = " << dist << endl;
              // cout << " transf_dist = " << transf_dist << endl;
              face_rhs_mat += curr_ip_face.Weight() * shape_only_face_ho * Trans(transf_dist);
              
            }

            CalcInverse(face_mass_mat);
            deform_contribution = face_mass_mat * face_rhs_mat;
             
            Array<int> dnums_deform;
            if (D==2)
              fes_deform->GetInnerDofNrs(elnr, dnums_deform);
            else
              fes_deform->GetFaceDofNrs(global_face_nr, dnums_deform);


            FlatMatrixFixWidth<D> deform_vec(dnums_deform.Size(),lh);
            FlatVector<> deform_vec_as_vec(D*dnums_deform.Size(),&deform_vec(0,0));
            deform->GetVector().GetIndirect(dnums_deform,deform_vec_as_vec);

            deform_vec += deform_contribution;
            deform->GetVector().SetIndirect(dnums_deform,deform_vec_as_vec);

            // count the number of times that a deformation value
            // has been added for the d.o.f.
            FlatVector<> values(D*dnums_deform.Size(),lh);
            factor->GetIndirect(dnums_deform,values);
            for (int k = 0; k < values.Size(); ++k)
              values(k) += 1.0;
            factor->SetIndirect(dnums_deform,values);
            
            // cout << " mip_center.GetPoint() = " << mip_center.GetPoint() << endl;
            // cout << " deform_contribution = " << deform_contribution << endl;
            // cout << " dnums_deform = " << dnums_deform << endl;
            // getchar();
            
          }

          if (!no_cells)
          if (D==3) // cell inner 
          {
            // ... 

          }


          
        }

        // averaging of the (summed) deformation
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

        //statistics:
        cout << " deformpoints(edge) = " << *n_deformed_points_edge << endl;
        cout << " accepted_points(edge) = " << *n_accepted_points_edge << endl;
        cout << " corrected_points(edge) = " << *n_corrected_points_edge << endl;
        cout << " rejected_points(edge) = " << *n_rejected_points_edge << endl;

        cout << " deformpoints(face) = " << *n_deformed_points_face << endl;
        cout << " accepted_points(face) = " << *n_accepted_points_face << endl;
        cout << " corrected_points(face) = " << *n_corrected_points_face << endl;
        cout << " rejected_points(face) = " << *n_rejected_points_face << endl;

        cout << " deformpoints(cell) = " << *n_deformed_points_cell << endl;
        cout << " accepted_points(cell) = " << *n_accepted_points_cell << endl;
        cout << " corrected_points(cell) = " << *n_corrected_points_cell << endl;
        cout << " rejected_points(cell) = " << *n_rejected_points_cell << endl;

        const double n_deformed_points_total
          = *n_deformed_points_edge+*n_deformed_points_face+*n_deformed_points_cell;
        const double n_accepted_points_total
          = *n_accepted_points_edge+*n_accepted_points_face+*n_accepted_points_cell;
        const double n_corrected_points_total
          = *n_corrected_points_edge+*n_corrected_points_face+*n_corrected_points_cell;
        const double n_rejected_points_total
          = *n_rejected_points_edge+*n_rejected_points_face+*n_rejected_points_cell;

        cout << " deformpoints(total) = " << n_deformed_points_total << endl;
        cout << " accepted_points(total) = " << n_accepted_points_total << endl;
        cout << " corrected_points(total) = " << n_corrected_points_total << endl;
        cout << " rejected_points(total) = " << n_rejected_points_total << endl;

        cout << " totalits = " << *n_totalits << endl;
        cout << " totalits/deformpoints = " << *n_totalits/ n_deformed_points_total << endl;
        cout << " maxits = " << *n_maxits << endl;
      }

      ma->SetDeformation(deform);

      // if (D==2)
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
static RegisterNumProc<NumProcGeometryTest<3> > npxgeomtestasdf3d("xgeomtest3d");
