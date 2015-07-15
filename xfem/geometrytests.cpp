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
  

  template<int D>
  class LsetEvaluator
  {
    const ScalarFiniteElement<D> * scafe = NULL;
    FlatVector<> scavalues;
    shared_ptr<CoefficientFunction> coef = NULL;
    const ElementTransformation * eltrans = NULL;
  public:
    LsetEvaluator(const ScalarFiniteElement<D> & sca_fe, FlatVector<> sca_values) :
      scafe(&sca_fe), scavalues(sca_values)
    { ; }

    LsetEvaluator(shared_ptr<CoefficientFunction> acoef, const ElementTransformation & aeltrans) :
      coef(acoef), eltrans(&aeltrans)
    { ; }

    double Evaluate(const IntegrationPoint & ip, LocalHeap & lh) const
    {
      if (scafe)
      {
        HeapReset hr (lh);
        FlatVector<> shape(scafe->GetNDof(), lh);
        scafe->CalcShape(ip, shape);
        return InnerProduct(shape, scavalues);
      }
      else
      {
        MappedIntegrationPoint<D,D> mip(ip,*eltrans);
        return coef->Evaluate(mip);
      }
    }

    Vec<D> EvaluateGrad(const IntegrationPoint & ip, LocalHeap & lh) const 
    {
      if (scafe)
      {
        HeapReset hr (lh);
        FlatMatrixFixWidth<D> dshape(scafe->GetNDof(), lh);
        scafe->CalcDShape(ip, dshape);
        return Trans(dshape) * scavalues;
      }
      else
      {
        MappedIntegrationPoint<D,D> mip(ip,*eltrans);
        Vec<D> der;
        CalcGradientOfCoeff(coef, mip, der, lh);
        return Trans(mip.GetJacobian()) * der;
      }
    }

  };

  template<int D>
  void SearchCorrespondingPoint (
    const LsetEvaluator<D> & lseteval,                               //<- lset_ho
    const Vec<D> & init_point, double goal_val,                      //<- init.point and goal val
    const Mat<D> & trafo_of_normals, const Vec<D> & init_search_dir, //<- search direction
    bool dynamic_search_dir,
    Vec<D> & final_point, LocalHeap & lh,                            //<- result and localheap
    double * n_totalits = nullptr,
    double * n_maxits = nullptr
    )                            
  {
    HeapReset hr(lh);
      
    IntegrationPoint curr_ip;
    for (int d = 0; d < D; ++d) curr_ip(d) = init_point(d);

    Vec<D> search_dir = init_search_dir;

    int it = 0;
    for (it = 0; it < 100; ++it)
    {
      const double curr_val = lseteval.Evaluate(curr_ip,lh); // InnerProduct(shape, sca_values);
      const Vec<D> curr_grad = lseteval.EvaluateGrad(curr_ip,lh); //Trans(dshape) * sca_values;
      const double curr_defect = goal_val - curr_val;
      if (abs(curr_defect) < 1e-14)
        break;

      if (dynamic_search_dir)
      {
        search_dir = trafo_of_normals * curr_grad;
        search_dir /= L2Norm(search_dir);
      }
        
      const double dphidn = InnerProduct(curr_grad,search_dir);

      for (int d = 0; d < D; ++d)
        curr_ip(d) += curr_defect / dphidn * search_dir(d);
    }

    if (n_totalits)
#pragma omp critical (totalits)
      *n_totalits += it;
#pragma omp critical (maxits)
    if (n_maxits)
      *n_maxits = max((double)it,*n_maxits);

    if (it == 100){
      std::cout << " SearchCorrespondingPoint:: did not converge " << std::endl;
      // getchar();
      final_point = init_point;
    }
    else
      for (int d = 0; d < D; ++d) final_point(d) = curr_ip(d);
  }



  template <int D>
  class PsiStarIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_lset_p1;
    shared_ptr<CoefficientFunction> coef_lset_ho;
    double max_deform = -1;
    double lower_lset_bound = 0.0;
    double upper_lset_bound = 0.0;
  public:
    PsiStarIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_lset_p1(coeffs[0]),coef_lset_ho(coeffs[1]) 
    {
      if (coeffs.Size() > 2)
        max_deform = coeffs[2]->EvaluateConst();
      if (coeffs.Size() > 3)
        lower_lset_bound = coeffs[3]->EvaluateConst();
      if (coeffs.Size() > 4)
        upper_lset_bound = coeffs[4]->EvaluateConst();
      // if (D==3)
      //   throw Exception("Implementation only 2D for now");
    }
    virtual ~PsiStarIntegrator(){ ; };

    virtual string Name () const { return "PsiStarIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }

    // Calculates the element matrix
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> elvec,
                       LocalHeap & lh) const
    {
      elvec = 0.0;
      const ScalarFiniteElement<D> & scafe = dynamic_cast<const ScalarFiniteElement<D> &>(fel);
      
      FlatMatrixFixWidth<D> elvecmat(scafe.GetNDof(),&elvec(0));
      elvecmat = 0.0;

      ELEMENT_TYPE et = eltrans.GetElementType();
      bool has_neg = false;
      bool has_pos = false;
      for (int s = 0; s < ElementTopology::GetNVertices(et); ++s)
      {
        const double * v = ElementTopology::GetVertices(et)[s];
        IntegrationPoint ip(v[0],v[1],v[2]);
        MappedIntegrationPoint<D,D> mip(ip, eltrans);
        const double val = coef_lset_p1->Evaluate(mip);
        if (val==0.0) has_neg = has_pos = true;
        if (val > lower_lset_bound)
          has_pos = true;
        if (val < upper_lset_bound)
          has_neg = true;
      }
      if (!has_neg || !has_pos)
        return;
      
      FlatVector<> shape (scafe.GetNDof(),lh);
      
      IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*scafe.Order());
      for (int l = 0 ; l < ir.GetNIP(); l++)
      {
        MappedIntegrationPoint<D,D> mip(ir[l], eltrans);
        scafe.CalcShape(ir[l],shape);

        Vec<D> grad;
        CalcGradientOfCoeff(coef_lset_p1, mip, grad, lh);
        Mat<D> trafo_of_normals = mip.GetJacobianInverse() * Trans(mip.GetJacobianInverse());
        
        Vec<D> normal = mip.GetJacobianInverse() * grad;
        double len = L2Norm(normal);
        normal /= len;

        Vec<D> orig_point;
        for (int d = 0; d < D; ++d)
          orig_point(d) = ir[l](d);

        double goal_val = coef_lset_p1->Evaluate(mip);
        Vec<D> final_point;
        SearchCorrespondingPoint<D>(LsetEvaluator<D>(coef_lset_ho, eltrans),
                                    orig_point, goal_val, 
                                    trafo_of_normals, normal, false,
                                    final_point, lh);
        Vec<D> ref_dist = (final_point - orig_point);
        const double ref_dist_size = L2Norm(ref_dist);
        if ((max_deform >= 0.0) && (ref_dist_size > max_deform))
        {
          ref_dist *= max_deform / ref_dist_size; 
        }

        
        Vec<D> deform = mip.GetJacobian() * ref_dist;


        elvecmat += mip.GetWeight() * shape * Trans(deform);
      }      
      
      
      
    }
  };

  template class PsiStarIntegrator<2>;
  template class PsiStarIntegrator<3>;

  static RegisterLinearFormIntegrator<PsiStarIntegrator<2> > initpsistar2d ("psistar", 2, 2);
  static RegisterLinearFormIntegrator<PsiStarIntegrator<3> > initpsistar3d ("psistar", 3, 2);

  



  template <int D>
  class RestrictedMassIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
    shared_ptr<CoefficientFunction> coef_lset_p1;
    double lower_lset_bound = 0.0;
    double upper_lset_bound = 0.0;
  public:
    RestrictedMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0]), coef_lset_p1(coeffs[1])
    {
      if (coeffs.Size() > 2)
        lower_lset_bound = coeffs[2]->EvaluateConst();
      if (coeffs.Size() > 3)
        upper_lset_bound = coeffs[3]->EvaluateConst();
      // if (D==3)
      //   throw Exception("Implementation only 2D for now");
    }
    virtual ~RestrictedMassIntegrator(){ ; };

    virtual string Name () const { return "RestrictedMassIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }
    virtual bool IsSymmetric () const { return true; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const
    {
      elmat = 0.0;
      const ScalarFiniteElement<D> & scafe = dynamic_cast<const ScalarFiniteElement<D> &>(fel);
      
      elmat = 0.0;

      ELEMENT_TYPE et = eltrans.GetElementType();
      bool has_neg = false;
      bool has_pos = false;
      for (int s = 0; s < ElementTopology::GetNVertices(et); ++s)
      {
        const double * v = ElementTopology::GetVertices(et)[s];
        IntegrationPoint ip(v[0],v[1],v[2]);
        MappedIntegrationPoint<D,D> mip(ip, eltrans);
        const double val = coef_lset_p1->Evaluate(mip);
        if (val==0.0) has_neg = has_pos = true;
        if (val > lower_lset_bound)
          has_pos = true;
        if (val < upper_lset_bound)
          has_neg = true;
      }
      if (!has_neg || !has_pos)
        return;
      
      FlatVector<> shape (scafe.GetNDof(),lh);
      
      IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*scafe.Order());
      for (int l = 0 ; l < ir.GetNIP(); l++)
      {
        MappedIntegrationPoint<D,D> mip(ir[l], eltrans);
        const double coef_val = coef->Evaluate(mip);
        scafe.CalcShape(ir[l],shape);
        elmat += coef_val * mip.GetWeight() * shape * Trans(shape);
      }      
    }
  };

  template class RestrictedMassIntegrator<2>;
  template class RestrictedMassIntegrator<3>;

  static RegisterBilinearFormIntegrator<RestrictedMassIntegrator<2> > initmassstar2d ("massstar", 2, 2);
  static RegisterBilinearFormIntegrator<RestrictedMassIntegrator<3> > initmassstar3d ("massstar", 3, 2);
  


  void PrintConvergenceTable(const Array<double> & tab, string label="")
  {

    ofstream fout("conv_"+label+".out");
    fout << tab;
    cout << endl;
    cout << label << ":" << endl;
    for (int k = 0; k < tab.Size(); ++k)
    {
      cout << setw(16) << tab[k];
      if(k>0)
        cout << "\t" << -log(tab[k]/tab[k-1])/log(2);
      else if (tab.Size()>1)
        cout << "\teoc:";
      cout << endl;
    }
    if (tab.Size()>1)
    {
      cout << setw(16) << "av. eoc:";
      cout << "\t" << -log(tab[tab.Size()-1]/tab[0])/(log(2)*(tab.Size()-1));
    }
    cout << endl;
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
    double surface_ctrl=-1;
    Array<double> volume_errors;
    Array<double> surface_errors;
    Array<double> max_lset_errors;
    Array<double> lset_l1_errors;
    Array<double> diff_phi_l2_errors;
    Array<double> diff_phi_max_errors;
    
    bool dynamic_search_dir=false; //take normal of accurate level set 

    double lower_lset_bound=0.0; //domain of interest for deformation: lower bound
    double upper_lset_bound=0.0; //domain of interest for deformation: lower bound

    bool no_edges = false;
    bool no_faces = false;
    bool no_cells = false;

    bool adaptive = false;
    
    int order = -1;

    
  public:
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


    NumProcGeometryTest (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset", ""), true);
      gf_lset_p1  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p1", ""), false);
      gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_ho", ""), true);
      if (!gf_lset_ho)
        gf_lset_ho  = apde->GetGridFunction (flags.GetStringFlag ("gf_levelset_p2", ""), false);
      deform  = apde->GetGridFunction (flags.GetStringFlag ("deformation", ""));

      order = deform->GetFESpace()->GetOrder();
      
      no_edges = flags.GetDefineFlag("no_edges");
      no_faces = flags.GetDefineFlag("no_faces");
      no_cells = flags.GetDefineFlag("no_cells");

      adaptive = flags.GetDefineFlag("adaptive");
        
      lower_lset_bound = flags.GetNumFlag("lower_lset_bound",0.0);
      upper_lset_bound = flags.GetNumFlag("upper_lset_bound",0.0);
      
      dynamic_search_dir = flags.GetDefineFlag("dynamic_search_dir");
      
      accept_threshold = flags.GetNumFlag("accept_threshold",flags.GetNumFlag("threshold",0.1));
      reject_threshold = flags.GetNumFlag("reject_threshold",accept_threshold * 4.0);
      
      volume_ctrl = flags.GetNumFlag("volume",-1);
      surface_ctrl = flags.GetNumFlag("surface",-1);

      double radius = flags.GetNumFlag("radius",-1);
      if (radius > 0)
      {
        if (D==2)
        {
          volume_ctrl = M_PI * radius * radius;
          surface_ctrl = 2.0 * M_PI * radius;
        }
        else
        {
          volume_ctrl = 4.0/3.0 * M_PI * radius * radius * radius;
          surface_ctrl = 4.0 * M_PI * radius * radius;
        }
      }
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



    

    void Test (LocalHeap & lh)
    {
      cout << setprecision(18) << endl;
      double volume = 0.0;
      double surface = 0.0;
      double max_lset = 0.0;
      double lset_l1 = 0.0;
      double diff_phi_l2 = 0.0;
      double diff_phi_max = 0.0;
      double diff_phi_vol = 0.0;

      Vec<D> diff_phi_max_location;
      Vec<D> diff_phi_max_location2;
      
      int ne=ma->GetNE();

      if (adaptive && (D == 3))
      {
	int nse = ma->GetNSE();
	for (int i = 0; i < nse; i++)
	  Ng_SetSurfaceRefinementFlag (i+1, 0);
      }

      
      for (int elnr = 0; elnr < ne; ++elnr)
      {
        HeapReset hr(lh);
        Ngs_Element ngel = ma->GetElement(elnr);
        ELEMENT_TYPE eltype = ngel.GetType();
        Array<int> dofs;
        gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dofs);
        FlatVector<> vals(dofs.Size(),lh);
        gf_lset_p1->GetVector().GetIndirect(dofs,vals);

        ma->SetDeformation(deform);
        ElementTransformation & eltrans_curved = ma->GetTrafo (ElementId(VOL,elnr), lh);

        ma->SetDeformation(nullptr);
        ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);

        IntegrationPoint ipzero(0.0,0.0,0.0);
        MappedIntegrationPoint<D,D> mx0(ipzero,eltrans);
        
        {

          // element only measure error if "at the interface"
          Array<int> dnums_lset_p1;
          gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dnums_lset_p1);
          FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),lh);
          gf_lset_p1->GetVector().GetIndirect(dnums_lset_p1,lset_vals_p1);
            
          bool has_pos = (lset_vals_p1[D] > lower_lset_bound);
          bool has_neg = (lset_vals_p1[D] < upper_lset_bound);
          for (int d = 0; d < D; ++d)
          {
            if (lset_vals_p1[d] > lower_lset_bound)
              has_pos = true;
            if (lset_vals_p1[d] < upper_lset_bound)
              has_neg = true;
          }
          if (has_pos && has_neg)
          {
            IntegrationRule pir = SelectIntegrationRule (eltrans_curved.GetElementType(), 4*order + 2);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans_curved);
              Vec<D> y = mx0.GetJacobianInverse() * (mip.GetPoint() - mx0.GetPoint()); //point such that level set is approximately that of x_hathat
              IntegrationPoint ipy(y);
              MappedIntegrationPoint<D,D> mipy(ipy,eltrans);

              const double lset_val_transf_p1 = gf_lset_p1->Evaluate(mip);
              const double lset_val = gf_lset_ho->Evaluate(mipy);

              diff_phi_l2 += mip.GetWeight() * sqr(lset_val_transf_p1-lset_val);
              if (abs(lset_val_transf_p1-lset_val) > diff_phi_max)
              {
                diff_phi_max_location = mipy.GetPoint();
                MappedIntegrationPoint<D,D> mip_orig(pir[i], eltrans);
                diff_phi_max_location2 = mip_orig.GetPoint();
                diff_phi_max = abs(lset_val_transf_p1-lset_val);

                // cout << " lset_val_transf_p1 = " << lset_val_transf_p1 << endl;
                // cout << " lset_val = " << lset_val << endl;
              }
              diff_phi_vol += mip.GetWeight();
            }
            if (adaptive)
            {
              Ng_SetRefinementFlag (elnr+1, 1);
            }
          }
          else
            if (adaptive)
            {
              Ng_SetRefinementFlag (elnr+1, 0);
            }
        }

        /*
        {

          IntegrationPoint ip(&fquad_if.points(i,0),0.0); // x_hathat
          MappedIntegrationPoint<D,D> mip(ip, eltrans_curved); // x

          Vec<D> y = mx0.GetJacobianInverse() * (mip.GetPoint() - mx0.GetPoint()); //point such that level set is approximately that of x_hathat
          IntegrationPoint ipy(y);
          MappedIntegrationPoint<D,D> mipy(ipy,eltrans);
          // now mip.GetPoint() == mipy.GetPoint()

          const double lset_val_transf_p1 = gf_lset_p1->Evaluate(mip);
          
          const double lset_val = gf_lset_ho->Evaluate(mipy);


        }
        */



        
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
            
            IntegrationRule pir = SelectIntegrationRule (eltrans_curved.GetElementType(), 2*order + 2);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans_curved);
              volume += mip.GetWeight();
            }

          }
          continue;
        }

        ScalarFieldEvaluator * lset_eval_p
          = ScalarFieldEvaluator::Create(D,*gf_lset_ho,eltrans,lh);

        auto cquad = new CompositeQuadratureRule<D>() ;

        ELEMENT_TYPE et_time = ET_POINT;

        auto xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                       *cquad, lh, 
                                                       2*order, 0, 
                                                       0, 0);

        xgeom->MakeQuadRule();
        FlatXLocalGeometryInformation fxgeom(*xgeom,lh);
        const FlatCompositeQuadratureRule<D> & fcompr(fxgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(NEG));

        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans_curved);
          volume += mip.GetWeight();
        }

        const FlatQuadratureRuleCoDim1<D> & fquad_if(fcompr.GetInterfaceRule());
        for (int i = 0; i < fquad_if.Size(); ++i)
        {
          IntegrationPoint ip(&fquad_if.points(i,0),0.0); // x_hathat
          MappedIntegrationPoint<D,D> mip(ip, eltrans_curved); // x

          Vec<D> y = mx0.GetJacobianInverse() * (mip.GetPoint() - mx0.GetPoint()); //point such that level set is approximately that of x_hathat
          IntegrationPoint ipy(y);
          MappedIntegrationPoint<D,D> mipy(ipy,eltrans);
          // now mip.GetPoint() == mipy.GetPoint()

          // const double lset_val_transf_p1 = gf_lset_p1->Evaluate(mip);
          
          const double lset_val = gf_lset_ho->Evaluate(mipy);

          max_lset = max(abs(lset_val),max_lset);
          
          
          Mat<D,D> Finv = mip.GetJacobianInverse();
          const double absdet = mip.GetMeasure();

          Vec<D> nref = fquad_if.normals.Row(i);
          Vec<D> normal = absdet * Trans(Finv) * nref ;
          double len = L2Norm(normal);
          normal /= len;

          const double weight = fquad_if.weights(i) * len;

          lset_l1 += weight * abs(lset_val);
          surface += weight;
        }        
        
      }
      
      diff_phi_l2 =sqrt(diff_phi_l2/diff_phi_vol);

      if (surface_ctrl>0)
      {
        surface_errors.Append(abs(surface-surface_ctrl));
        PrintConvergenceTable(surface_errors,"surface_error");
      }
      else
        cout << " surface = " << surface << endl;
      if (volume_ctrl>0)
      {
        volume_errors.Append(abs(volume-volume_ctrl));
        PrintConvergenceTable(volume_errors,"volume_error");
      }
      else
        cout << " volume = " << volume << endl;

      lset_l1_errors.Append(lset_l1);
      PrintConvergenceTable(lset_l1_errors,"lset_l1_error");
      
      max_lset_errors.Append(max_lset);
      PrintConvergenceTable(max_lset_errors,"lset_max_error");

      diff_phi_l2_errors.Append(diff_phi_l2);
      PrintConvergenceTable(diff_phi_l2_errors,"phi_l2_error");
      
      diff_phi_max_errors.Append(diff_phi_max);
      PrintConvergenceTable(diff_phi_max_errors,"phi_max_error");

      cout << " diff_phi_max_location = " << diff_phi_max_location << endl;
      cout << " diff_phi_max_location2 = " << diff_phi_max_location2 << endl;
      
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

      double max_dist = 0.0;
      
      shared_ptr<BaseVector> factor = deform->GetVector().CreateVector();
      *factor = 0.0;
      deform->GetVector() = 0.0;
      
      int ne = ma->GetNE();
      // int ned = ma->GetNEdges();

      shared_ptr<FESpace> fes_deform = deform->GetFESpace();

      SetInterpolants(clh);

      if (true)
      {

        // BitArray non_dir_edges(deform->GetFESpace()->GetNDof());
        // non_dir_edges.Clear();

        // for (int ednr = 0; ednr < ne; ++ednr)
        // {
        //   if (deform->GetFESpace()->IsDirichletEdge(ednr))
        //     continue;
        //   deform->GetFESpace()->GetEdgeDofNrs(ednr,dnums);
        //   for( auto dof : dnums )
        //     non_dir_edges.Set(dof);
        // }

        BitArray if_dofs(deform->GetFESpace()->GetNDof());
        if_dofs.Clear();

        Array<int> dnums;
        
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          Array<int> dnums_lset_p1;
          gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dnums_lset_p1);
          FlatVector<> lset_vals_p1(dnums_lset_p1.Size(),clh);
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
          
          deform->GetFESpace()->GetDofNrs(elnr,dnums);
          for( auto dof : dnums)
            if_dofs.Set(dof);
        }

        if_dofs.And(*deform->GetFESpace()->GetFreeDofs());

        // if_dofs.Or(*deform->GetFESpace()->GetFreeDofs());
        

        
        Flags flags;
        flags.SetFlag("symmetric");
        auto blfm = CreateBilinearForm(deform->GetFESpace(),"proj",flags);
        Array<shared_ptr<CoefficientFunction>> coefs_rm;
        coefs_rm.Append(make_shared<ConstantCoefficientFunction>(1.0));
        coefs_rm.Append(gf_lset_p1);
        coefs_rm.Append(make_shared<ConstantCoefficientFunction>(lower_lset_bound));
        coefs_rm.Append(make_shared<ConstantCoefficientFunction>(upper_lset_bound));
        auto bfi_mass = make_shared<RestrictedMassIntegrator<D> > (coefs_rm);
        auto bfi_mass1 = make_shared<BlockBilinearFormIntegrator> (bfi_mass,D,0);
        blfm -> AddIntegrator ( bfi_mass1 );
        auto bfi_mass2 = make_shared<BlockBilinearFormIntegrator> (bfi_mass,D,1);
        blfm -> AddIntegrator ( bfi_mass2);
        if (D==3)
        {
          auto bfi_mass3 = make_shared<BlockBilinearFormIntegrator> (bfi_mass,D,2);
          blfm -> AddIntegrator ( bfi_mass3);
        }

        Flags lfflags;
        auto lfpsi = CreateLinearForm(deform->GetFESpace(),"lfpsistar",lfflags);
        Array<shared_ptr<CoefficientFunction>> coefs;
        coefs.Append(gf_lset_p1);
        coefs.Append(gf_lset_ho);
        coefs.Append(make_shared<ConstantCoefficientFunction>(accept_threshold));
        coefs.Append(make_shared<ConstantCoefficientFunction>(lower_lset_bound));
        coefs.Append(make_shared<ConstantCoefficientFunction>(upper_lset_bound));
        auto lfi_psi = make_shared<PsiStarIntegrator<D>> (coefs);
        lfpsi -> AddIntegrator (lfi_psi);
        // auto lfi_psi = make_shared<SourceIntegrator<D>> (make_shared<ConstantCoefficientFunction>(1.0));        lfpsi -> AddIntegrator (make_shared<BlockLinearFormIntegrator>(lfi_psi,D ,0));

        blfm -> Assemble(clh);
        lfpsi -> Assemble(clh);
        // blfm->GetMatrix().SetInverseType("sparsecholesky");
        deform->GetVector() = *(blfm->GetMatrix().InverseMatrix(&if_dofs)) * lfpsi->GetVector();
      }
      else
      for (int q = 0; q < 2; ++q)
      {
        if (q==0)
        {
          no_faces = true;
          no_edges = false;
        }
        else
        {
          no_faces = false;
          no_edges = true;
        }

        *factor=0.0;
        
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
          // const double lset_center_p1 = InnerProduct(shape_p1,lset_vals_p1);
          // const double h = pow(mip_center.GetJacobiDet(),1.0/D);
          

          sca_fe_p1.CalcDShape(ip_center,dshape_p1);
          Vec<D> grad = Trans(dshape_p1) * lset_vals_p1;

          Mat<D> trafo_of_normals = mip_center.GetJacobianInverse() * Trans(mip_center.GetJacobianInverse());
          Vec<D> normal = trafo_of_normals *  grad;
          double len = L2Norm(normal);
          normal /= len;
          
          Array<int> edges;
          Array<int> faces;
          Array<int> verts;

          ma->GetElVertices(elnr,verts);
          ma->GetElEdges(elnr,edges);
          ma->GetElFaces(elnr,faces);

          const int order = deform->GetFESpace()->GetOrder();
          Array<int> el_dnums_deform;
          deform->GetFESpace()->GetDofNrs(elnr,el_dnums_deform);

          FlatMatrixFixWidth<D> element_coefs(el_dnums_deform.Size(),lh);
          FlatMatrixFixWidth<D> element_count(el_dnums_deform.Size(),lh);
          element_coefs = 0.0;
          element_count = 0.0;

          if (no_faces && !no_edges)
          {
            FlatVector<> el_coefs_as_vec(D*el_dnums_deform.Size(),&element_coefs(0,0));
            deform->GetVector().GetIndirect(el_dnums_deform,el_coefs_as_vec);
          }
          

          // classify degrees of freedoms into verts/edge/face/cell dofs:
          Array<IntRange> dof_range_of_vert;
          Array<IntRange> dof_range_of_edge;
          Array<IntRange> dof_range_of_face;
          IntRange dof_range_of_cell;
          {
            int offset = 0;
            for (int d = 0; d < D+1; ++d)
            {
              dof_range_of_vert.Append(IntRange(offset,offset+1));
              offset++;
            }

            for (int e = 0; e < ElementTopology::GetNEdges(eltype); ++e)
            {
              dof_range_of_edge.Append(IntRange(offset,offset+order-1));
              offset+=order-1;
            }

            for (int f = 0; f < ElementTopology::GetNFaces(eltype); ++f)
            {
              dof_range_of_face.Append(IntRange(offset,offset+(order-2)*(order-1)/2));
              offset+=(order-2)*(order-1)/2;
            }

            if (D==2)
              dof_range_of_cell = IntRange(offset,offset);
            else
            {
              dof_range_of_cell = IntRange(offset,offset+(order-3)*(order-2)*(order-1)/6);
              offset+=(order-3)*(order-2)*(order-1)/6;
            }
            
            if (offset != el_dnums_deform.Size())
            {
              std::cout << " miscounted 1 " << std::endl;
              throw Exception(" I did not count right ");
            }

            // for (int k = 0; k < dof_range_of_vert.Size(); ++k)
            //   cout << " dof_range_of_vert[k] = " << dof_range_of_vert[k] << endl;
            // for (int k = 0; k < dof_range_of_edge.Size(); ++k)
            //   cout << " dof_range_of_edge[k] = " << dof_range_of_edge[k] << endl;
            // for (int k = 0; k < dof_range_of_face.Size(); ++k)
            //   cout << " dof_range_of_face[k] = " << dof_range_of_face[k] << endl;
            // cout << " dof_range_of_cell = " << dof_range_of_cell << endl;

          }
              

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
            
            FlatMatrixFixWidth<D> edge_coefs(inner_edge_dofs,
                                             &element_coefs(dof_range_of_edge[ref_edge_nr].First(),0));

            const IntegrationRule & ir_edge = SelectIntegrationRule (ET_SEGM, 2*edge_order + 2);
            IntegrationPoint curr_vol_ip;

            for (int l = 0; l < ir_edge.GetNIP(); l++)
            {
              const IntegrationPoint & curr_ip_edge(ir_edge[l]); //one integration point case...
              
              for (int d = 0; d < D; ++d)
                curr_vol_ip(d) = curr_ip_edge(0) * ElementTopology::GetVertices(eltype)[v1][d]
                  + (1.0-curr_ip_edge(0)) * ElementTopology::GetVertices(eltype)[v2][d];

              IntegrationPoint old_ip(curr_vol_ip);

              sca_fe_p1.CalcShape(curr_vol_ip,shape_p1);
              double lset_lin = InnerProduct(shape_p1,lset_vals_p1);
            
              Vec<D> orig_point;
              for (int d = 0; d < D; ++d) orig_point(d) = curr_vol_ip(d);
              Vec<D> final_point;
            
              //statistics:
              *n_deformed_points_edge+=1.0;
              SearchCorrespondingPoint<D>(LsetEvaluator<D>(sca_fe_ho, lset_vals_ho),
                                          orig_point,
                                          lset_lin, trafo_of_normals, normal, dynamic_search_dir,
                                          final_point, lh, n_totalits, n_maxits);
              // for (int d = 0; d < D; ++d) curr_vol_ip(d) = final_point(d);
            
              Vec<D> dist = final_point - orig_point;

              double distnorm = L2Norm(dist);
              max_dist = max(max_dist,distnorm);
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
            edge_coefs = edge_mass_mat * edge_rhs_mat;

            for (int l : dof_range_of_edge[ref_edge_nr])
              element_count(l,0) += 1.0;
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

            if (deform->GetFESpace()->IsDirichletFace(global_face_nr))
              continue;
            
            // vertex numbers (w.r.t. to ref. element) of current face:  
            const int * local_vert = ElementTopology::GetFaces(eltype)[ref_face_nr];

            Array<IntRange> dof_range_of_vert_of_face;
            Array<IntRange> dof_range_of_edge_of_face;
            int ndof_all_edges_of_face = 0;

            { // fill local numbering (of verts and edges)
              for (int v = 0; v < ElementTopology::GetNVertices(etface); ++v)
              {
                dof_range_of_vert_of_face.Append(dof_range_of_vert[local_vert[v]]);;
                ndof_all_edges_of_face ++;
              }
            
              for (int e = 0; e < ElementTopology::GetNEdges(etface); ++e)
              {
                int edge_vertex1 = ElementTopology::GetEdges(etface)[e][0];
                int edge_vertex2 = ElementTopology::GetEdges(etface)[e][1];
                int edge_of_cell = ElementTopology::GetEdgeNr(eltype,local_vert[edge_vertex1],
                                                              local_vert[edge_vertex2]);
                dof_range_of_edge_of_face.Append(dof_range_of_edge[edge_of_cell]);
                ndof_all_edges_of_face += dof_range_of_edge[edge_of_cell].Size();
              }
            }

            FlatMatrixFixWidth<D> coefs_edges (ndof_all_edges_of_face,lh);
            int offset = 0;
            { // fill local coefs corresponding to verts and edges
              for (int v = 0; v < ElementTopology::GetNVertices(etface); ++v)
              {
                coefs_edges.Rows(offset,offset+1)
                  = element_coefs.Rows(dof_range_of_vert_of_face[v]);
                offset ++;
              }
              for (int e = 0; e < ElementTopology::GetNEdges(etface); ++e)
              {
                coefs_edges.Rows(offset,offset + dof_range_of_edge_of_face[e].Size())
                  = element_coefs.Rows(dof_range_of_edge_of_face[e]);
                offset += dof_range_of_edge_of_face[e].Size();
              }            
            }

            ma->GetFacePNums(global_face_nr, face_verts);
            H1HighOrderFE<ET_TRIG> & face_fe = *(new (lh) H1HighOrderFE<ET_TRIG>(face_order));
            face_fe.SetVertexNumbers(face_verts);

            FlatVector<> shape_face_ho(face_fe.GetNDof(),lh);
            int inner_face_dofs = face_fe.GetNDof() - offset;

            FlatVector<> shape_only_edges_ho(offset,&shape_face_ho(0));
            FlatVector<> shape_only_face_ho(inner_face_dofs,&shape_face_ho(offset));
            
            FlatMatrix<> face_mass_mat(inner_face_dofs,inner_face_dofs,lh);
            FlatMatrixFixWidth<D> face_rhs_mat(inner_face_dofs,lh);
            face_mass_mat = 0.0;
            face_rhs_mat = 0.0;
            
            const IntegrationRule & ir_face = SelectIntegrationRule (etface, 2 * face_order + 2);

            IntegrationPoint curr_vol_ip;
            FlatMatrixFixWidth<D> face_coefs(inner_face_dofs, &element_coefs(dof_range_of_face[ref_face_nr].First(),0));
            
            for (int l = 0; l < ir_face.GetNIP(); l++)
            {
              const IntegrationPoint & curr_ip_face(ir_face[l]); //one integration point case...
              
              for (int d = 0; d < D; ++d)
              {
                curr_vol_ip(d) = 0.0;
                curr_vol_ip(d) += curr_ip_face(0) * ElementTopology::GetVertices(eltype)[local_vert[0]][d];
                curr_vol_ip(d) += curr_ip_face(1) * ElementTopology::GetVertices(eltype)[local_vert[1]][d];
                curr_vol_ip(d) += (1-curr_ip_face(0)-curr_ip_face(1)) * ElementTopology::GetVertices(eltype)[local_vert[2]][d];
              }

              IntegrationPoint old_ip(curr_vol_ip);
              // Vec<D> old_ref_point = old_ip.Point();

              sca_fe_p1.CalcShape(curr_vol_ip,shape_p1);
              double lset_lin = InnerProduct(shape_p1,lset_vals_p1);
            
              Vec<D> orig_point;
              for (int d = 0; d < D; ++d) orig_point(d) = curr_vol_ip(d);
              Vec<D> final_point;
            
              //statistics:
              *n_deformed_points_face+=1.0;
              SearchCorrespondingPoint<D>(LsetEvaluator<D>(sca_fe_ho, lset_vals_ho),
                                          orig_point,
                                          lset_lin, trafo_of_normals, normal, dynamic_search_dir,
                                          final_point, lh, n_totalits, n_maxits);
              for (int d = 0; d < D; ++d) curr_vol_ip(d) = final_point(d);
            
              Vec<D> dist = final_point - orig_point;

              double distnorm = L2Norm(dist);
              max_dist = max(max_dist,distnorm);
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

              // correction with deformation of edges..
              transf_dist -= Trans(coefs_edges) * shape_only_edges_ho;

              // cout << " mip_center.GetJacobian() = " << mip_center.GetJacobian() << endl;
              // cout << " normal = " << normal << endl;
              // cout << " dist = " << dist << endl;
              // cout << " transf_dist = " << transf_dist << endl;
              face_rhs_mat += curr_ip_face.Weight() * shape_only_face_ho * Trans(transf_dist);
              
            }
            CalcInverse(face_mass_mat);
            element_coefs = 0.0;
            face_coefs = face_mass_mat * face_rhs_mat;

            // cout << " face_coefs = " << face_coefs << endl;
            for (int l : dof_range_of_face[ref_face_nr])
              element_count(l,0) += 1.0;
            
          }
          
          if (!no_cells)
          if (D==3) // cell inner 
          {
            // ... 

          }


          // Write element_coefs into global vector:
          {

            FlatMatrixFixWidth<D> deform_vec(el_dnums_deform.Size(),lh);
            FlatVector<> deform_vec_as_vec(D*el_dnums_deform.Size(),&deform_vec(0,0));
            deform->GetVector().GetIndirect(el_dnums_deform,deform_vec_as_vec);
            deform_vec_as_vec += element_coefs;
            deform->GetVector().SetIndirect(el_dnums_deform,deform_vec_as_vec);

            // count the number of times that a deformation value
            // has been added for the d.o.f.
            FlatVector<> factors(D*el_dnums_deform.Size(),lh);
            factor->GetIndirect(el_dnums_deform,factors);
            for (int k = 0; k < factors.Size(); k+=D)
              factors(k) += element_count(k) ;
            factor->SetIndirect(el_dnums_deform,factors);
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

      cout << " max_dist = " << max_dist << endl;
          
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
