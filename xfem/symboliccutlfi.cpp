/*********************************************************************/
/* File:   symboliccutlfi.cpp                                        */
/* Author: Christoph Lehrenfeld based on symbolicintegrator.cpp      */
/*         from Joachim Schoeberl (in NGSolve)                       */
/* Date:   September 2016                                            */
/*********************************************************************/
/* 
   Symbolic cut integrators
*/

#include <fem.hpp>
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../cutint/xintegration.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"
#include "../cutint/spacetimecutrule.hpp"
namespace ngfem
{


  SymbolicCutLinearFormIntegrator ::
  SymbolicCutLinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                   shared_ptr<CoefficientFunction> acf,
                                   VorB vb)
    : SymbolicLinearFormIntegrator(acf,vb,VOL),  
      lsetintdom(lsetintdom_in)
  {
    ;

  }

  
  void 
  SymbolicCutLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  void 
  SymbolicCutLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  template <typename SCAL>
  void
  SymbolicCutLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    static Timer timer("symbolicCutLFI - CalcElementVector");
    RegionTimer reg (timer);
    HeapReset hr(lh);
    
    // tstart.Start();
    
    if (element_vb != VOL)
    {
      switch (trafo.SpaceDim())
      {
        // case 1:
        //   T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
        //   return;
        // case 2:
        //   T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
        //   return;
        // case 3:
        //   T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
        //   return;
        default:
          throw Exception ("symbolicCutLFI, EB not yet implemented");
          // throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
        }
    }
  
    // only 1D, 2D and some 3D integration rules are supported yet
    auto et = trafo.GetElementType();
    if (!(et == ET_SEGM || et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX))
      throw Exception("SymbolicCutlfi can only treat simplices right now");

    LevelsetIntegrationDomain lsetintdom_local(lsetintdom);
    if (lsetintdom_local.GetIntegrationOrder() < 0) // integration order shall not be enforced by lsetintdom
      lsetintdom_local.SetIntegrationOrder(2 * fel.Order());

    ProxyUserData ud(proxies.Size(), gridfunction_cfs.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel; // is this important (or needed at all)?

    elvec = 0;
    const IntegrationRule *ir;
    Array<double> wei_arr;
    tie(ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    if (ir == nullptr)
      return;

    // case: simd enabled; if an error occurs, switch back to normal evaluation
    if (simd_evaluate && globxvar.SIMD_EVAL){
      try
      {
        
        SIMD_IntegrationRule & simd_ir = *(new (lh) SIMD_IntegrationRule(*ir, lh));
        FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(wei_arr, lh);        
        auto &simd_mir = trafo(simd_ir, lh);

        PrecomputeCacheCF(cache_cfs, simd_mir, lh);
        for (CoefficientFunction *cf : gridfunction_cfs)
          ud.AssignMemory(cf, simd_ir.GetNIP(), cf->Dimension(), lh);

		    for (auto proxy : proxies)
		    {
		      // NgProfiler::StartThreadTimer(telvec_dvec, tid);
              FlatMatrix<SIMD<SCAL>> proxyvalues(proxy->Dimension(), simd_ir.Size(), lh);
              for (size_t k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;

                cf->Evaluate(simd_mir, proxyvalues.Rows(k, k + 1));
                for (size_t i = 0; i < simd_mir.Size(); i++)
                  proxyvalues(k, i) *= simd_mir[i].GetMeasure() * simd_wei_arr[i];
              }

              proxy->Evaluator()->AddTrans(fel, simd_mir, proxyvalues, elvec);
            }
      } catch (ExceptionNOSIMD e) {
        cout << IM(6) << e.What() << endl
             << "switching back to standard evaluation" << endl;
        simd_evaluate = false;
        T_CalcElementVector(fel, trafo, elvec, lh);
      }
      return;
	}

  BaseMappedIntegrationRule & mir = trafo(*ir, lh);

  PrecomputeCacheCF(cache_cfs, mir, lh);
  for (CoefficientFunction *cf : gridfunction_cfs)
    ud.AssignMemory(cf, ir->GetNIP(), cf->Dimension(), lh);

	FlatVector<SCAL> elvec1(elvec.Size(), lh);
	elvec1 = 0.0;
	FlatMatrix<SCAL> values(ir->Size(), 1, lh);
  
	/// WHAT FOLLOWS IN THIS FUNCTION IS COPY+PASTE FROM NGSOLVE !!!

	for (auto proxy : proxies){
	  // td.Start();
      FlatMatrix<SCAL> proxyvalues(mir.Size(), proxy->Dimension(), lh);
      for (int k = 0; k < proxy->Dimension(); k++){
        ud.testfunction = proxy;
        ud.test_comp = k;
        cf -> Evaluate (mir, values);
        // for (int i=0; i < mir.Size(); i++)
        //   values(i,0) = cf->Evaluate(mir[i]);

        for (int i = 0; i < mir.Size(); i++)
          proxyvalues(i,k) = mir[i].GetMeasure() * wei_arr[i] * values(i,0);
      }
      // td.Stop();
      // tb.Start();
      proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
      // tb.Stop();
      elvec += elvec1;
    }
    
    return;
  }


  SymbolicCutFacetLinearFormIntegrator ::
  SymbolicCutFacetLinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                   shared_ptr<CoefficientFunction> acf,
                                   VorB vb)
    : SymbolicFacetLinearFormIntegrator(acf,vb) 
  {
   
    lsetintdom = make_shared<LevelsetIntegrationDomain>(lsetintdom_in);
    time_order = lsetintdom_in.GetTimeIntegrationOrder();
  }

  template <typename TSCAL>
  void
  SymbolicCutFacetLinearFormIntegrator ::
  T_CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr,
                            const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                            const ElementTransformation & strafo,
                            FlatVector<TSCAL> elvec,
                            LocalHeap & lh) const
  {
    //throw Exception("SymbolicCutFacetLinearFormIntegrator::T_CalcFacetVector not yet implemented");   
  
    static Timer t("SymbolicCutFacetLFI::CalcFacetVector - boundary", NoTracing);
    HeapReset hr(lh);

    elvec = 0;
    
    FlatVector<TSCAL> elvec1(elvec.Size(), lh);

    int maxorder = fel1.Order();

    auto etvol = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (etvol, LocalFacetNr);

    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(etvol, ElVertices);  

    const IntegrationRule * ir_scr = nullptr;
   
    Array<double> wei_arr;
    auto gflset = lsetintdom->GetLevelsetGF();
    if(gflset == nullptr) throw Exception("No gf in SymbolicCutFacetBilinearFormIntegrator::T_CalcFacetVector :(");

    Array<DofId> dnums(0,lh);
    gflset->GetFESpace()->GetDofNrs(trafo1.GetElementId(),dnums);
    FlatVector<> elvec_lset(dnums.Size(),lh);
    gflset->GetVector().GetIndirect(dnums,elvec_lset);

    shared_ptr<SpaceTimeFESpace> st_FE = nullptr;
    shared_ptr<NodalTimeFE> time_FE = nullptr;
    if(time_order > -1){
      st_FE = dynamic_pointer_cast<SpaceTimeFESpace >(gflset->GetFESpace());
      if(st_FE == nullptr) 
        throw Exception("Unable to cast SpaceTimeFESpace in SymbolicCutFacetLinearFormIntegrator::T_CalcFacetVector");
      time_FE = dynamic_pointer_cast< NodalTimeFE>(st_FE->GetTimeFE());
      if(time_FE == nullptr) 
        throw Exception("Unable to cast time finite element in SymbolicCutFacetLinearFormIntegrator::T_CalcFacetVector");
    }

    int NV = ElementTopology::GetNVertices(etfacet); // number of vertices of et_facet

    if (Dim(etfacet) == 2) 
    {
      
      static Timer t("symbolicCutBFI - CoDim2", NoTracing);
      RegionTimer reg (t);

      IVec<4> int_tuple = 
        SwitchET<ET_HEX,ET_TET,ET_PRISM,ET_PYRAMID> (etvol, [&LocalFacetNr, &ElVertices] (auto et)
         { return ET_trait<et>::GetFaceSort(LocalFacetNr,ElVertices); });

      FlatVector<> lset_fv(NV,lh);
      if(time_order < 0) {
        for (int i = 0; i < NV; i++)
          lset_fv(i) = elvec_lset[int_tuple[i]];
        ir_scr = StraightCutIntegrationRuleUntransformed(lset_fv, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
      }
      else {
        int M = time_FE->GetNDof(); //time nodes
        int L = dnums.Size()/M; // dofs (or verts) per vol element 
        FlatVector<> cf_lset_at_element(NV*M, lh);
        FlatMatrix<> cf_lset_at_element_as_mat(M,NV,&cf_lset_at_element(0));

        for (int i = 0; i < NV; i++)
          cf_lset_at_element_as_mat.Col(i) = elvec_lset.Slice(int_tuple[i],L);

        tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(cf_lset_at_element, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
      }
      if(ir_scr == nullptr) return;
      //cout << "ir_scr: " << *ir_scr << endl;
      if (lsetintdom->GetDomainType() == IF) // correct measure for 3D -> 1D case ( codim1 - adjustment is done automatically)
      {
        if (time_order > -1)
          throw Exception("correcting the IF-measure in space-time case is still missing");
        vector<double> lset(NV);
        for(int i=0; i<NV; i++) 
          lset[i] = lset_fv(i);
        LevelsetWrapper lsw(lset, etfacet);
        Vec<3> tetdiffvec2(0.);

        IntegrationRule & ir_scr_intet2 = transform1( LocalFacetNr, (*ir_scr), lh);
        MappedIntegrationRule<3,3> mir3(ir_scr_intet2,trafo1,lh);
        int npoints = ir_scr->Size();

        for (int i = 0; i < npoints; i++)
        {
            IntegrationPoint & ip = (*ir_scr)[i];
            Vec<3> normal = lsw.GetNormal(ip.Point());
            Vec<2> tang = {normal[1],-normal[0]};

            tetdiffvec2 = transform1.GetJacobian( LocalFacetNr, lh) * tang;
            auto F = mir3[i].GetJacobian();
            auto mapped_tang = F * tetdiffvec2;
            const double ratio_meas1D = L2Norm(mapped_tang);
            ip.SetWeight((*ir_scr)[i].Weight() * ratio_meas1D);
        }
      }
    }
    else if (Dim(etfacet) == 1) 
    {

        IVec<2> int_tuple = 
          SwitchET<ET_TRIG,ET_QUAD> (etvol, [&LocalFacetNr, &ElVertices] (auto et) { return ET_trait<et>::GetEdgeSort(LocalFacetNr,ElVertices); });

        if(time_order < 0) {
            Vec<2> lset_vals_edge = {elvec_lset[int_tuple[0]],elvec_lset[int_tuple[1]]}; 
            ir_scr = StraightCutIntegrationRuleUntransformed(lset_vals_edge, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
        }
        else {
            int M = time_FE->GetNDof(); //time nodes
            int L = dnums.Size()/M; // dofs (or verts) per vol element 
            FlatVector<> cf_lset_at_element(NV*M, lh);
            FlatMatrix<> cf_lset_at_element_as_mat(M,NV,&cf_lset_at_element(0));

            for (int i = 0; i < NV; i++)
              cf_lset_at_element_as_mat.Col(i) = elvec_lset.Slice(int_tuple[i],L);

            tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(cf_lset_at_element, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
        }
        if (ir_scr == nullptr) return;
    }
    else if (Dim(etfacet) == 0) {
      FlatVector <> lset_on_facet(elvec_lset.Size()/2,lh);
      lset_on_facet = elvec_lset.Slice(LocalFacetNr,2);
      //LocalFacetNr1
        if(time_order < 0) {
            ir_scr = StraightCutIntegrationRuleUntransformed(lset_on_facet, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
        }
        else {
          tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(lset_on_facet, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
        }
        if (ir_scr == nullptr) return;

        //throw Exception("no 1D cut facet integration provided yet");
    }
    else
      throw Exception("no 1D cut facet integration provided yet");

    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, (*ir_scr), lh);
    MarkAsSpaceTimeIntegrationRule(ir_facet_vol1);

    if(transform1.FacetType(LocalFacetNr) == ET_POINT){
        Vec<3> right_pnt;
        right_pnt = ir_facet_vol1[0].Point();
        for(int i=1; i<ir_facet_vol1.Size(); i++) ir_facet_vol1[i].Point() = right_pnt;
    }

    // evaluate proxy-values
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    //const_cast<ElementTransformation&>(strafo).userdata = &ud;

    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);

    //PrecomputeCacheCF(cache_cfs, mir1, lh);

    RegionTimer reg(t);
     
    FlatMatrix<TSCAL> val(mir1.Size(), 1,lh);
    for (auto proxy : proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<TSCAL> proxyvalues(mir1.Size(), proxy->Dimension(), lh);
     
        mir1.ComputeNormalsAndMeasure (etvol, LocalFacetNr);
        
	for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);  // needed for grad(u), mesh_size, but: index = material index
            // cf -> Evaluate (smir, val);
            proxyvalues.Col(k) = val.Col(0);
	    //cout << "proxyvalues << " << proxyvalues << endl; 
          }

          if (lsetintdom->GetDomainType() == IF) {
            for (int i = 0; i < mir1.Size(); i++) {
             proxyvalues.Row(i) *= ( time_order < 0 ? (*ir_scr)[i].Weight() : wei_arr[i]); 
            }    
	  }
	  else // codim 1
          {
             // throw Exception("Foo!");
	    for (int i = 0; i < mir1.Size(); i++){
	      proxyvalues.Row(i) *= mir1[i].GetMeasure() * ( time_order < 0 ? (*ir_scr)[i].Weight() : wei_arr[i]);
             }
          }

        elvec1 = 0.0;
        proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, elvec1, lh);
	elvec += elvec1;
      }
  
  
  }


}
