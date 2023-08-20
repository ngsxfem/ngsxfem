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

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    elvec = 0;

    // case: simd enabled; if an error occurs, switch back to normal evaluation
    if (simd_evaluate && globxvar.SIMD_EVAL){
      try
      {
        const IntegrationRule *ns_ir;
        Array<double> ns_wei_arr;
        tie(ns_ir, ns_wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
        if (ns_ir == nullptr)
          return;
        
        SIMD_IntegrationRule simd_ir(*ns_ir, lh);
        
        const int simd_blocks = (ns_ir->Size() + SIMD<IntegrationPoint>::Size() - 1) / SIMD<IntegrationPoint>::Size();
        SIMD<double> *simd_wei_arr = new (lh) SIMD<double>[simd_blocks];
        for (int i = 0; i < simd_blocks; i++)
        {
          simd_wei_arr[i] = [&](int j)
          {
            const int nr = i * SIMD<IntegrationPoint>::Size() + j;
            if (nr < ns_ir->Size())
              return ns_wei_arr[nr];
            else
              return 0.;
          };
        }

        /*
        const SIMD_IntegrationRule *simd_ir;
        Array<SIMD<double>> wei_arr;
        tie(simd_ir, wei_arr) = CreateSIMD_CutIntegrationRule(make_tuple(ns_ir, ns_wei_arr), lh);

        if (simd_ir == nullptr)
          return;
        
        // if (simd_wei_arr == nullptr)
        //   return;

        SIMD_BaseMappedIntegrationRule &mir = trafo(*simd_ir, lh);
        auto &mir = trafo(*simd_ir, lh);
        
        */
        auto &mir = trafo(simd_ir, lh);
        
        for (CoefficientFunction *cf : gridfunction_cfs)
          ud.AssignMemory(cf, simd_ir.GetNIP(), cf->Dimension(), lh);

        PrecomputeCacheCF(cache_cfs, mir, lh);

		    //elvec = 0;

		    for (auto proxy : proxies)
		    {
		      // NgProfiler::StartThreadTimer(telvec_dvec, tid);
              FlatMatrix<SIMD<SCAL>> proxyvalues(proxy->Dimension(), simd_ir.Size(), lh);
              for (size_t k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;

                cf->Evaluate(mir, proxyvalues.Rows(k, k + 1));
                for (size_t i = 0; i < mir.Size(); i++)
                  proxyvalues(k, i) *= mir[i].GetMeasure() * simd_wei_arr[i];
              }

              proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, elvec);
            }

		    //delete[] wei_arr;
		    /// WHAT FOLLOWS IN THIS FUNCTION IS COPY+PASTE FROM NGSOLVE !!!
		    /*
		    for (auto proxy : proxies)
		    {
            // td.Start();
            FlatMatrix<SCAL> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
            {
              ud.testfunction = proxy;
              ud.test_comp = k;
              // cf -> Evaluate (mir, values);
              for (int i = 0; i < mir.Size(); i++)
                values(i, 0) = cf->Evaluate(mir[i]);
              for (int i = 0; i < mir.Size(); i++)
                proxyvalues(i, k) = mir[i].GetMeasure() * ns_wei_arr[i] * values(i, 0);
              }
              // td.Stop();
              // tb.Start();
              proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
              // tb.Stop();
              elvec += elvec1;
            }
		    }*/
      } catch (ExceptionNOSIMD e) {
        cout << IM(6) << e.What() << endl
             << "switching back to standard evaluation" << endl;
        simd_evaluate = false;
        T_CalcElementVector(fel, trafo, elvec, lh);
      }

      return;
	}

    // case: simd disabled
    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    
    if (ir == nullptr)
      return;


    BaseMappedIntegrationRule & mir = trafo(*ir, lh);

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

}