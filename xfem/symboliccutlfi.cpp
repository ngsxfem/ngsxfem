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

}