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
    

    auto et = trafo.GetElementType();
    if (! (et == ET_SEGM || et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX) )
      throw Exception("SymbolicCutlfi can only treat simplices right now");

    LevelsetIntegrationDomain lsetintdom_local(lsetintdom);    
    if (lsetintdom_local.GetIntegrationOrder() < 0) // integration order shall not be enforced by lsetintdom
      lsetintdom_local.SetIntegrationOrder(2*fel.Order());
        
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    elvec = 0;

    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    
    if (ir == nullptr)
      return;


    BaseMappedIntegrationRule & mir = trafo(*ir, lh);

    FlatVector<SCAL> elvec1(elvec.Size(), lh);

    FlatMatrix<SCAL> values(ir->Size(), 1, lh);

    /// WHAT FOLLOWS IN THIS FUNCTION IS COPY+PASTE FROM NGSOLVE !!!

    for (auto proxy : proxies)
      {
        // td.Start();
        FlatMatrix<SCAL> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            // cf -> Evaluate (mir, values);
            for (int i=0; i < mir.Size(); i++)
              values(i,0) = cf->Evaluate(mir[i]);

            for (int i = 0; i < mir.Size(); i++)
              proxyvalues(i,k) = mir[i].GetMeasure() * wei_arr[i] * values(i,0);
          }
        // td.Stop();
        // tb.Start();
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
        // tb.Stop();
        elvec += elvec1;
      }

  }
  
}
