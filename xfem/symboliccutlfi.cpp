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
#include "../cutint/xintegration.hpp"
namespace ngfem
{

  SymbolicCutLinearFormIntegrator ::
  SymbolicCutLinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                   shared_ptr<CoefficientFunction> acf,
                                   DOMAIN_TYPE adt,
                                   int aforce_intorder,
                                   int asubdivlvl,
                                   SWAP_DIMENSIONS_POLICY apol,
                                   VorB vb)
    : SymbolicLinearFormIntegrator(acf,vb,VOL), cf_lset(acf_lset), dt(adt),
      force_intorder(aforce_intorder), subdivlvl(asubdivlvl), pol(apol)
  {
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(cf_lset,subdivlvl);
    lsetintdom = make_shared<LevelsetIntegrationDomain>(cf_lset,gf_lset,adt,-1,-1,subdivlvl,pol);
  }

  SymbolicCutLinearFormIntegrator ::
  SymbolicCutLinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                   shared_ptr<CoefficientFunction> acf,
                                   VorB vb)
    : SymbolicLinearFormIntegrator(acf,vb,VOL),
      cf_lset(lsetintdom_in.GetLevelsetCF()),
      gf_lset(lsetintdom_in.GetLevelsetGF()),
      dt(lsetintdom_in.GetDomainType()),
      force_intorder(-1),
      subdivlvl(lsetintdom_in.GetNSubdivisionLevels()),
      pol(lsetintdom_in.GetSwapDimensionPolicy())
  {
    lsetintdom = make_shared<LevelsetIntegrationDomain>(lsetintdom_in);
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
    static Timer t("symbolicCutLFI - CalcElementVector", 2);
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
    
    RegionTimer reg(t);

    auto et = trafo.GetElementType();
    if (! (et == ET_SEGM || et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX) )
      throw Exception("SymbolicCutBFI can only treat simplices right now");

    if (lsetintdom->GetIntegrationOrder() < 0)
      lsetintdom->SetIntegrationOrder(2*fel.Order());
        
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    elvec = 0;

    const IntegrationRule * ir;
    Array<double> wei_arr;
    // tie (ir1, wei_arr) = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, intorder, time_order, lh, subdivlvl, pol);
    tie (ir, wei_arr) = CreateCutIntegrationRule(*lsetintdom, trafo, lh);
    
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
