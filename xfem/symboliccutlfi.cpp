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
                                   int asubdivlvl)
    : SymbolicLinearFormIntegrator(acf,VOL,false), cf_lset(acf_lset), dt(adt),
    force_intorder(aforce_intorder), subdivlvl(asubdivlvl)
  {
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(cf_lset,subdivlvl);
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
    static Timer t("symbolicCutLFI - CalcElementMatrix", 2);
    HeapReset hr(lh);
    
    // tstart.Start();
    
    if (element_boundary)
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

    int intorder = 2*fel.Order();

    auto et = trafo.GetElementType();
    if (! (et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX) )
      throw Exception("SymbolicCutBFI can only treat simplices right now");

    if (force_intorder >= 0)
      intorder = force_intorder;
    
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    elvec = 0;

    const IntegrationRule * ir = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, intorder, lh, subdivlvl);
    if (ir == nullptr)
      return;
    ///
    const IntegrationRule * ir2 = nullptr;
    if (time_order > -1) //simple tensor product rule (no moving cuts with this..) ...
    {
       static bool warned = false;
       if (!warned)
       {
         cout << "WARNING: This is a pretty simple tensor product rule in space-time.\n";
         cout << "         A mapped integration rule of this will not see the time,\n";
         cout << "         but the underlying integration rule will." << endl;
         warned = true;
       }
       auto ir1D = SelectIntegrationRule (ET_SEGM, 3);
       ir2 = new (lh) IntegrationRule(ir->Size()*ir1D.Size(),lh);
       for (int i = 0; i < ir1D.Size(); i ++)
         for (int j = 0; j < ir->Size(); j ++)
           (*ir2)[i*ir->Size()+j] = IntegrationPoint((*ir)[j](0),(*ir)[j](1),ir1D[i](0),(*ir)[j].Weight()*ir1D[i].Weight());
       // cout << *ir2<< endl;
    }
    else
      ir2 = ir;

    BaseMappedIntegrationRule & mir = trafo(*ir2, lh);

    FlatVector<SCAL> elvec1(elvec.Size(), lh);
    
    FlatMatrix<SCAL> values(ir->Size(), 1, lh);

    /// WHAT FOLLOWS IN THIS FUNCTION IS COPY+PASTE FROM NGSOLVE !!!

    for (auto proxy : proxies)
      {
        // td.Start();
        FlatMatrix<SCAL> proxyvalues(ir->Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            
            cf -> Evaluate (mir, values);
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
          }
        // td.Stop();
        // tb.Start();
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
        // tb.Stop();
        elvec += elvec1;
      }
      
  }
  
}
