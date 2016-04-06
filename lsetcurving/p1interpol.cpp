/*********************************************************************/
/* File:   p1interpol.cpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   July 17th 2015                                            */
/*********************************************************************/

#include "p1interpol.hpp"

namespace ngcomp
{ 

/* ---------------------------------------- 
   Interpolate a coefficient function 
   or an h1ho function into the space of
   piecewise linears
   ---------------------------------------- */
  InterpolateP1::InterpolateP1 (shared_ptr<CoefficientFunction> a_coef, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(a_coef), gf(nullptr), gf_p1(a_gf_p1)
  { ; }
  
  InterpolateP1::InterpolateP1 (shared_ptr<GridFunction> a_gf, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(nullptr), gf(a_gf), gf_p1(a_gf_p1) 
  { ; }

  void InterpolateP1::Do(LocalHeap & lh)
  {
    static Timer time_fct ("LsetCurv::InterpolateP1::Do");
    RegionTimer reg (time_fct);
    
    int nv=ma->GetNV();
    gf_p1->GetVector() = 0.0;
      
    for (int vnr = 0; vnr < nv; ++vnr)
    {
      HeapReset hr(lh);
      // Vec<D> point;
      // ma->GetPoint<D>(vnr,point);
      // Mat<1,D> pointmat;
      // pointmat.Row(0) = point;
      // IntegrationPoint ip(0.0);
      // FE_ElementTransformation<0,D> eltrans(ET_POINT,pointmat);
      // MappedIntegrationPoint<0,D> mip(ip,eltrans);

      double val_lset;
      if (coef)
        throw Exception (" not without dimension hack for MIP<0,D>...");
        // val_lset = coef->Evaluate(mip);
      else
      {
        Array<int> dof;
        gf->GetFESpace()->GetVertexDofNrs(vnr, dof);
        FlatVector<> fval(1,&val_lset);
        gf->GetVector().GetIndirect(dof,fval);
      }
        
      Array<int> dof;
      gf_p1->GetFESpace()->GetVertexDofNrs(vnr,dof);
      FlatVector<> val(1,&val_lset);
      gf_p1->GetVector().SetIndirect(dof,val);
    }


  }

  NumProcInterpolateP1::NumProcInterpolateP1(shared_ptr<PDE> apde, const Flags & flags)
  {
    auto coef = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient", ""), true);
    auto gf_p1 = apde->GetGridFunction (flags.GetStringFlag ("gridfunction_p1", ""), true);
    auto gf_ho = apde->GetGridFunction (flags.GetStringFlag ("gridfunction_ho", ""), true);

    if (!coef && !gf_ho)
      throw Exception("please provide gridfunction_ho or coefficient");
    
    if (!gf_p1)
      throw Exception("please provide gridfunction_p1");

    interpol = make_shared<InterpolateP1>(gf_ho,gf_p1);

  }
  void NumProcInterpolateP1::Do (LocalHeap & lh)
  {
    interpol->Do(lh);
  }
  
}

static RegisterNumProc<NumProcInterpolateP1> npinterpolp1("interpolatep1");
