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
  template <int D> 
  InterpolateP1<D>::InterpolateP1 (shared_ptr<CoefficientFunction> a_coef, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(a_coef), gf(nullptr), gf_p1(a_gf_p1)
  { ; }
  
  template <int D> 
  InterpolateP1<D>::InterpolateP1 (shared_ptr<GridFunction> a_gf, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(nullptr), gf(a_gf), gf_p1(a_gf_p1) 
  { ; }

  template <int D> 
  void InterpolateP1<D>::Do(LocalHeap & lh)
  {
    int nv=ma->GetNV();
    gf_p1->GetVector() = 0.0;
      
    for (int vnr = 0; vnr < nv; ++vnr)
    {
      HeapReset hr(lh);
      Vec<D> point;
      ma->GetPoint<D>(vnr,point);
      Mat<1,D> pointmat;
      pointmat.Row(0) = point;
      IntegrationPoint ip(0.0);
      FE_ElementTransformation<0,D> eltrans(ET_POINT,pointmat);
      MappedIntegrationPoint<0,D> mip(ip,eltrans);

      double val_lset;
      if (coef)
        val_lset = coef->Evaluate(mip);
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
    const int D = apde->GetMeshAccess()->GetDimension();
    auto coef = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient", ""), true);
    cout << " coef = " << coef << endl;
    auto gf_p1 = apde->GetGridFunction (flags.GetStringFlag ("gridfunction_p1", ""), true);
    cout << " gf_p1 = " << gf_p1 << endl;
    auto gf_ho = apde->GetGridFunction (flags.GetStringFlag ("gridfunction_ho", ""), true);
    cout << " gf_ho = " << gf_ho << endl;

    if (!coef && !gf_ho)
      throw Exception("please provide gridfunction_ho or coefficient");
    
    if (!gf_p1)
      throw Exception("please provide gridfunction_p1");

    if (D==2)
      if (coef)
        interpol2d = make_shared<InterpolateP1<2>>(coef,gf_p1);
      else
        interpol2d = make_shared<InterpolateP1<2>>(gf_ho,gf_p1);
    else
      if (coef)
        interpol3d = make_shared<InterpolateP1<3>>(coef,gf_p1);
      else
        interpol3d = make_shared<InterpolateP1<3>>(gf_ho,gf_p1);

  }
  void NumProcInterpolateP1::Do (LocalHeap & lh)
  {
    if (interpol2d)
      interpol2d->Do(lh);
    if (interpol3d)
      interpol3d->Do(lh);
  }
  
}

static RegisterNumProc<NumProcInterpolateP1> npinterpolp1("interpolatep1");
