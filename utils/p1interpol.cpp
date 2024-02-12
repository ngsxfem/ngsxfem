/*********************************************************************/
/* File:   p1interpol.cpp                                            */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   July 17th 2015                                            */
/*********************************************************************/

#include "p1interpol.hpp"
#include "../cutint/spacetimecutrule.hpp"
#include "../utils/ngsxstd.hpp"

namespace ngcomp
{

/* ----------------------------------------
   Interpolate a coefficient function
   or an h1ho function into the space of
   piecewise linears
   ---------------------------------------- */
  InterpolateP1::InterpolateP1 (shared_ptr<CoefficientFunction> a_coef, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(a_coef), gf(nullptr), gf_p1(a_gf_p1)
  {; }

  InterpolateP1::InterpolateP1 (shared_ptr<GridFunction> a_gf, shared_ptr<GridFunction> a_gf_p1)
    : ma(a_gf_p1->GetMeshAccess()), coef(nullptr), gf(a_gf), gf_p1(a_gf_p1)
  {; }

  void InterpolateP1::Do(LocalHeap & lh, double eps_perturbation, double tref_val)
  {
    static int timer = NgProfiler::CreateTimer ("LsetCurv::InterpolateP1::Do"); NgProfiler::RegionTimer reg (timer);

    int nv=ma->GetNV();
    gf_p1->GetVector() = 0.0;

    for (int vnr = 0; vnr < nv; ++vnr)
    {
      HeapReset hr(lh);

      double val_lset;
      if (coef)
      {
        /*
        Array<int> elnums;
        ma -> GetVertexElements (vnr, elnums);
        Ngs_Element ngel = ma -> GetElement (ElementId(VOL,elnums[0]));
        */
        Ngs_Element ngel = ma -> GetElement (ElementId(VOL,ma->GetVertexElements (vnr)[0]));
        auto & eltrans = ma -> GetTrafo (ngel, lh);

        if( ma -> GetDimension() == 2)
        {
          Vec<2> point;
          ma->GetPoint<2>(vnr,point);
          IntegrationPoint ip(0,0,0,0);
          MappedIntegrationPoint<2,2> mip(ip,eltrans);
          Vec<2> refpoint = mip.GetJacobianInverse() * (point - mip.GetPoint());
          IntegrationPoint ip_f(refpoint[0],refpoint[1],0.0,tref_val);
          if (tref_val >= 0)
            MarkAsSpaceTimeIntegrationPoint(ip_f);
          MappedIntegrationPoint<2,2> mip_f(ip_f,eltrans);
          val_lset = coef->Evaluate(mip_f);
        }
        else if( ma -> GetDimension() == 1)
        {
          Vec<1> point;
          ma->GetPoint<1>(vnr,point);
          IntegrationPoint ip(0,0,0,0);
          MappedIntegrationPoint<1,1> mip(ip,eltrans);
          Vec<1> refpoint = mip.GetJacobianInverse() * (point - mip.GetPoint());
          IntegrationPoint ip_f(refpoint[0],0.0,0.0,tref_val);
          if (tref_val >= 0)
            MarkAsSpaceTimeIntegrationPoint(ip_f);
          MappedIntegrationPoint<1,1> mip_f(ip_f,eltrans);
          val_lset = coef->Evaluate(mip_f);
        }
        else if ( ma -> GetDimension() == 3)
        {
          Vec<3> point;
          ma->GetPoint<3>(vnr,point);
          IntegrationPoint ip(0,0,0,0);
          
          MappedIntegrationPoint<3,3> mip(ip,eltrans);
          Vec<3> refpoint = mip.GetJacobianInverse() * (point - mip.GetPoint());
          IntegrationPoint ip_f(refpoint,tref_val);
          if (tref_val >= 0)
            MarkAsSpaceTimeIntegrationPoint(ip_f);
          MappedIntegrationPoint<3,3> mip_f(ip_f,eltrans);
          val_lset = coef->Evaluate(mip_f);
        }
        else
          throw Exception ("D==0 not yet implemnted");
      }
      else
      {
        Array<int> dof;
        // gf->GetFESpace()->GetVertexDofNrs(vnr, dof);
        gf->GetFESpace()->GetDofNrs(NodeId(NT_VERTEX,vnr), dof); //vnr, dof);
        FlatVector<> fval(1,&val_lset);
        gf->GetVector().GetIndirect(dof,fval);
      }

      Array<int> dof;
      gf_p1->GetFESpace()->GetVertexDofNrs(vnr,dof);
      FlatVector<> val(1,&val_lset);
      // avoid vertex cuts by introducing a small perturbation:
      if (abs(val(0)) < eps_perturbation)
        val(0) = eps_perturbation;
      if (dof[0] != -1)
        gf_p1->GetVector().SetIndirect(dof,val);
    }


  }

}
