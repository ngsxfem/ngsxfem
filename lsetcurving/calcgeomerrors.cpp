/*********************************************************************/
/* File:   calcgeomerrors.cpp                                        */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. July 2015                                             */
/*********************************************************************/

#include "calcgeomerrors.hpp"
#include "../cutint/xintegration.hpp"

using namespace xintegration;

namespace ngcomp
{ 

  void PrintConvergenceTable(const Array<double> & tab, string label)
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

  
  template<int D>
  void CalcDistances (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, StatisticContainer & cont, LocalHeap & lh){
    auto ma = deform->GetMeshAccess();
    
    int ne=ma->GetNE();

    int order = deform->GetFESpace()->GetOrder();

    double lset_error_l1 = 0.0;
    double lset_error_max = 0.0;
    Vec<D> point_of_max_lset_error;
    IntegrationPoint ip_of_max_lset_error;
    ofstream pointsout("pointsout");
    for (int elnr = 0; elnr < ne; ++elnr)
    {
      HeapReset hr(lh);
      Ngs_Element ngel = ma->GetElement(elnr);
      ELEMENT_TYPE eltype = ngel.GetType();
      Array<int> dofs;
      gf_lset_p1->GetFESpace()->GetDofNrs(elnr,dofs);
      FlatVector<> lset_vals_p1(dofs.Size(),lh);
      gf_lset_p1->GetVector().GetIndirect(dofs,lset_vals_p1);

      ma->SetDeformation(deform);
      ElementTransformation & eltrans_curved = ma->GetTrafo (ElementId(VOL,elnr), lh);

      ma->SetDeformation(nullptr);
      ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);

      IntegrationPoint ipzero(0.0,0.0,0.0);
      MappedIntegrationPoint<D,D> mx0(ipzero,eltrans);



      bool has_pos = (lset_vals_p1[D] > 0.0);
      bool has_neg = (lset_vals_p1[D] < 0.0);
      for (int d = 0; d < D; ++d)
      {
        if (lset_vals_p1[d] > 0.0)
          has_pos = true;
        if (lset_vals_p1[d] < 0.0)
          has_neg = true;
      }

      if (has_neg && has_pos)
      {

        ScalarFieldEvaluator * lset_eval_p
          = ScalarFieldEvaluator::Create(D,*gf_lset_p1,eltrans,lh);

        auto cquad = new CompositeQuadratureRule<D>() ;

        ELEMENT_TYPE et_time = ET_POINT;

        auto xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                       *cquad, lh, 
                                                       2*order, 0, 
                                                       0, 0);

        xgeom->MakeQuadRule();
        FlatXLocalGeometryInformation fxgeom(*xgeom,lh);
        const FlatCompositeQuadratureRule<D> & fcompr(fxgeom.GetCompositeRule<D>());
        // const FlatQuadratureRule<D> & fquad(fcompr.GetRule(NEG));
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
          
          const double lset_val = lset_ho->Evaluate(mipy);

          pointsout << mip.GetPoint()[0] << " " << mip.GetPoint()[1] << " " <<  abs(lset_val) << endl;
          if (abs(lset_val) > lset_error_max)
          {
            lset_error_max = abs(lset_val);
            point_of_max_lset_error = mip.GetPoint();
            ip_of_max_lset_error = ip;
          }
          
          Mat<D,D> Finv = mip.GetJacobianInverse();
          const double absdet = mip.GetMeasure();

          Vec<D> nref = fquad_if.normals.Row(i);
          Vec<D> normal = absdet * Trans(Finv) * nref ;
          double len = L2Norm(normal);
          normal /= len;

          const double weight = fquad_if.weights(i) * len;

          lset_error_l1 += weight * abs(lset_val);
          // surface += weight;
        }        
        pointsout << endl;
      }
    }
    cont.ErrorL1Norm.Append(lset_error_l1);
    cont.ErrorMaxNorm.Append(lset_error_max);
    cout << " point_of_max_lset_error = " << point_of_max_lset_error << endl;
    cout << " ip_of_max_lset_error = " << ip_of_max_lset_error << endl;
  }



  NumProcCalcErrors::NumProcCalcErrors (shared_ptr<PDE> apde, const Flags & flags)
  {
    // lower_lset_bound = flags.GetNumFlag("lower_lset_bound",0.0);
    // upper_lset_bound = flags.GetNumFlag("upper_lset_bound",0.0);
    gf_lset_p1 = apde->GetGridFunction(flags.GetStringFlag("levelset_p1","gf_lset_p1"));
    lset = apde->GetCoefficientFunction(flags.GetStringFlag("levelset","lset"));
    deform = apde->GetGridFunction(flags.GetStringFlag("deform","deform"));
  }
  
  void NumProcCalcErrors::Do (LocalHeap & lh)
  {
    auto ma = deform->GetMeshAccess();
    if (ma->GetDimension() == 2)
      CalcDistances<2>(lset, gf_lset_p1, deform,
                       // lower_lset_bound, upper_lset_bound,
                       error_container, lh);
    else
      CalcDistances<3>(lset, gf_lset_p1, deform,
                       // lower_lset_bound, upper_lset_bound,
                       error_container, lh);



    PrintConvergenceTable(error_container.ErrorL1Norm, "lset_on_gamma_l1");
    PrintConvergenceTable(error_container.ErrorMaxNorm, "lset_on_gamma_max");
  }


}

static RegisterNumProc<NumProcCalcErrors> npcalcerr("calcerrors");

