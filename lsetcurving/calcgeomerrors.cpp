/*********************************************************************/
/* File:   calcgeomerrors.cpp                                        */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. July 2015                                             */
/*********************************************************************/

#include "calcgeomerrors.hpp"
#include "../cutint/xintegration.hpp"
#include "shiftintegrators.hpp"

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
  void CalcDistances (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, StatisticContainer & cont, LocalHeap & lh, double refine_threshold, bool abs_ref_threshold){
    static Timer time_fct ("CalcDistances");
    RegionTimer reg (time_fct);

    auto ma = deform->GetMeshAccess();

    if (refine_threshold > 0)
    {
      int ne = ma->GetNE();
      for (int i = 0; i < ne; i++)
        Ng_SetRefinementFlag (i+1, 0);
      if (D==3)
      {
        int nse = ma->GetNSE();
        for (int i = 0; i < nse; i++)
          Ng_SetSurfaceRefinementFlag (i+1, 0);
      }
    }

    int ne=ma->GetNE();

    int order = deform->GetFESpace()->GetOrder();

    double lset_error_l1 = 0.0;
    double lset_error_max = 0.0;
    Vec<D> point_of_max_lset_error;
    IntegrationPoint ip_of_max_lset_error;
    // ofstream pointsout("pointsout");

    ProgressOutput progress (ma, "calc distance on element", ma->GetNE());

    int marked = 0;

    // IterateElements
    //   (*(gf_lset_p1->GetFESpace()), VOL, clh,  [&] (FESpace::Element el, LocalHeap & lh)
    //    {
    for (auto el : gf_lset_p1->GetFESpace()->Elements(VOL,lh))
    {
      int elnr = el.Nr();
      HeapReset hr(lh);
      progress.Update ();

      Ngs_Element ngel = ma->GetElement(el);
      ELEMENT_TYPE eltype = ngel.GetType();
      Array<int> dofs;
      gf_lset_p1->GetFESpace()->GetDofNrs(el,dofs);
      FlatVector<> lset_vals_p1(dofs.Size(),lh);
      gf_lset_p1->GetVector().GetIndirect(dofs,lset_vals_p1);

      ElementTransformation * eltrans_curved;
      ElementTransformation * eltrans;

// #pragma omp critical(deform)
      {
        ma->SetDeformation(deform);
        eltrans_curved = &ma->GetTrafo (el, lh);
      }
// #pragma omp critical(deform)
      {
        ma->SetDeformation(nullptr);
        eltrans = &ma->GetTrafo (el, lh);
      }
      IntegrationPoint ipzero(0.0,0.0,0.0);
      MappedIntegrationPoint<D,D> mx0(ipzero,*eltrans);


      if (ElementInRelevantBand(lset_vals_p1, 0.0, 0.0))
      {
        const IntegrationRule * ir = CreateCutIntegrationRule(nullptr, gf_lset_p1, *eltrans,
                                                              IF, 2*order, -1, lh, 0);
        const IntegrationRule & fquad_if(*ir);
        
        bool mark_this_el = false;

        for (int i = 0; i < fquad_if.Size(); ++i)
        {
          IntegrationPoint & ip(fquad_if[i]); // x_hathat
          MappedIntegrationPoint<D,D> mip(ip, *eltrans_curved); // x

          Vec<D> y = mx0.GetJacobianInverse() * (mip.GetPoint() - mx0.GetPoint()); //point such that level set is approximately that of x_hathat
          IntegrationPoint ipy(y);
          MappedIntegrationPoint<D,D> mipy(ipy,*eltrans);
          // now mip.GetPoint() == mipy.GetPoint()

          // const double lset_val_transf_p1 = gf_lset_p1->Evaluate(mip);

          const double lset_val = lset_ho->Evaluate(mipy);

          // pointsout << mip.GetPoint()[0] << " " << mip.GetPoint()[1] << " " <<  abs(lset_val) << endl;
          const double h = std::pow(mipy.GetJacobiDet(),1.0/D);
// #pragma omp critical(max)
          if (abs(lset_val) > lset_error_max)
          {
            lset_error_max = abs(lset_val);
            point_of_max_lset_error = mip.GetPoint();
            ip_of_max_lset_error = ip;
          }


          if (refine_threshold > 0)
            if ( (abs_ref_threshold && (abs(lset_val) > refine_threshold))
                 || (!abs_ref_threshold && (abs(lset_val) > refine_threshold * h)) )
            {
              mark_this_el = true;
            }

          const double weight = mip.GetWeight();

// #pragma omp atomic
          lset_error_l1 += weight * abs(lset_val);
          // surface += weight;
        }
        // pointsout << endl;
        // delete cquad;

        if (mark_this_el)
        {
          Ng_SetRefinementFlag (elnr+1, 1);
// #pragma omp atomic
          marked++;
        }
      }
    }
    // });

    progress.Done();

    if (refine_threshold > 0)
      cout << " marked " << marked << " elements for refinement " << endl;

    cont.ErrorL1Norm.Append(lset_error_l1);
    cont.ErrorMaxNorm.Append(lset_error_max);
    // cout << " point_of_max_lset_error = " << point_of_max_lset_error << endl;
    // cout << " ip_of_max_lset_error = " << ip_of_max_lset_error << endl;
  }


  template<int D>
  void CalcDeformationError (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> gf_lset_p1, shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn, StatisticContainer & cont, LocalHeap & clh, double lower_lset_bound, double upper_lset_bound){
    static Timer time_fct ("CalcDeformationError");
    RegionTimer reg (time_fct);

    auto ma = deform->GetMeshAccess();

    int ne=ma->GetNE();

    BitArray el_curved(ne);
    el_curved.Clear();


    int order = deform->GetFESpace()->GetOrder();

    double deform_l2 = 0.0;
    double domain_l2 = 0.0;
    double deform_max = 0.0;
    ProgressOutput progress (ma, "calc deformation on element", ma->GetNE());


    IterateElements
      (*(deform->GetFESpace()), VOL, clh,  [&] (FESpace::Element el, LocalHeap & lh)
    {
      int elnr = el.Nr();
      HeapReset hr(lh);
      progress.Update ();

      Ngs_Element ngel = ma->GetElement(el);
      ELEMENT_TYPE eltype = ngel.GetType();
      Array<int> dofs;
      gf_lset_p1->GetFESpace()->GetDofNrs(el,dofs);
      FlatVector<> lset_vals_p1(dofs.Size(),lh);
      gf_lset_p1->GetVector().GetIndirect(dofs,lset_vals_p1);

      ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);

      IntegrationPoint ipzero(0.0,0.0,0.0);
      MappedIntegrationPoint<D,D> mx0(ipzero,eltrans);

      if (ElementInRelevantBand(lset_vals_p1, lower_lset_bound, upper_lset_bound))
      {
        el_curved.Set(elnr);

        const ScalarFiniteElement<D> & scafe
          = dynamic_cast<const ScalarFiniteElement<D> &>(el.GetFE());
        FlatVector<> shape(scafe.GetNDof(),lh);
        FlatMatrixFixWidth<D> elvec(scafe.GetNDof(),lh);
        FlatVector<> elvec_as_vec(D*scafe.GetNDof(),&elvec(0,0));
        Array<int> dnums;
        deform->GetFESpace()->GetDofNrs(el,dnums);
        deform->GetVector().GetIndirect(dnums,elvec_as_vec);

        IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*scafe.Order());
        for (int l = 0; l < ir.GetNIP(); l++)
        {
          MappedIntegrationPoint<D,D> mip(ir[l], eltrans);
          scafe.CalcShape(ir[l],shape);
          Vec<D> deform_h = Trans(elvec) * shape;

          Vec<D> deform;
          {    // point search:
            Vec<D> grad;
            CalcGradientOfCoeff(gf_lset_p1, mip, grad, lh);
            Mat<D> trafo_of_normals = mip.GetJacobianInverse() * Trans(mip.GetJacobianInverse());

            Vec<D> normal = mip.GetJacobianInverse() * grad;
            double len = L2Norm(normal);
            normal /= len;

            Vec<D> qnormal;
            if (qn)
            {
              qn->Evaluate(mip,qnormal);
              normal = mip.GetJacobianInverse() * qnormal;
              // double len = L2Norm(normal);
              // normal /= len;
              // qnormal /= L2Norm(qnormal);
            }
            Vec<D> orig_point;
            for (int d = 0; d < D; ++d)
              orig_point(d) = ir[l](d);

            double goal_val = gf_lset_p1->Evaluate(mip);
            Vec<D> final_point;
            SearchCorrespondingPoint<D>(LsetEvaluator<D>(lset_ho, eltrans),
                                        orig_point, goal_val,
                                        trafo_of_normals, normal, false,
                                        final_point, lh);
            Vec<D> ref_dist = (final_point - orig_point);
            deform = mip.GetJacobian() * ref_dist;
            // if (qn)
            //   deform = InnerProduct(deform,qnormal) * qnormal;

          }

          deform_h -= deform;

          const double weight = mip.GetWeight();
          const double weight_deform_sqr = weight * sqr(L2Norm(deform_h));
#pragma omp atomic
          deform_l2 += weight_deform_sqr;
#pragma omp atomic
          domain_l2 += weight;
          const double new_max = L2Norm(deform_h);
#pragma omp critical(max)
          if (new_max > deform_max)
            deform_max = new_max;
        }
      }


    });

    progress.Done ();

    deform_l2 /= domain_l2;

    cont.ErrorL2Norm.Append(sqrt(deform_l2));
    cont.ErrorMaxNorm.Append(deform_max);




    double deform_jump_integral = 0.0;
    double facet_integral = 0.0;

    int nf=ma->GetNFacets();

    int dim = ma->GetDimension();
    BitArray fine_facet(nf);
    fine_facet.Clear();
    Array<int> elfacets;
    for (int i = 0; i < ne; ++i)
    {
      ma->GetElFacets(i,elfacets);
      for (int j=0; j<elfacets.Size(); j++)
        fine_facet.Set(elfacets[j]);
    }


    ProgressOutput progress_f (ma, "calc jump of deformation on facet", nf);

    LocalHeap & lh(clh);

    for (int facnr = 0; facnr < nf; ++facnr)
    {
      progress_f.Update();
      if (!fine_facet.Test(facnr)) continue;
      HeapReset hr(lh);

      int el1 = -1;
      int el2 = -1;
      int facnr1 = -1;
      int facnr2 = -1;

      Array<int> elnums, fnums;
      ma->GetFacetElements(facnr,elnums);
      el1 = elnums[0];

      if(elnums.Size()<2) continue;
      el2 = elnums[1];

      if(!el_curved.Test(el1) || !el_curved.Test(el2)) continue;

      ma->GetElFacets(el1,fnums);
      for (int k=0; k<fnums.Size(); k++)
        if(facnr==fnums[k]) facnr1 = k;

      ma->GetElFacets(el2,fnums);
      for (int k=0; k<fnums.Size(); k++)
        if(facnr==fnums[k]) facnr2 = k;

      ElementTransformation & eltrans1 = ma->GetTrafo (ElementId(VOL,el1), lh);
      ElementTransformation & eltrans2 = ma->GetTrafo (ElementId(VOL,el2), lh);
      Array<int> vnums1, vnums2;

      ma->GetElVertices (el1, vnums1);
      ma->GetElVertices (el2, vnums2);

      ELEMENT_TYPE eltype1 = deform->GetFESpace()->GetFE(ElementId(VOL,el1),lh).ElementType();
      ELEMENT_TYPE eltype2 = deform->GetFESpace()->GetFE(ElementId(VOL,el2),lh).ElementType();

      Facet2ElementTrafo transform1(eltype1,vnums1);
      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);

      Vec<D> normal_ref1; //, normal_ref2;
      for (int i=0; i<D; i++) {
        normal_ref1(i) = normals1[facnr1][i];
        // normal_ref2(i) = normals2[LocalFacetNr2][i];
      }

      Facet2ElementTrafo transform2(eltype2,vnums2);
      // const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, facnr1);

      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*deform->GetFESpace()->GetOrder());


      for (int l = 0; l < ir_facet.GetNIP(); l++)
      {
        Vec<D> deform_at_point [2];
        IntegrationPoint ip[2];
        ip[0] = transform1(facnr1, ir_facet[l]);
        ip[1] = transform2(facnr2, ir_facet[l]);

        double weight = 0.0;

        ElementTransformation * eltransj [] = {&eltrans1, &eltrans2};
        for (int j = 0; j < 2; ++j)
        {
          MappedIntegrationPoint<D,D> mip (ip[j], *eltransj[j]);

          { // point search:
            Vec<D> grad;
            CalcGradientOfCoeff(gf_lset_p1, mip, grad, lh);
            Mat<D> trafo_of_normals = mip.GetJacobianInverse() * Trans(mip.GetJacobianInverse());

            Vec<D> normal = mip.GetJacobianInverse() * grad;
            double len = L2Norm(normal);
            normal /= len;

            Vec<D> qnormal;
            if (qn)
            {
              qn->Evaluate(mip,qnormal);
              normal = mip.GetJacobianInverse() * qnormal;
              // double len = L2Norm(normal);
              // normal /= len;
              // qnormal /= L2Norm(qnormal);
            }

            Vec<D> orig_point;
            for (int d = 0; d < D; ++d)
              orig_point(d) = ip[j](d);

            double goal_val = gf_lset_p1->Evaluate(mip);
            Vec<D> final_point;
            SearchCorrespondingPoint<D>(LsetEvaluator<D>(lset_ho, *eltransj[j]),
                                        orig_point, goal_val,
                                        trafo_of_normals, normal, false,
                                        final_point, lh);
            Vec<D> ref_dist = (final_point - orig_point);

            deform_at_point[j] = mip.GetJacobian() * ref_dist;

            // if (qn)
            // {
            //   deform_at_point[j] = InnerProduct(deform_at_point[j],qnormal) * qnormal;
            // }
          }


          Vec<D> normal1 = mip.GetJacobiDet() * Trans (mip.GetJacobianInverse()) * normal_ref1;
          weight = L2Norm (normal1);

        }
        Vec<D> deform_jump = deform_at_point[1] - deform_at_point[0];
        deform_jump_integral += weight * sqr(L2Norm(deform_jump));
        facet_integral += weight;
      }
    }
    progress_f.Update ();
    cont.ErrorMisc.Append(sqrt(deform_jump_integral/facet_integral));
  }

  template void CalcDistances<2>(shared_ptr<CoefficientFunction> , shared_ptr<GridFunction> ,
                                 shared_ptr<GridFunction> , StatisticContainer & , LocalHeap & ,
                                 double , bool );
  template void CalcDistances<3>(shared_ptr<CoefficientFunction> , shared_ptr<GridFunction> ,
                                 shared_ptr<GridFunction> , StatisticContainer & , LocalHeap & ,
                                 double , bool );
  template void CalcDeformationError<2> (shared_ptr<CoefficientFunction> , shared_ptr<GridFunction> ,
                                         shared_ptr<GridFunction> , shared_ptr<CoefficientFunction> ,
                                         StatisticContainer & , LocalHeap & , double , double );
  template void CalcDeformationError<3> (shared_ptr<CoefficientFunction> , shared_ptr<GridFunction> ,
                                         shared_ptr<GridFunction> , shared_ptr<CoefficientFunction> ,
                                         StatisticContainer & , LocalHeap & , double , double );
  
}

