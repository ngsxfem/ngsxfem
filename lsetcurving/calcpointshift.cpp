#include "calcpointshift.hpp"

// using namespace ngsolve;
using namespace ngfem;

namespace ngfem
{ 

  template<int D>
  void CalcGradientOfCoeff(shared_ptr<CoefficientFunction> coef, const MappedIntegrationPoint<D,D>& mip,
                           Vec<D>& der, LocalHeap& lh)
  {
    static Timer time_fct ("CalcGradientOfCoeff");
    RegionTimer reg (time_fct);

    HeapReset hr(lh);
    // bmatu = 0;
    // evaluate dshape by numerical diff
    //fel_u, eltrans, sip, returnval, lh

    const IntegrationPoint& ip = mip.IP();//volume_ir[i];
    const ElementTransformation & eltrans = mip.GetTransformation();
    
    Vec<D> der_ref;
  
    double eps = 1e-7;
    for (int j = 0; j < D; j++)   // d / dxj
    {
      IntegrationPoint ipl(ip);
      ipl(j) -= eps;
      MappedIntegrationPoint<D,D> mipl(ipl, eltrans);

      IntegrationPoint ipr(ip);
      ipr(j) += eps;
      MappedIntegrationPoint<D,D> mipr(ipr, eltrans);

      const double valright = coef->Evaluate(mipr);
      const double valleft = coef->Evaluate(mipl);
      
      der_ref[j] = (1.0/(2*eps)) * (valright-valleft);
    }
    der = Trans(mip.GetJacobianInverse()) * der_ref;
  }

  
  template<int D>
  double LsetEvaluator<D>::Evaluate(const IntegrationPoint & ip, LocalHeap & lh) const
  {
    if (scafe)
    {
      HeapReset hr (lh);
      FlatVector<> shape(scafe->GetNDof(), lh);
      scafe->CalcShape(ip, shape);
      return InnerProduct(shape, scavalues);
    }
    else
    {
      MappedIntegrationPoint<D,D> mip(ip,*eltrans);
      return coef->Evaluate(mip);
    }
  }

  template<int D>
  Vec<D> LsetEvaluator<D>::EvaluateGrad(const IntegrationPoint & ip, LocalHeap & lh) const 
  {
    if (scafe)
    {
      HeapReset hr (lh);
      FlatMatrixFixWidth<D> dshape(scafe->GetNDof(), lh);
      scafe->CalcDShape(ip, dshape);
      return Trans(dshape) * scavalues;
    }
    else
    {
      MappedIntegrationPoint<D,D> mip(ip,*eltrans);
      Vec<D> der;
      CalcGradientOfCoeff(coef, mip, der, lh);
      return Trans(mip.GetJacobian()) * der;
    }
  }


  bool ElementInRelevantBand (shared_ptr<CoefficientFunction> lset_p1,
                              const ElementTransformation & eltrans,
                              double lower_lset_bound, 
                              double upper_lset_bound)
  {
    ELEMENT_TYPE et = eltrans.GetElementType();
    bool has_neg = false;
    bool has_pos = false;
    for (int s = 0; s < ElementTopology::GetNVertices(et); ++s)
    {
      const double * v = ElementTopology::GetVertices(et)[s];
      IntegrationPoint ip(v[0],v[1],v[2]);
      double val;
      if (ElementTopology::GetSpaceDim(et)==2)
      {
        MappedIntegrationPoint<2,2> mip(ip, eltrans);
        val = lset_p1->Evaluate(mip);
      }
      else
      {
        MappedIntegrationPoint<3,3> mip(ip, eltrans);
        val = lset_p1->Evaluate(mip);
      }
      if (val==0.0) has_neg = has_pos = true;
      if (val > lower_lset_bound)
        has_pos = true;
      if (val < upper_lset_bound)
        has_neg = true;
    }
    
    if (!has_neg || !has_pos)
      return false;
    else
      return true;
  }  

  bool ElementInRelevantBand (FlatVector<> lset_p1,
                              double lower_lset_bound, 
                              double upper_lset_bound)
  {
    bool has_neg = false;
    bool has_pos = false;
    for (int s = 0; s < lset_p1.Size(); ++s)
    {
      const double val = lset_p1[s];
      if (val==0.0) has_neg = has_pos = true;
      if (val > lower_lset_bound)
        has_pos = true;
      if (val < upper_lset_bound)
        has_neg = true;
    }
    
    if (!has_neg || !has_pos)
      return false;
    else
      return true;
  }  


  
  template<int D>
  void SearchCorrespondingPoint (
    const LsetEvaluator<D> & lseteval,                               //<- lset_ho
    const Vec<D> & init_point, double goal_val,                      //<- init.point and goal val
    const Mat<D> & trafo_of_normals, const Vec<D> & init_search_dir, //<- search direction
    bool dynamic_search_dir,
    Vec<D> & final_point, LocalHeap & lh,                            //<- result and localheap
    double * n_totalits,
    double * n_maxits)
  {
    static Timer timer("SearchCorrespondingPoint"); 
    RegionTimer reg (timer);
    // static Timer time_not_conv ("SearchCorrespondingPoint::not converged");
    // static Timer time_conv ("SearchCorrespondingPoint::converged");
    // static Timer time_its ("SearchCorrespondingPoint::iterations");
    
    HeapReset hr(lh);
      
    IntegrationPoint curr_ip;
    for (int d = 0; d < D; ++d) curr_ip(d) = init_point(d);

    Vec<D> search_dir = init_search_dir;
    int it = 0;
    for (it = 0; it < 20; ++it)
    {
      //RegionTimer reg_its (time_its);
      const double curr_val = lseteval.Evaluate(curr_ip,lh); // InnerProduct(shape, sca_values);
      const Vec<D> curr_grad = lseteval.EvaluateGrad(curr_ip,lh); //Trans(dshape) * sca_values;
      const double curr_defect = goal_val - curr_val;
      if (abs(curr_defect) < 1e-14) // && it > 2)
        break;

      if (dynamic_search_dir)
      {
        search_dir = trafo_of_normals * curr_grad;
        // search_dir /= L2Norm(search_dir);
      }
        
      const double dphidn = InnerProduct(curr_grad,search_dir);

      for (int d = 0; d < D; ++d)
        curr_ip(d) += curr_defect / dphidn * search_dir(d);
    }

    if (n_totalits)
#pragma omp critical (totalits)
      *n_totalits += it;
#pragma omp critical (maxits)
    if (n_maxits)
      *n_maxits = max2((double)it,*n_maxits);

    if (it == 20){
      //RegionTimer reg (time_not_conv);
      
      std::cout << IM(2) << " SearchCorrespondingPoint:: did not converge " << std::endl;
      // getchar();
      final_point = init_point;
    }
    else
    {
      //RegionTimer reg (time_conv);
      for (int d = 0; d < D; ++d) final_point(d) = curr_ip(d);
    }
  }


  template void CalcGradientOfCoeff<1>
  (shared_ptr<CoefficientFunction>, const MappedIntegrationPoint<1,1>&, Vec<1>&, LocalHeap&);
  template void CalcGradientOfCoeff<2>
  (shared_ptr<CoefficientFunction>, const MappedIntegrationPoint<2,2>&, Vec<2>&, LocalHeap&);
  template void CalcGradientOfCoeff<3>
  (shared_ptr<CoefficientFunction>, const MappedIntegrationPoint<3,3>&, Vec<3>&, LocalHeap&);

  template class LsetEvaluator<1>;
  template class LsetEvaluator<2>;
  template class LsetEvaluator<3>;
  
  template void SearchCorrespondingPoint<1> (const LsetEvaluator<1> &, const Vec<1> &, double, const Mat<1> &, const Vec<1> &, bool, Vec<1> &, LocalHeap &, double *, double *);
  template void SearchCorrespondingPoint<2> (const LsetEvaluator<2> &, const Vec<2> &, double, const Mat<2> &, const Vec<2> &, bool, Vec<2> &, LocalHeap &, double *, double *);
  template void SearchCorrespondingPoint<3> (const LsetEvaluator<3> &, const Vec<3> &, double, const Mat<3> &, const Vec<3> &, bool, Vec<3> &, LocalHeap &, double *, double *);
  
}
