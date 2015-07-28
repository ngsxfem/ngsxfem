#pragma once

#include <fem.hpp>

// using namespace ngsolve;
using namespace ngfem;

namespace ngfem
{ 

  template<int D>
  void CalcGradientOfCoeff(shared_ptr<CoefficientFunction> coef, const MappedIntegrationPoint<D,D>& mip,
                           Vec<D>& der, LocalHeap& lh);
  
  template<int D>
  class LsetEvaluator
  {
    const ScalarFiniteElement<D> * scafe = NULL;
    FlatVector<> scavalues;
    shared_ptr<CoefficientFunction> coef = NULL;
    const ElementTransformation * eltrans = NULL;
  public:
    LsetEvaluator(const ScalarFiniteElement<D> & sca_fe, FlatVector<> sca_values) :
      scafe(&sca_fe), scavalues(sca_values)
    { ; }

    LsetEvaluator(shared_ptr<CoefficientFunction> acoef, const ElementTransformation & aeltrans) :
      coef(acoef), eltrans(&aeltrans)
    { ; }

    double Evaluate(const IntegrationPoint & ip, LocalHeap & lh) const;
    Vec<D> EvaluateGrad(const IntegrationPoint & ip, LocalHeap & lh) const;
  };




  bool ElementInRelevantBand (shared_ptr<CoefficientFunction> lset_p1,
                              const ElementTransformation & eltrans,
                              double lower_lset_bound, 
                              double upper_lset_bound);

  bool ElementInRelevantBand (FlatVector<> lset_p1,
                              double lower_lset_bound, 
                              double upper_lset_bound);
  
  
  template<int D>
  void SearchCorrespondingPoint (
    const LsetEvaluator<D> & lseteval,                               //<- lset_ho
    const Vec<D> & init_point, double goal_val,                      //<- init.point and goal val
    const Mat<D> & trafo_of_normals, const Vec<D> & init_search_dir, //<- search direction
    bool dynamic_search_dir,
    Vec<D> & final_point, LocalHeap & lh,                            //<- result and localheap
    double * n_totalits = nullptr,
    double * n_maxits = nullptr
    );
  
}
