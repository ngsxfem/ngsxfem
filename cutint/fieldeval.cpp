/// from ngsxfem
#include "fieldeval.hpp"

namespace ngfem
{

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const EvalFunction & evalf, const ElementTransformation & eltrans,  LocalHeap & a_lh)
  {
    throw Exception(" No evalfunction-evaluator anymore");
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const EvalFunction & evalf, const ElementTransformation & eltrans, double t, LocalHeap & a_lh)
  {
    throw Exception(" No evalfunction-evaluator anymore");
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const EvalFunction & evalf, const ElementTransformation & eltrans, const TimeInterval & ti, LocalHeap & a_lh)
  {
    throw Exception(" No spacetime for now ");
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const CoefficientFunction & evalf, const ElementTransformation & eltrans,  LocalHeap & a_lh)
  {
    switch (dim)
    {
    case 1 :
      // throw Exception(" dimension 1 does not make sense ... ");
      return new (a_lh) CoefficientFunctionEvaluator<1>(evalf, eltrans);
    case 2 :
      return new (a_lh) CoefficientFunctionEvaluator<2>(evalf, eltrans);
    case 3 :
      return new (a_lh) CoefficientFunctionEvaluator<3>(evalf, eltrans);
    default :
      throw Exception(" ScalarFieldEvaluator::Create - Dimension > 3");
      break;
    }
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const CoefficientFunction & evalf, const ElementTransformation & eltrans, const TimeInterval & ti, LocalHeap & a_lh)
  {
    throw Exception(" No spacetime for now ");
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const CoefficientFunction & evalf, const ElementTransformation & eltrans, double t, LocalHeap & a_lh)
  {
    switch (dim)
    {
    case 1 :
      // throw Exception(" dimension 1 does not make sense ... ");
      return new (a_lh) CoefficientFunctionEvaluator<1>(evalf, eltrans, t);
    case 2 :
      return new (a_lh) CoefficientFunctionEvaluator<2>(evalf, eltrans, t);
    case 3 :
      cout << IM(1) << " ScalarFieldEvaluator::Create - eval functions only evaluate in 3 dimensions"
           << " - prescribing the 4th dimension does not make sense" << endl;
      return new (a_lh) CoefficientFunctionEvaluator<3>(evalf, eltrans, t);
    default :
      throw Exception(" ScalarFieldEvaluator::Create - Dimension > 3");
      break;
    }
  }

  ScalarFieldEvaluator* ScalarFieldEvaluator::Create(int dim, const FiniteElement & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh)
  {
    switch (dim)
    {
    case 1 :
      // throw Exception(" dimension 1 does not make sense ... ");
      return new (a_lh) ScalarFEEvaluator<1>(a_fe,a_linvec,a_lh);
    case 2 :
      return new (a_lh) ScalarFEEvaluator<2>(a_fe,a_linvec,a_lh);
    case 3 :
      return new (a_lh) ScalarFEEvaluator<3>(a_fe,a_linvec,a_lh);
    default :
      throw Exception(" ScalarFieldEvaluator::Create - Dimension > 3");
      break;
    }
  }


  template <int D>
  double ScalarFEEvaluator<D> :: operator()(const Vec<D>& point) const
  {
    for (int i = 0; i < D; ++i)
      ip(i) = point(i);

    double ret = 0;
    {
      HeapReset hr(lh);
      FlatVector<> shape(linvec.Size(),lh);
      // if (st_fe)
      // {
      //   if (timefixed)
      //     st_fe->CalcShapeSpaceTime(ip,fixedtime,shape,lh);
      //   else
      //     throw Exception(" you have a spacetime finite element of dim D and evaluate on D-1. Please fix a time level first!");
      // }
      // else
      s_fe->CalcShape(ip,shape);
      ret = InnerProduct(shape,linvec);
    }
    return ret;
  }

  template <int D>
  double ScalarFEEvaluator<D> :: operator()(const Vec<D+1>& point) const
  {
    for (int i = 0; i < D; ++i)
      ip(i) = point(i);
    double ret = 0;
    {
      HeapReset hr(lh);
      FlatVector<> shape(linvec.Size(),lh);
      // if (st_fe)
      //   st_fe->CalcShapeSpaceTime(ip,point(D),shape,lh);
      // else
      throw Exception(" you evaluate in D+1 although you are not a space-time FE!");
      ret = InnerProduct(shape,linvec);
    }
    return ret;
  }

  template class ScalarFEEvaluator<1>;
  template class ScalarFEEvaluator<2>;
  template class ScalarFEEvaluator<3>;



} // end of namespace
