#ifndef FILE_FIELDEVAL_HPP
#define FILE_FIELDEVAL_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
// #include "../spacetime/spacetimefe.hpp"

// #include <solve.hpp>
// using namespace ngsolve;
using namespace ngcomp;

namespace ngfem
{
  typedef std::pair<double,double> TimeInterval;

  class ScalarFieldEvaluator
  {
  public:
    virtual double Evaluate_SD(const FlatVector<>& point) const
    {
      throw Exception(" nonono " + to_string(point.Size()));
    }

    template <int SD>
    double operator()(const Vec<SD>& point) const
    {
      return Evaluate_SD(FlatVector<>(SD,const_cast<double *>(&point(0))));
    }

    static ScalarFieldEvaluator* Create(int dim, const FiniteElement & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh);

    static ScalarFieldEvaluator* Create(int dim, const EvalFunction & evalf, const ElementTransformation& eltrans, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const EvalFunction & evalf, const ElementTransformation& eltrans, const TimeInterval & ti, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const EvalFunction & evalf, const ElementTransformation& eltrans, double t, LocalHeap & a_lh);

    static ScalarFieldEvaluator* Create(int dim, const CoefficientFunction & coeff, const ElementTransformation& eltrans, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const CoefficientFunction & coeff, const ElementTransformation& eltrans, const TimeInterval & ti, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const CoefficientFunction & coeff, const ElementTransformation& eltrans, double t, LocalHeap & a_lh);

  };

  template <int D>
  class ScalarFEEvaluator : public ScalarFieldEvaluator
  {
  protected:
    // const ScalarSpaceTimeFiniteElement<D> * st_fe;
    const ScalarFiniteElement<D> * s_fe;
    FlatVector<> linvec;
    mutable IntegrationPoint ip;
    LocalHeap & lh;
    mutable double fixedtime = 0;
    mutable bool timefixed = false;
  public:
    ScalarFEEvaluator(const FiniteElement & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh)
      : linvec(a_linvec),
      lh(a_lh)
    {
      // st_fe = dynamic_cast< const ScalarSpaceTimeFiniteElement<D> * >(&a_fe);
      s_fe = dynamic_cast< const ScalarFiniteElement<D> * >(&a_fe);

      // if (st_fe == NULL && s_fe == NULL)
      if (s_fe == NULL)
      {
        cout << IM(1) << " D = " << D << endl;
        throw Exception("ScalarFEEvaluator - constructor: cast failed...");
      }
    }

    void FixTime(double a_fixedtime ) const
    {
      timefixed = true;
      fixedtime = a_fixedtime;
    }

    void UnFixTime() const
    {
      timefixed = false;
    }

    virtual double operator()(const Vec<D>& point) const;
    virtual double operator()(const Vec<D+1>& point) const;
  };

  template <int D> // D : resulting space dimension..
  class CoefficientFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    const CoefficientFunction & eval;
    const ElementTransformation & eltrans;
    bool use_fixedtime = false;
    double fixedtime = 0.0;

  public:
    CoefficientFunctionEvaluator( const CoefficientFunction & a_eval, const ElementTransformation & a_eltrans)
      : eval(a_eval), eltrans(a_eltrans) {; }

    CoefficientFunctionEvaluator( const CoefficientFunction & a_eval, const ElementTransformation & a_eltrans, double a_fixedtime)
      : eval(a_eval), eltrans(a_eltrans), use_fixedtime(true), fixedtime(a_fixedtime) {; }


    virtual double Evaluate_SD(const FlatVector<>& point) const
    {
      IntegrationPoint ip(point,1);
      shared_ptr<BaseMappedIntegrationPoint> mip = NULL;
      if (D==point.Size())
        mip = make_shared< MappedIntegrationPoint<D,D> >(ip, eltrans);
      else if (D==point.Size()+1)
        mip = make_shared< MappedIntegrationPoint<D-1,D> >(ip, eltrans);
      else
        throw Exception (" Dimensions do not match");
      if (!fixedtime)
      {
        return eval.Evaluate(*mip);
      }
      else
      {
        throw Exception (" Is this still used somewhere ? ");
      }
    }
  };

} // end of namespace


#endif
