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
    virtual double operator()(const Vec<1>& point) const
    {
      throw Exception(" nonono 1");
    }

    virtual double operator()(const Vec<2>& point) const
    {
      throw Exception(" nonono 2");
    }

    virtual double operator()(const Vec<3>& point) const
    {
      throw Exception(" nonono 3");
    }

    virtual double operator()(const Vec<4>& point) const
    {
      throw Exception(" nonono 4");
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
        cout << " D = " << D << endl;
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
      : eval(a_eval), eltrans(a_eltrans) { ; }

    CoefficientFunctionEvaluator( const CoefficientFunction & a_eval, const ElementTransformation & a_eltrans, double a_fixedtime) 
      : eval(a_eval), eltrans(a_eltrans), use_fixedtime(true), fixedtime(a_fixedtime) { ; }


    virtual double operator()(const Vec<1>& point) const
    {
      IntegrationPoint ip(point(0));
      MappedIntegrationPoint<1,D> mip(ip, eltrans);
      if (!fixedtime)
      {
        return eval.Evaluate(mip);
      }
      else
      {
        DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
        mipp.Point().Range(0,D) = mip.GetPoint();
        mipp.Point()[D] = fixedtime;
        return eval.Evaluate(mipp);
      }
    }

    virtual double operator()(const Vec<2>& point) const
    {
      IntegrationPoint ip(point(0), point(1));
      MappedIntegrationPoint<2,D> mip(ip, eltrans);
      if (!fixedtime)
      {
        return eval.Evaluate(mip);
      }
      else
      {
        DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
        mipp.Point().Range(0,D) = mip.GetPoint();
        mipp.Point()[D] = fixedtime;
        return eval.Evaluate(mipp);
      }
    }

    virtual double operator()(const Vec<3>& point) const
    {
      IntegrationPoint ip(point(0),point(1),point(2));
      shared_ptr<BaseMappedIntegrationPoint> mip = NULL;
      if (D==3)
        mip = make_shared< MappedIntegrationPoint<3,3> >(ip, eltrans);
      else
        throw Exception (" Dimensions do not match, D < 3 ? ");
      if (!fixedtime)
      {
        return eval.Evaluate(*mip);
      }
      else
      {
        DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
        mipp.Point().Range(0,D) = mip->GetPoint();
        mipp.Point()[D] = fixedtime;
        return eval.Evaluate(mipp);
      }
    }
  };


  template <int D> // D : resulting space dimension..
  class SpaceTimeCoefficientFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    const CoefficientFunction & eval;
    const ElementTransformation & eltrans;
    TimeInterval ti;
  public:
    SpaceTimeCoefficientFunctionEvaluator( const CoefficientFunction & a_eval, const ElementTransformation & a_eltrans, const TimeInterval & a_ti) 
      : eval(a_eval), eltrans(a_eltrans), ti(a_ti) 
    { 
      // static bool first = true;
      // if (first)
      // {
      //   cout << " WARNING SpaceTimeCoefficientFunctionEvaluator only evaluates as a time-independent lset " << endl;
      //   first = false;
      // }
    }

    virtual double operator()(const Vec<2>& point) const
    {
      IntegrationPoint ip(point(0));
      MappedIntegrationPoint<1,D> mip(ip, eltrans);
      DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
      mipp.Point().Range(0,D) = mip.GetPoint();
      mipp.Point()[D] = (1.0-point(1)) * ti.first + point(1) * ti.second;
      double ret = eval.Evaluate(mipp);
      return ret;
    }

    virtual double operator()(const Vec<3>& point) const
    {
      IntegrationPoint ip(point(0),point(1));
      MappedIntegrationPoint<2,D> mip(ip, eltrans);
      DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
      mipp.Point().Range(0,D) = mip.GetPoint();
      // mipp.Point()[D] = (1.0-point(2)) * ti.first + point(2) * ti.second;
      mipp.Point()[D] = (1.0-point(2)) * ti.first + point(2) * ti.second;
      // std::cout << " mipp.Point() = " << mipp.Point() << std::endl;
      // std::cout << " eval.Evaluate(mipp) = " << eval.Evaluate(mipp) << std::endl;
      double res = eval.Evaluate(mipp);

      return res;
    }

    virtual double operator()(const Vec<4>& point) const
    {
      IntegrationPoint ip(point(0),point(1),point(2));

      shared_ptr<BaseMappedIntegrationPoint> mip = NULL;
      if (D==3)
        mip = make_shared< MappedIntegrationPoint<3,3> >(ip, eltrans);
      else
        throw Exception (" Dimensions do not match, D < 3 ? ");
      DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
      mipp.Point().Range(0,D) = mip->GetPoint();
      mipp.Point()[D] = (1.0-point(3)) * ti.first + point(3) * ti.second;
      return eval.Evaluate(mipp);
    }


  };



} // end of namespace


#endif
