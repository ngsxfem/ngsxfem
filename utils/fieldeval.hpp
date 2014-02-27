#ifndef FILE_FIELDEVAL_HPP
#define FILE_FIELDEVAL_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
#include "../spacetime/spacetimefe.hpp"

using namespace ngsolve;

namespace ngfem
{
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
    const ScalarSpaceTimeFiniteElement<D> * st_fe;
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
      st_fe = dynamic_cast< const ScalarSpaceTimeFiniteElement<D> * >(&a_fe);
      s_fe = dynamic_cast< const ScalarFiniteElement<D> * >(&a_fe);

      if (st_fe == NULL && s_fe == NULL)
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
  class EvalFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    const EvalFunction & eval;
    const ElementTransformation & eltrans;
    bool use_fixedtime = false;
    double fixedtime = 0.0;
  public:
    EvalFunctionEvaluator( const string & str, const ElementTransformation & a_eltrans) 
      : eval(str), eltrans(a_eltrans) { ; }

    EvalFunctionEvaluator( const string & str, const ElementTransformation & a_eltrans, double a_fixedtime) 
      : eval(str), eltrans(a_eltrans), use_fixedtime(true), fixedtime(a_fixedtime) { ; }

    EvalFunctionEvaluator( const EvalFunction & str, const ElementTransformation & a_eltrans) 
      : eval(str), eltrans(a_eltrans) { ; }

    EvalFunctionEvaluator( const EvalFunction & str, const ElementTransformation & a_eltrans, double a_fixedtime) 
      : eval(str), eltrans(a_eltrans), use_fixedtime(true), fixedtime(a_fixedtime) { ; }

    virtual double operator()(const Vec<1>& point) const
    {
      if (!fixedtime)
      {
        IntegrationPoint ip(point(0));
        MappedIntegrationPoint<1,D> mip(ip, eltrans);
        return eval.Eval(& mip.GetPoint()(0));
      }
      else
      {
        IntegrationPoint ip(point(0));
        MappedIntegrationPoint<1,D> mip(ip, eltrans);
        Vec<D+1> p; 
        for (int d = 0; d < D; ++d)
          p(d) = mip.GetPoint()(d);
        p(D) = fixedtime;
        return eval.Eval(&p(0));
      }
    }

    virtual double operator()(const Vec<2>& point) const
    {
      if (!fixedtime)
      {
        IntegrationPoint ip(point(0),point(1));
        MappedIntegrationPoint<2,D> mip(ip, eltrans);
        return eval.Eval(& mip.GetPoint()(0));
      }
      else
      {
        IntegrationPoint ip(point(0),point(1));
        MappedIntegrationPoint<2,D> mip(ip, eltrans);
        Vec<D+1> p; 
        for (int d = 0; d < D; ++d)
          p(d) = mip.GetPoint()(d);
        p(D) = fixedtime;
        return eval.Eval(&p(0));
      }
    }

    virtual double operator()(const Vec<3>& point) const
    {
      if (!fixedtime)
      {
        IntegrationPoint ip(point(0),point(1),point(2));
        MappedIntegrationPoint<3,D> mip(ip, eltrans);
        return eval.Eval(& mip.GetPoint()(0));
      }
      else
      {
        IntegrationPoint ip(point(0),point(1),point(2));
        MappedIntegrationPoint<3,D> mip(ip, eltrans);
        Vec<D+1> p; 
        for (int d = 0; d < D; ++d)
          p(d) = mip.GetPoint()(d);
        p(D) = fixedtime;
        return eval.Eval(&p(0));
      }
    }

  };

  template <int D> // D : resulting space dimension..
  class SpaceTimeEvalFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    const EvalFunction & eval;
    const ElementTransformation & eltrans;
    TimeInterval ti;
  public:
    SpaceTimeEvalFunctionEvaluator( const string & str, 
                                    const ElementTransformation & a_eltrans, 
                                    const TimeInterval & a_ti) 
      : eval(str), eltrans(a_eltrans), ti(a_ti) { ; }

    SpaceTimeEvalFunctionEvaluator( const EvalFunction & str, 
                                    const ElementTransformation & a_eltrans, 
                                    const TimeInterval & a_ti) 
      : eval(str), eltrans(a_eltrans), ti(a_ti) { ; }

    virtual double operator()(const Vec<1>& point) const
    {
      Vec<3> p(0.0); 
      p(0) = ( (1.0-point(0)) * ti.first + point(0) * ti.second );
      return eval.Eval(& p(0));
    }

    virtual double operator()(const Vec<2>& point) const
    {
      IntegrationPoint ip(point(0));
      MappedIntegrationPoint<1,D> mip(ip, eltrans);
      Vec<3> p = mip.GetPoint();
      p(D) = ( (1.0-point(1)) * ti.first + point(1) * ti.second );
      return eval.Eval(& p(0));
    }

    virtual double operator()(const Vec<3>& point) const
    {
      IntegrationPoint ip(point(0),point(1));
      MappedIntegrationPoint<2,D> mip(ip, eltrans);
      Vec<3> p = mip.GetPoint();
      p(D) = ( (1.0-point(2)) * ti.first + point(2) * ti.second );
      // std::cout << "p = " << p << std::endl;
      return eval.Eval(& p(0));
    }

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
      MappedIntegrationPoint<3,D> mip(ip, eltrans);
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
      mipp.Point()[D] = (1.0-point(2)) * ti.first + point(2) * ti.second;
      // std::cout << " mipp.Point() = " << mipp.Point() << std::endl;
      // std::cout << " eval.Evaluate(mipp) = " << eval.Evaluate(mipp) << std::endl;
      return eval.Evaluate(mipp);
    }

    virtual double operator()(const Vec<4>& point) const
    {
      IntegrationPoint ip(point(0),point(1),point(2));
      MappedIntegrationPoint<3,D> mip(ip, eltrans);
      DimMappedIntegrationPoint<D+1> mipp(ip,eltrans);
      mipp.Point().Range(0,D) = mip.GetPoint();
      mipp.Point()[D] = (1.0-point(3)) * ti.first + point(3) * ti.second;
      return eval.Evaluate(mipp);
    }


  };



} // end of namespace


#endif
