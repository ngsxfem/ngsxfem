#pragma once

#include <fem.hpp>
#include "../utils/ngsxstd.hpp"

using ngcomp::CoefficientFunction;

namespace ngfem
{
  /// The coefficient evaluates random between two bounds (default: 0 and 1) pointwise
  class TimeVariableCoefficientFunction : public CoefficientFunction
  {
  private:
    double time;
    bool time_is_fixed;
  public:
    ///
    TimeVariableCoefficientFunction ();
    ///
    virtual ~TimeVariableCoefficientFunction () {}
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    virtual double EvaluateConst () const;
//    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
//    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;
//    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
//                           FlatArray<AFlatMatrix<double>*> input,
//                           AFlatMatrix<double> values) const;
    virtual void PrintReport (ostream & ost) const;
    virtual void FixTime(double t);
    virtual void UnfixTime();

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const;
    
  };

  /// Evaluates a different CoefficientFunction at a fixed time
  class FixTimeCoefficientFunction : public CoefficientFunction
  {
  private:
    shared_ptr<CoefficientFunction> coef;
    shared_ptr<ParameterCoefficientFunction<double>> time;
  public:
    ///
    FixTimeCoefficientFunction (shared_ptr<CoefficientFunction>, shared_ptr<ParameterCoefficientFunction<double>> );
    ///
    virtual ~FixTimeCoefficientFunction () {}
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                          FlatVector<> result) const;
    virtual double EvaluateConst () const;
//    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
//    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;
//    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
//                           FlatArray<AFlatMatrix<double>*> input,
//                           AFlatMatrix<double> values) const;
    virtual void PrintReport (ostream & ost) const;

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const;
    
  };
  
}
