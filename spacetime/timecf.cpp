#include "timecf.hpp"

namespace ngfem
{
  TimeVariableCoefficientFunction::TimeVariableCoefficientFunction ()
    : CoefficientFunction(1)
  { time_is_fixed = false ; }

  ///
  double TimeVariableCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    //if(abs(mip.IP()(2) - mip.IP().Weight() ) > 1e-8) cout << "TimeVariableCoefficientFunction::Evaluate IP:" << mip.IP() << endl;
    if(time_is_fixed) return time;
    else {
      if (IsSpaceTimeIntegrationPoint(mip.IP()))
        return mip.IP().Weight();
      else
        throw Exception("TimeVariableCoefficientFunction::Evaluate called with a mere space IR");
        //return mip.IP()(2);
    }
  }

  double TimeVariableCoefficientFunction::EvaluateConst () const
  {
    throw Exception("not constant");
  }

  void TimeVariableCoefficientFunction::FixTime(double t){
      time = t; time_is_fixed = true;
  }

  void TimeVariableCoefficientFunction::UnfixTime(){
      time_is_fixed = false;
  }


  shared_ptr<CoefficientFunction> TimeVariableCoefficientFunction ::
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    if (var == this)
      return dir;
    else
      return make_shared<ConstantCoefficientFunction> (0.0);
    // throw Exception(string("Deriv not implemented for CF ")+typeid(*this).name());
  }
  
/*
   void TimeVariableCoefficientFunction::Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
   {
    for (int i = 0; i < values.Height(); ++i)
      for (int j = 0; j < values.Width(); ++j)
        values(i,j) = unif(re);
   }

   void TimeVariableCoefficientFunction::Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
   {
    for (int i = 0; i < values.Height(); ++i)
      for (int j = 0; j < values.Width(); ++j)
        values(i,j) = unif(re);
   }

   void TimeVariableCoefficientFunction::Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                                            FlatArray<AFlatMatrix<double>*> input,
                                            AFlatMatrix<double> values) const
   {
    for (int i = 0; i < values.Height(); ++i)
      for (int j = 0; j < values.Width(); ++j)
        values(i,j) = unif(re);
   }
 */
  void TimeVariableCoefficientFunction::PrintReport (ostream & ost) const
  {
    ost << "CoefficientFunction for time (in space time quad. rules)" << endl;
  }

  FixTimeCoefficientFunction::FixTimeCoefficientFunction (shared_ptr<CoefficientFunction> coef_, 
                                                          shared_ptr<ParameterCoefficientFunction<double>> t)
    : CoefficientFunction(coef_->Dimension(),coef_->IsComplex()), coef(coef_), time(t)
  { 
    bool hasproxy = false;
    coef->TraverseTree ([&hasproxy] (CoefficientFunction & cf)
    {
      if (dynamic_cast<ProxyFunction*> (&cf))
        hasproxy = true;
    });
    if (hasproxy)
      throw Exception("FixTimeCoef is called on a CoefficientFunction that contains a ProxyFunction.\n\
It is suggested to do the following instead:\n\
  * Use fix_tref_proxy(..,tref) directly on the involved proxies (if possible) or\n\
  * fix the integration domain. E.g. dCut(..,tref=1) yields space-time instead of spatial integration\n\
    points with fixed (reference) time value tref=1.");
  }

  ///
  double FixTimeCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    IntegrationPoint ip(mip.IP());
    MarkAsSpaceTimeIntegrationPoint(ip);
    ip.SetWeight(time->GetValue());
    double ret = 0;
    Switch<3> (mip.DimSpace()-1, [&] (auto DIM1) {
      if (mip.GetTransformation().VB())
      {
        MappedIntegrationPoint<DIM1.value,DIM1.value+1,double> newmip(ip, mip.GetTransformation());
        ret = coef->Evaluate(newmip);
      }
      else
      {
        MappedIntegrationPoint<DIM1.value+1,DIM1.value+1,double> newmip(ip, mip.GetTransformation());
        ret = coef->Evaluate(newmip);
      }
    });
    return ret;
  }

  void FixTimeCoefficientFunction::Evaluate(const BaseMappedIntegrationPoint & mip,
                                            FlatVector<> result) const
  {
    IntegrationPoint ip(mip.IP());
    MarkAsSpaceTimeIntegrationPoint(ip);
    ip.SetWeight(time->GetValue());
    Switch<3> (mip.DimSpace()-1, [&] (auto DIM1) {
      if (mip.GetTransformation().VB())
      {
        MappedIntegrationPoint<DIM1.value,DIM1.value+1,double> newmip(ip, mip.GetTransformation());
        coef->Evaluate(newmip,result);
      }
      else
      {
        MappedIntegrationPoint<DIM1.value+1,DIM1.value+1,double> newmip(ip, mip.GetTransformation());
        coef->Evaluate(newmip,result);
      }
    });
  }

  
  double FixTimeCoefficientFunction::EvaluateConst () const
  {
    return coef->EvaluateConst();
  }


  shared_ptr<CoefficientFunction> FixTimeCoefficientFunction ::
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    if (dynamic_pointer_cast<TimeVariableCoefficientFunction>(dir) != nullptr)
      return make_shared<ConstantCoefficientFunction> (0.0);
    else 
      return coef->Diff(var,dir);
  }
  
  void FixTimeCoefficientFunction::PrintReport (ostream & ost) const
  {
    ost << "CoefficientFunction for fixing the time variable" << endl;
  }
  
}
