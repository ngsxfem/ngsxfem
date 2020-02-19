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
}
