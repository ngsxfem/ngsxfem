#include "timecf.hpp"

namespace ngfem
{
  TimeVariableCoefficientFunction::TimeVariableCoefficientFunction ()
    : CoefficientFunction(1)
  {; }

  ///
  double TimeVariableCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    return mip.IP()(2);
  }

  double TimeVariableCoefficientFunction::EvaluateConst () const
  {
    throw Exception("not constant");
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
