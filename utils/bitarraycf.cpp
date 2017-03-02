#include "bitarraycf.hpp"

namespace ngfem
{

  BitArrayCoefficientFunction::BitArrayCoefficientFunction (shared_ptr<BitArray> aba)
    : CoefficientFunction(1.0,false), ba(aba)
  {}

  BitArrayCoefficientFunction::~BitArrayCoefficientFunction ()
  {}

  double BitArrayCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    if (ba->Test(ip.GetTransformation().GetElementNr()))
      return 1.0;
    else
      return 0.0;
  }

  void BitArrayCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                FlatMatrix<double> values) const
  {
    if (ba->Test(ir.GetTransformation().GetElementNr()))
      values = 1.0;
    else
      values = 0.0;
  }

  void BitArrayCoefficientFunction :: Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    if (ba->Test(ir.GetTransformation().GetElementNr()))
      values.AddSize(Dimension(), ir.Size()) = 1.0;
    else
      values.AddSize(Dimension(), ir.Size()) = 0.0;
  }

}
