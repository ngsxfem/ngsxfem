#pragma once

#include <bla.hpp>
#include <comp.hpp>
// #include <python_ngstd.hpp>

namespace ngfem
{

  class BitArrayCoefficientFunction : public CoefficientFunction
  {
    shared_ptr<BitArray> ba;
  public:
    BitArrayCoefficientFunction (shared_ptr<BitArray> acf);
    virtual ~BitArrayCoefficientFunction ();
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;
  };

}
