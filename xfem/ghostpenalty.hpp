#ifndef FILE_GHOSTPENALTY_HPP
#define FILE_GHOSTPENALTY_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"
#include "../spacetime/spacetimefe.hpp"
#include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{


  template <int D>
  class LowOrderGhostPenaltyIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_lam_neg;
    CoefficientFunction *coef_lam_pos;
    double told = 0.0;
    double tnew = 1.0;
    double tau = 1.0;
    double delta = 1.0;
  public:
    LowOrderGhostPenaltyIntegrator (Array<CoefficientFunction*> & coeffs) 
      : FacetBilinearFormIntegrator(coeffs)
    { 
      coef_lam_neg  = coeffs[0];
      coef_lam_pos  = coeffs[1];
      if (coeffs.Size()>3)
      {
        told = coeffs[2]->EvaluateConst();
        tnew = coeffs[3]->EvaluateConst();
        tau = tnew-told;
        delta = coeffs[4]->EvaluateConst();
      }
      else
        delta = coeffs[2]->EvaluateConst();
    }

    virtual ~LowOrderGhostPenaltyIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 0; }
    
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new LowOrderGhostPenaltyIntegrator (coeffs);
    }

    virtual void SetTimeInterval (const TimeInterval & ti)
    { told = ti.first; tnew = ti.second; tau = tnew-told; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("LowOrderGhostPenaltyIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                  const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                  const FiniteElement & volumefel2, int LocalFacetNr2,
                                  const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh// ,
                                  // BitArray* twice
                                 ) const;

  };


}

#endif
