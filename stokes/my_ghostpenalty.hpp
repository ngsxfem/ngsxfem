#ifndef FILE_GHOSTPENALTY_HPP
#define FILE_GHOSTPENALTY_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
//#include "../spacetime/spacetimefe.hpp"
//#include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{

  template <int D>
  class MyGhostPenaltyIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    double gamma = 1.0;
  public:
    MyGhostPenaltyIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetBilinearFormIntegrator(coeffs)
    { 
      gamma = coeffs[0]->EvaluateConst();
    }

    virtual ~MyGhostPenaltyIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual bool IsSymmetric () const 
    { return true; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new MyGhostPenaltyIntegrator (coeffs);
    }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("MyGhostPenaltyIntegrator::CalcElementMatrix - not implemented!");
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
