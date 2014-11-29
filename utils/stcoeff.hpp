#ifndef FILE_STCOEFF_HPP
#define FILE_STCOEFF_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include <coefficient.hpp>
//#include "../spacetime/spacetimefe.hpp"

//using namespace ngsolve;
using namespace ngcomp;

namespace ngfem
{

  // template <int D>
  // inline bool MakeHigherDimensionCoefficientFunction(shared_ptr<CoefficientFunction> coef_in, const CoefficientFunction *& coef_out)
  // {

  //   //TODO: not necessary any more: DomainVariableCoefficientFunction not a template any more
  //   DomainVariableCoefficientFunction * coef = dynamic_cast<DomainVariableCoefficientFunction * > (coef_in);
  //   if (coef != NULL)
  //   {
  //     int numreg = coef->NumRegions();
  //     if (numreg == INT_MAX) numreg = 1;
  //     Array< shared_ptr<EvalFunction> > evals;
  //     evals.SetSize(numreg);
  //     for (int i = 0; i < numreg; ++i)
  //     {
  //       evals[i] = coef->GetEvalFunction(i);
  //     }
  //     coef_out = new DomainVariableCoefficientFunction(evals); 
  //   }
  //   else
  //     coef_out = coef_in;
  //   return (coef!=NULL);
  // }


} // end of namespace


#endif
