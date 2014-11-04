
#include "stcoeff.hpp"

using namespace ngsolve;

namespace ngfem
{

  template <int D>
  SpaceTimeDomainVariableCoefficientFunction<D>::SpaceTimeDomainVariableCoefficientFunction( const Array<EvalFunction*> & afun )
    : DomainVariableCoefficientFunction<D>(afun)
  {
    for (int i = 0; i < fun.Size(); ++i)
      if (fun[i])
        fun[i]->DefineArgument ("t", 3);
    numarg = 4;
  }

  template <int D>
  void SpaceTimeDomainVariableCoefficientFunction<D> ::
  Evaluate(const BaseMappedIntegrationPoint & ip,
           FlatVector<> result) const
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (fun.Size() == 1) elind = 0;

    if (! fun[elind] -> IsComplex ())
    {
      VectorMem<10> args(numarg);
      args.Range(0,D+1) = static_cast<const DimMappedIntegrationPoint<D+1>&>(ip).GetPoint();
	
      for (int i = 0, an = 3; i < depends_on.Size(); i++)
	  {
	    int dim = depends_on[i]->Dimension();
	    depends_on[i] -> Evaluate (ip, args.Range(an,an+dim));
	    an += dim;
	  }
      fun[elind]->Eval (&args(0), &result(0), result.Size());      
    }
    else
    {
      VectorMem<10, Complex> args(numarg);
      args.Range(0,D+1) = static_cast<const DimMappedIntegrationPoint<D+1>&>(ip).GetPoint();
	
      for (int i = 0, an = 3; i < depends_on.Size(); i++)
	  {
	    int dim = depends_on[i]->Dimension();
	    depends_on[i] -> Evaluate (ip, args.Range(an,an+dim));
	    an += dim;
	  }
      fun[elind]->Eval (&args(0), &result(0), result.Size());      
    }
  }



  template class SpaceTimeDomainVariableCoefficientFunction<1>;
  template class SpaceTimeDomainVariableCoefficientFunction<2>;
  template class SpaceTimeDomainVariableCoefficientFunction<3>;


} // end of namespace
