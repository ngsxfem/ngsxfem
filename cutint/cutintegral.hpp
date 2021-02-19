#pragma once

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
#include <integratorcf.hpp> // 
#include <variant> // 
#include <memory> // 

#include "../cutint/xintegration.hpp"
using namespace xintegration;

// #include "xfiniteelement.hpp"
// #include "../spacetime/spacetimefe.hpp"
// #include "../spacetime/spacetimeintegrators.hpp"
// #include "../utils/stcoeff.hpp"

namespace ngfem
{

  class CutDifferentialSymbol;
  class CutIntegral : public Integral
  {
  public:
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;
    CutIntegral (shared_ptr<CoefficientFunction> _cf,
                 shared_ptr<CutDifferentialSymbol> _dx);
    virtual ~CutIntegral() { }

    template <typename TSCAL>
    TSCAL T_CutIntegrate (const ngcomp::MeshAccess & ma,
                       FlatVector<TSCAL> element_wise);

    virtual double Integrate (const ngcomp::MeshAccess & ma,
                              FlatVector<double> element_wise);

    virtual Complex Integrate (const ngcomp::MeshAccess & ma,
                               FlatVector<Complex> element_wise);


    virtual shared_ptr<BilinearFormIntegrator> MakeBilinearFormIntegrator();
    virtual shared_ptr<LinearFormIntegrator> MakeLinearFormIntegrator();
  };



  class CutDifferentialSymbol : public DifferentialSymbol
  {
  public:
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;
    
    CutDifferentialSymbol (VorB _vb) : DifferentialSymbol(_vb) { ; }
    CutDifferentialSymbol (shared_ptr<LevelsetIntegrationDomain> _lsetdom, VorB _vb, VorB _element_vb, bool _skeleton)
      : DifferentialSymbol(_vb, _element_vb, _skeleton, 0), lsetintdom(_lsetdom) { ; }

    virtual ~CutDifferentialSymbol() { ; }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      if (lsetintdom)
        return make_shared<CutIntegral> (cf, make_shared<CutDifferentialSymbol>(*this));
      else
        throw Exception("no level set domain prescribed. Cannot define a CutIntegral.");
    }


  };

  class FacetPatchDifferentialSymbol;
  class FacetPatchIntegral : public Integral
  {
  public:
    using Integral::dx;
    int time_order;
    FacetPatchIntegral (shared_ptr<CoefficientFunction> _cf,
                 shared_ptr<FacetPatchDifferentialSymbol> _dx);
    virtual ~FacetPatchIntegral() { }

    virtual double Integrate (const ngcomp::MeshAccess & ma,
                              FlatVector<double> element_wise)
    {
      throw Exception("Integrate not Implemented for FacetPatchIntegral");
    }

    virtual Complex Integrate (const ngcomp::MeshAccess & ma,
                               FlatVector<Complex> element_wise)
    {
      throw Exception("Integrate not Implemented for FacetPatchIntegral");
    }

    virtual shared_ptr<BilinearFormIntegrator> MakeBilinearFormIntegrator();
    virtual shared_ptr<LinearFormIntegrator> MakeLinearFormIntegrator()
    {
      throw Exception("MakeLinearFormIntegrator not Implemented for FacetPatchIntegral");
    }
  };



  class FacetPatchDifferentialSymbol : public DifferentialSymbol
  {
  public:
    int time_order;
    FacetPatchDifferentialSymbol (VorB _vb) : DifferentialSymbol(_vb) { ; }
    FacetPatchDifferentialSymbol (VorB _vb, VorB _element_vb, bool _skeleton, int _time_order)
      : DifferentialSymbol(_vb, _element_vb, _skeleton, 0), time_order(_time_order){ ; }

    virtual ~FacetPatchDifferentialSymbol() { ; }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      return make_shared<FacetPatchIntegral> (cf, make_shared<FacetPatchDifferentialSymbol>(*this));
    }


  };
}