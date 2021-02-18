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
    // using Integral::cf;
    // using Integral::dx;
    shared_ptr<CutDifferentialSymbol> cdx;
    CutIntegral (shared_ptr<CoefficientFunction> _cf,
                 shared_ptr<CutDifferentialSymbol> _dx)
      : Integral(_cf, VOL), cdx(_dx) { ; }
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
    VorB vb;
    VorB element_vb = VOL;
    bool skeleton = false;
    optional<variant<BitArray,string>> definedon;
    int bonus_intorder = 0;
    shared_ptr<ngcomp::GridFunction> deformation;
    std::map<ELEMENT_TYPE,shared_ptr<IntegrationRule>> userdefined_intrules;
    shared_ptr<BitArray> definedonelements;
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;
    
    CutDifferentialSymbol (VorB _vb) : DifferentialSymbol(_vb) { ; }
    CutDifferentialSymbol (shared_ptr<LevelsetIntegrationDomain> _lsetdom, VorB _vb, VorB _element_vb, bool _skeleton)
      : DifferentialSymbol(_vb, _element_vb, _skeleton, 0), lsetintdom(_lsetdom) { ; }

    virtual ~CutDifferentialSymbol() { }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      if (lsetintdom)
        return make_shared<CutIntegral> (cf, make_shared<CutDifferentialSymbol>(*this));
      else
        throw Exception("no level set domain prescribed. Cannot define a CutIntegral.");
    }


  };
}