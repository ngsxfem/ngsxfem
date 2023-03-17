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


    virtual shared_ptr<BilinearFormIntegrator> MakeBilinearFormIntegrator() const;
    virtual shared_ptr<LinearFormIntegrator> MakeLinearFormIntegrator() const;
  };



  class CutDifferentialSymbol : public DifferentialSymbol
  {
  public:
    shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;
    double scale = 1;
    
    CutDifferentialSymbol () : DifferentialSymbol(VOL) { ; }
    CutDifferentialSymbol (VorB _vb) : DifferentialSymbol(_vb) { ; }
    CutDifferentialSymbol (shared_ptr<LevelsetIntegrationDomain> _lsetdom, VorB _vb, VorB _element_vb, bool _skeleton)
      : DifferentialSymbol(_vb, _element_vb, _skeleton, 0), lsetintdom(_lsetdom) { ; }
    CutDifferentialSymbol (CutDifferentialSymbol & _cds, double _scale)
      : DifferentialSymbol(_cds), lsetintdom(_cds.lsetintdom), scale(_cds.scale*_scale) { ; }

    virtual ~CutDifferentialSymbol() { ; }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      if (lsetintdom)
      {
        if (scale != 1.0)
          return make_shared<CutIntegral> (scale * cf, make_shared<CutDifferentialSymbol>(*this));
        else
          return make_shared<CutIntegral> (cf, make_shared<CutDifferentialSymbol>(*this));
      }
      else
        throw Exception("no level set domain prescribed. Cannot define a CutIntegral.");
    }


  };

  class FacetPatchDifferentialSymbol;
  class FacetPatchIntegral : public Integral
  {
  public:
    int time_order;
    optional<double> tref = nullopt;
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

    virtual shared_ptr<BilinearFormIntegrator> MakeBilinearFormIntegrator() const;
    virtual shared_ptr<LinearFormIntegrator> MakeLinearFormIntegrator() const
    {
      throw Exception("MakeLinearFormIntegrator not Implemented for FacetPatchIntegral");
    }
  };



  class FacetPatchDifferentialSymbol : public DifferentialSymbol
  {
  public:
    int time_order;
    double scale = 1;
    optional<double> tref = nullopt;
    FacetPatchDifferentialSymbol (VorB _vb) : DifferentialSymbol(_vb) { ; }
    FacetPatchDifferentialSymbol (VorB _vb, VorB _element_vb, bool _skeleton, int _time_order, optional<double> _tref)
      : DifferentialSymbol(_vb, _element_vb, _skeleton, 0), time_order(_time_order), tref(_tref){ ; }
    FacetPatchDifferentialSymbol (FacetPatchDifferentialSymbol & _cds, double _scale)
      : DifferentialSymbol(_cds), time_order(_cds.time_order), scale(_cds.scale*_scale), tref(_cds.tref) { ; }

    virtual ~FacetPatchDifferentialSymbol() { ; }
    virtual shared_ptr<Integral> MakeIntegral(shared_ptr<CoefficientFunction> cf) const
    {
      if (scale!=1.0)
        return make_shared<FacetPatchIntegral> (scale * cf, make_shared<FacetPatchDifferentialSymbol>(*this));
      else
        return make_shared<FacetPatchIntegral> (cf, make_shared<FacetPatchDifferentialSymbol>(*this));
    }


  };
}