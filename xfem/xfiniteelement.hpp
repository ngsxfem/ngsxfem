#ifndef FILE_XFINITEELEMENT_HPP
#define FILE_XFINITEELEMENT_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

/// from ngxfem
#include "../cutint/xintegration.hpp"  // for MasterElement
// #include "xfemIntegrators.hpp"



using namespace ngsolve;
using namespace xintegration;

namespace ngfem
{
  /**
     a placeholder finite element
   */
  class XDummyFE : public FiniteElement
  {
  protected:
    const DOMAIN_TYPE sign;
    const ELEMENT_TYPE et;
  public:
    XDummyFE (DOMAIN_TYPE a_sign, ELEMENT_TYPE et);
    DOMAIN_TYPE GetDomainType() const { return sign;}
    virtual ELEMENT_TYPE ElementType() const { return et; }
  };

  /**
     a placeholder finite element
   */
  class LevelsetContainerFE : public FiniteElement
  {
  protected:
    const CoefficientFunction * coef_lset = NULL;
  public:
    double tnew;
    double told;
    LevelsetContainerFE (const CoefficientFunction *, double ta = 0.0, double tb = 0.0);
    virtual ELEMENT_TYPE ElementType() const { return ET_POINT; }
    const CoefficientFunction * GetLevelsetCoefficient() const { return coef_lset; }
  };


  /**
     surrounds a FiniteElement and adds information about signs of dofs and local geometry
   */
  class XFiniteElement : public FiniteElement
  {
  protected:
    const FiniteElement & base;
    const FlatArray<DOMAIN_TYPE> localsigns;
    FlatXLocalGeometryInformation fxgeom;
    FlatXLocalGeometryInformation fxgeom_downtrace;
    FlatXLocalGeometryInformation fxgeom_uptrace;
    bool empty = false;
  public:
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   const XLocalGeometryInformation* a_localgeom,
                   const XLocalGeometryInformation* a_localgeom_downtrace,
                   const XLocalGeometryInformation* a_localgeom_uptrace,
                   LocalHeap & lh);
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   const XLocalGeometryInformation* a_localgeom,
                   const XLocalGeometryInformation* a_localgeom_downtrace,
                   LocalHeap & lh);
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   const XLocalGeometryInformation* a_localgeom,
                   LocalHeap & lh);
    virtual ~XFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    const FiniteElement & GetBaseFE() const { return base; };

    const FlatArray<DOMAIN_TYPE>& GetSignsOfDof() const; 

    const FlatXLocalGeometryInformation & GetFlatLocalGeometry() const 
    { 
      if (fxgeom.empty)
        throw Exception(" no geometry ");
      else
        return fxgeom;
    } 

    const FlatXLocalGeometryInformation & GetFlatLocalGeometryUpTrace() const 
    { 
      if (fxgeom_uptrace.empty)
        throw Exception(" no geometry ");
      else
        return fxgeom_uptrace;
    } 

    const FlatXLocalGeometryInformation & GetFlatLocalGeometryDownTrace() const 
    { 
      if (fxgeom_downtrace.empty)
        throw Exception(" no geometry ");
      else
        return fxgeom_downtrace;
    } 

    virtual ELEMENT_TYPE ElementType() const { return base.ElementType(); }

    void SetEmpty(bool se = true){ empty = se; if (se) ndof = 0;}
    bool Empty() const{ return empty;}

  };

} // end of namespace

#endif
