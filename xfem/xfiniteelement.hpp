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
     surrounds a FiniteElement and adds information about signs of dofs and local geometry
   */
  class XFiniteElement : public FiniteElement
  {
  protected:
    const FiniteElement & base;
    const FlatArray<DOMAIN_TYPE> localsigns;
    const XLocalGeometryInformation* localgeom;
  public:
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   const XLocalGeometryInformation* a_localgeom,
                   LocalHeap & lh);
    virtual ~XFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    const FlatArray<DOMAIN_TYPE>& GetSignsOfDof() const; 

    const XLocalGeometryInformation * GetLocalGeometry() const; 

    virtual ELEMENT_TYPE ElementType() const { return base.ElementType(); }

  };

} // end of namespace

#endif
