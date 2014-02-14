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
  template< class FE >
  class XDummyFE : public FE
  {
  protected:
    const DOMAIN_TYPE sign;
  public:
    XDummyFE (DOMAIN_TYPE a_sign);
    DOMAIN_TYPE GetDomainType() const { return sign;}
  };

  /**
     surrounds a FiniteElement and adds information about signs of dofs and local geometry
   */
  template< class FE >
  class XFiniteElement : public FiniteElement
  {
  protected:
    const FE& base;
    const Array<DOMAIN_TYPE> localsigns;
    const XLocalGeometryInformation& localgeom;
  public:
    XFiniteElement(const FE& a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   const XLocalGeometryInformation& a_localgeom);
    virtual ~XFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    const Array<DOMAIN_TYPE>& GetSignsOfDof() const; 

    const XLocalGeometryInformation& GetLocalGeometry() const; 
  };

} // end of namespace

#endif
