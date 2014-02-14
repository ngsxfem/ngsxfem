
/// from ngxfem
#include "xfiniteelement.hpp"

namespace ngfem
{

  template< class FE >
  XDummyFE<FE>::XDummyFE (DOMAIN_TYPE a_sign)
    : FE(), sign(a_sign) { ; }

  template< class FE >
  XFiniteElement<FE>::XFiniteElement(const FE & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                     const XLocalGeometryInformation& a_localgeom)
    : base(a_base), localsigns(a_localsigns), localgeom(a_localgeom)
  { ; };


  template< class FE >
  XFiniteElement<FE>::~XFiniteElement() { ; };

  /// the name
  template< class FE >
  string XFiniteElement<FE>::ClassName(void) const {return "X-"+base.ClassName();};

  template< class FE >
  const Array<DOMAIN_TYPE>& XFiniteElement<FE>::GetSignsOfDof() const  
  {
    return localsigns;
  };

  template< class FE >
  const XLocalGeometryInformation& XFiniteElement<FE>::GetLocalGeometry() const
  {
    return localgeom;
  };

  // template class XDummyFE<1>;
  // template class XDummyFE<2>;
  // template class XDummyFE<3>;

  // template class XFiniteElement<1>;
  // template class XFiniteElement<2>;
  // template class XFiniteElement<3>;


} // end of namespace
