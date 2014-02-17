
/// from ngxfem
#include "xfiniteelement.hpp"

namespace ngfem
{

  XDummyFE::XDummyFE (DOMAIN_TYPE a_sign, ELEMENT_TYPE a_et)
    : // FiniteElement(),
    sign(a_sign), et(a_et) { ndof = 0; }

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 const XLocalGeometryInformation& a_localgeom)
    : base(a_base), localsigns(a_localsigns), localgeom(a_localgeom)
  { ndof = base.GetNDof(); };


  XFiniteElement::~XFiniteElement() { ; };

  /// the name
  string XFiniteElement::ClassName(void) const {return "X-"+base.ClassName();};

  const Array<DOMAIN_TYPE>& XFiniteElement::GetSignsOfDof() const  
  {
    return localsigns;
  };

  const XLocalGeometryInformation& XFiniteElement::GetLocalGeometry() const
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
