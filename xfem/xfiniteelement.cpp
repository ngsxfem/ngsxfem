
/// from ngxfem
#include "xfiniteelement.hpp"

namespace ngfem
{

  XDummyFE::XDummyFE (DOMAIN_TYPE a_sign, ELEMENT_TYPE a_et)
    : // FiniteElement(),
    sign(a_sign), et(a_et) { ndof = 0; }

  LevelsetContainerFE::LevelsetContainerFE(shared_ptr<CoefficientFunction> coeflset, double ta, double tb)
    : // FiniteElement(),
    coef_lset(coeflset), told(ta), tnew(tb) {ndof = 0; }


  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 const XLocalGeometryInformation* a_localgeom,
                                 const XLocalGeometryInformation* a_localgeom_downtrace,
                                 const XLocalGeometryInformation* a_localgeom_uptrace,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh),
      fxgeom_downtrace(*a_localgeom_downtrace,lh),
      fxgeom_uptrace(*a_localgeom_uptrace,lh)
  { 
    ndof = base.GetNDof(); 
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 const XLocalGeometryInformation* a_localgeom,
                                 const XLocalGeometryInformation* a_localgeom_downtrace,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh),
      fxgeom_downtrace(*a_localgeom_downtrace,lh)
  { 
    ndof = base.GetNDof(); 
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 const XLocalGeometryInformation* a_localgeom,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh)
  { 
    ndof = base.GetNDof(); 
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };


  XFiniteElement::~XFiniteElement() { ; };

  /// the name
  string XFiniteElement::ClassName(void) const {return "X-"+base.ClassName();};

  const FlatArray<DOMAIN_TYPE>& XFiniteElement::GetSignsOfDof() const  
  {
    return localsigns;
  };

  // template class XDummyFE<1>;
  // template class XDummyFE<2>;
  // template class XDummyFE<3>;

  // template class XFiniteElement<1>;
  // template class XFiniteElement<2>;
  // template class XFiniteElement<3>;


} // end of namespace
