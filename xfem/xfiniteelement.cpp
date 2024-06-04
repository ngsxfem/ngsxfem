/// from ngsxfem
#include "xfiniteelement.hpp"

namespace ngfem
{

  XDummyFE::XDummyFE (DOMAIN_TYPE a_sign, ELEMENT_TYPE a_et)
    : // FiniteElement(),
    sign(a_sign), et(a_et) { ndof = 0; order = 0; }

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns,
                                 Allocator & lh)
    : base(a_base),
    localsigns(a_localsigns.Size(),lh)
  {
    ndof = base.GetNDof();
    order = base.Order();
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };


  XFiniteElement::~XFiniteElement() {; };
  /// the name
  string XFiniteElement::ClassName(void) const {return "X-"+base.ClassName();};

  const FlatArray<DOMAIN_TYPE>& XFiniteElement::GetSignsOfDof() const
  {
    return localsigns;
  };



  void SFiniteElement::CalcShape (const IntegrationPoint & ip,
                                  BareSliceVector<> shape) const
  {
    Vec<2> x = ip.Point();
    Vec<2> d = x - cuts.Col(0);
    Vec<2> dref = cuts.Col(1) - cuts.Col(0);
    const double xhat = InnerProduct(d,dref)/(L2Norm(dref)*L2Norm(dref));
    IntegrationPoint ipref(xhat,0,0,0.0);
    basefe->CalcShape(ipref,shape);
  }

  SFiniteElement::SFiniteElement(Mat<2> acuts, int aorder, Allocator & lh)
    : cuts(acuts) //, order(aorder)
  {
    ndof = aorder+1;
    L2HighOrderFE<ET_SEGM> * hofe =  new (lh) L2HighOrderFE<ET_SEGM> ();
    hofe -> L2HighOrderFE<ET_SEGM>::SetOrder (aorder);
    hofe -> L2HighOrderFE<ET_SEGM>::ComputeNDof();
    basefe = hofe;
  };

  SFiniteElement::~SFiniteElement() {; };

  /// the name
  string SFiniteElement::ClassName(void) const {return "SFE";};

} // end of namespace
