
/// from ngxfem
#include "xfiniteelement.hpp"

namespace ngfem
{

  XDummyFE::XDummyFE (DOMAIN_TYPE a_sign, ELEMENT_TYPE a_et)
    : // FiniteElement(),
    sign(a_sign), et(a_et) { ndof = 0; order = 0; }

  LevelsetContainerFE::LevelsetContainerFE(shared_ptr<CoefficientFunction> coeflset, double ta, double tb)
    : // FiniteElement(),
    coef_lset(coeflset), tnew(tb), told(ta) {ndof = 0; }


  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 shared_ptr<XLocalGeometryInformation> a_localgeom,
                                 shared_ptr<XLocalGeometryInformation> a_localgeom_downtrace,
                                 shared_ptr<XLocalGeometryInformation> a_localgeom_uptrace,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh),
      fxgeom_downtrace(*a_localgeom_downtrace,lh),
      fxgeom_uptrace(*a_localgeom_uptrace,lh)
  { 
    ndof = base.GetNDof(); 
    order = base.Order(); 
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 shared_ptr<XLocalGeometryInformation> a_localgeom,
                                 shared_ptr<XLocalGeometryInformation> a_localgeom_downtrace,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh),
      fxgeom_downtrace(*a_localgeom_downtrace,lh)
  { 
    ndof = base.GetNDof(); 
    order = base.Order(); 
    for (int l = 0; l < localsigns.Size(); ++l)
      localsigns[l] = a_localsigns[l];
  };

  XFiniteElement::XFiniteElement(const FiniteElement & a_base, const Array<DOMAIN_TYPE>& a_localsigns, 
                                 shared_ptr<XLocalGeometryInformation> a_localgeom,
                                 LocalHeap & lh)
    : base(a_base), 
      localsigns(a_localsigns.Size(),lh), 
      fxgeom(*a_localgeom,lh)
  { 
    ndof = base.GetNDof(); 
    order = base.Order(); 
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

  void SFiniteElement::CalcShape (const IntegrationPoint & ip, 
                                  BareSliceVector<> shape) const
  {
    Vec<2> x = ip.Point();

    Vec<2> d = x - cuts.Col(0);
    Vec<2> dref = cuts.Col(1) - cuts.Col(0);
    // cout << " InnerProduct(d,dref) = " << InnerProduct(d,dref) << endl;
    // cout << " L2Norm(dref) = " << L2Norm(dref) << endl;
    const double xhat = InnerProduct(d,dref)/(L2Norm(dref)*L2Norm(dref));
    IntegrationPoint ipref(xhat,0,0,0.0);
    // shape(0) = 0.93;
    // cout << "shape" << shape(0) << endl;
    // cout << "shape" << shape.Width() << endl;
    // cout << "shape" << shape.Height() << endl;
    basefe->CalcShape(ipref,shape);
    // cout << "shape" << shape(0) << endl;

    // cout << "x := " << x << endl;
    // cout << "a := " << cuts.Col(0) << endl;
    // cout << "b := " << cuts.Col(1) << endl;
    // cout << "d := " << d << endl;
    // cout << "dref := " << dref << endl;
    // cout << "xhat := " << xhat << endl;
    // getchar();

  }



  SFiniteElement::SFiniteElement(Mat<2> acuts, int aorder, LocalHeap & lh)
    : cuts(acuts) //, order(aorder)
  { 
    ndof = aorder+1 ;
    // basefe = new (lh) DGFiniteElement<1>();

    L2HighOrderFE<ET_SEGM> * hofe =  new (lh) L2HighOrderFE<ET_SEGM> ();
    // hofe -> SetVertexNumbers (ngel.vertices);
    hofe -> L2HighOrderFE<ET_SEGM>::SetOrder (aorder);
    hofe -> L2HighOrderFE<ET_SEGM>::ComputeNDof();
    basefe = hofe;
    // basefe->SetOrder(aorder);
  };

  SFiniteElement::~SFiniteElement() { ; };

  /// the name
  string SFiniteElement::ClassName(void) const {return "SFE";};




} // end of namespace
