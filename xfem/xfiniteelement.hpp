#ifndef FILE_XFINITEELEMENT_HPP
#define FILE_XFINITEELEMENT_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

/// from ngxfem
#include "../cutint/xintegration.hpp"  // for MasterElement
// #include "xfemIntegrators.hpp"



//using namespace ngsolve;
using namespace ngcomp;
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
    FlatXLocalGeometryInformation fxgeom;
    bool empty = false;
  public:
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns, 
                   shared_ptr<XLocalGeometryInformation> a_localgeom,
                   LocalHeap & lh);
    virtual ~XFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    const FiniteElement & GetBaseFE() const { return base; };

    const FlatArray<DOMAIN_TYPE>& GetSignsOfDof() const; 

    bool HasFlatLocalGeometry() const 
    {
      return !fxgeom.empty;
    }

    const FlatXLocalGeometryInformation & GetFlatLocalGeometry() const 
    { 
      if (fxgeom.empty)
        throw Exception(" no geometry ");
      else
        return fxgeom;
    } 

    virtual ELEMENT_TYPE ElementType() const { return base.ElementType(); }

    void SetEmpty(bool se = true){ empty = se; if (se) ndof = 0;}
    bool Empty() const{ return empty;}

  };

  /**
     surrounds a FiniteElement and adds information about signs of dofs and local geometry
   */
  class SFiniteElement : public ScalarFiniteElement<2>
  {
  protected:
    Mat<2> cuts;
    BaseScalarFiniteElement * basefe;
  public:
    SFiniteElement(Mat<2> acuts, int order, LocalHeap & lh);
    virtual ~SFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    /// compute shape, row is shape nr, col is ip nr
    virtual void CalcShape (const IntegrationPoint & ip,
                            BareSliceVector<> shape) const;

    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const
    {
      throw Exception("noenoe, ich soll nich");
    }

  };

} // end of namespace

#endif
