#ifndef FILE_XFINITEELEMENT_HPP
#define FILE_XFINITEELEMENT_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

/// from ngxfem
#include "../cutint/xintegration.hpp"  // for MasterElement

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
  public:
    XFiniteElement(const FiniteElement & a_base,
                   const Array<DOMAIN_TYPE>& a_localsigns,
                   Allocator & lh);
    virtual ~XFiniteElement();
    /// the name
    virtual string ClassName(void) const;

    const FiniteElement & GetBaseFE() const { return base; };

    const FlatArray<DOMAIN_TYPE>& GetSignsOfDof() const;

    virtual ELEMENT_TYPE ElementType() const { return base.ElementType(); }
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
    SFiniteElement(Mat<2> acuts, int order, Allocator & lh);
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
