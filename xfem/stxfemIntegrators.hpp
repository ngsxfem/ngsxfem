#ifndef FILE_XFEMINTEGRATORS_HPP
#define FILE_XFEMINTEGRATORS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"

namespace ngfem
{
  template <int D>
  class SpaceTimeXMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    CoefficientFunction * coef_told;
    CoefficientFunction * coef_tnew;
  public:
    SpaceTimeXMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]),coef_told(coeffs[2]),coef_tnew(coeffs[3]) { ; }
    virtual ~SpaceTimeXMassIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXMassIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class SpaceTimeXLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    CoefficientFunction * coef_told;
    CoefficientFunction * coef_tnew;
  public:
    SpaceTimeXLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]),coef_told(coeffs[2]),coef_tnew(coeffs[3]) { ; }
    virtual ~SpaceTimeXLaplaceIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXLaplaceIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

  };


  template <int D>
  class SpaceTimeXSourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    CoefficientFunction * coef_told;
    CoefficientFunction * coef_tnew;
  public:
    SpaceTimeXSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]),coef_told(coeffs[2]),coef_tnew(coeffs[3]) { ; }
    virtual ~SpaceTimeXSourceIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXSourceIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> & elvec,
                       LocalHeap & lh) const;

  };


/*
  template <int D>
  class XLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XLaplaceIntegrator(){ ; };

    virtual string Name () const { return "XLaplaceIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XConvectionIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XConvectionIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XConvectionIntegrator(){ ; };

    virtual string Name () const { return "XConvectionIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XRobinIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XRobinIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1])
      { ; }
    virtual ~XRobinIntegrator(){ ; };

    virtual string Name () const { return "XRobinIntegrator"; }

    virtual int DimElement () const { return D-1; }
    virtual int DimSpace () const { return D; }
    // it is a boundary integral
    virtual bool BoundaryForm () const { return true; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

  };


  template <int D>
  class XNeumannIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XNeumannIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XNeumannIntegrator(){ ; };

    virtual string Name () const { return "XNeumannIntegrator"; }

    virtual int DimElement () const { return D-1; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return true; }


    // Calculates the element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> & elvec,
                       LocalHeap & lh) const;

  };
*/

}

#endif

