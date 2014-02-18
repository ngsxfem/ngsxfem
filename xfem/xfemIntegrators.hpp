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
  class XMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XMassIntegrator(){ ; };

    virtual string Name () const { return "XMassIntegrator"; }

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
  class XSourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    XSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XSourceIntegrator(){ ; };

    virtual string Name () const { return "XSourceIntegrator"; }

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


  template <int D>
  class SpaceTimeXMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    SpaceTimeXMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
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
  class SpaceTimeXSourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    SpaceTimeXSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
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

  namespace NITSCHE_VARIANTS{
    enum KAPPA_CHOICE{
      HALFHALF,
      HANSBO,
      BETA,
      ALPHA,
      ALPHABETA,
      HANSBOBETA
    };
  }


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  class XNitscheIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * lambda;
  public:
    XNitscheIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        lambda(coeffs[4])
      { ; }
    virtual ~XNitscheIntegrator(){ ; };

    virtual string Name () const { return "XNitscheIntegrator"; }

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

