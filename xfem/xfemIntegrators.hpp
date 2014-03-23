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
  class XLaplaceIntegrator : public BilinearFormIntegrator
  {
    double alpha_neg;
    double alpha_pos;
  public:
    XLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
    { 
      alpha_neg = coeffs[0]->EvaluateConst();
      alpha_pos = coeffs[1]->EvaluateConst();
    }
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


/*
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
*/

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



  /* use xfem-implementation + "empty"-keyword
  template <int D>
  class NoXLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_lset;
    double alpha_neg;
    double alpha_pos;
  public:
    NoXLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
    { 
      alpha_neg = coeffs[0]->EvaluateConst();
      alpha_pos = coeffs[1]->EvaluateConst();
      coef_lset = coeffs[2];
    }
    virtual ~NoXLaplaceIntegrator(){ ; };

    virtual string Name () const { return "NoXLaplaceIntegrator"; }

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
  */




  template <int D>
  class FictXSourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    FictXSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~FictXSourceIntegrator(){ ; };

    virtual string Name () const { return "FictXSourceIntegrator"; }

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
  class FictXMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    FictXMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~FictXMassIntegrator(){ ; };

    virtual string Name () const { return "FictXMassIntegrator"; }

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
  class FictXLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    FictXLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~FictXLaplaceIntegrator(){ ; };

    virtual string Name () const { return "FictXLaplaceIntegrator"; }

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




}

#endif

