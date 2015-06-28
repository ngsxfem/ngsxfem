#ifndef FILE_XFEMINTEGRATORS_HPP
#define FILE_XFEMINTEGRATORS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"

namespace ngfem
{

  template<int D>
  void CastXScalarFiniteElements (const FiniteElement & base_fel,
                                  const ScalarFiniteElement<D> * & scafe,
                                  const XFiniteElement * & xfe,
                                  const XDummyFE * & dummfe);               
  
  template <int D>
  class XMassIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    XMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~XMassIntegrator(){ ; };

    virtual string Name () const { return "XMassIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }
    virtual bool IsSymmetric () const { return true; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

  };


  template <int D>
  class XLaplaceIntegrator : public BilinearFormIntegrator
  {
    double alpha_neg;
    double alpha_pos;
  public:
    XLaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
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
    virtual bool IsSymmetric () const { return true; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XConvectionIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> conv_neg;
    shared_ptr<CoefficientFunction> conv_pos;
  public:
    XConvectionIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : conv_neg(coeffs[0]),conv_pos(coeffs[1]) { ; }

    virtual ~XConvectionIntegrator(){ ; };

    virtual string Name () const { return "XConvectionIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }
    virtual bool IsSymmetric () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XSourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    XSourceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
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
                       FlatVector<double> elvec,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XRobinIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    XRobinIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1])
      { ; }
    virtual ~XRobinIntegrator(){ ; };

    virtual string Name () const { return "XRobinIntegrator"; }

    virtual int DimElement () const { return D-1; }
    virtual int DimSpace () const { return D; }
    // it is a boundary integral
    virtual bool BoundaryForm () const { return true; }
    virtual bool IsSymmetric () const { return true; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

  };

  template <int D>
  class XNeumannIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    XNeumannIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
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
                       FlatVector<double> elvec,
                       LocalHeap & lh) const;

  };


}

#endif

