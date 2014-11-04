#ifndef FILE_SDXFEMINTEGRATORS_HPP
#define FILE_SDXFEMINTEGRATORS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"

namespace ngfem
{

  template <int D>
  class SDXIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * conv_neg;
    CoefficientFunction * conv_pos;
    CoefficientFunction * mass_neg;
    CoefficientFunction * mass_pos;
  public:
    SDXIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : beta_neg(coeffs[0]),
        beta_pos(coeffs[1]),
        alpha_neg(coeffs[2]),
        alpha_pos(coeffs[3]),
        conv_neg(coeffs[4]),
        conv_pos(coeffs[5]),
        mass_neg(coeffs[6]),
        mass_pos(coeffs[7]) { ; }
    virtual ~SDXIntegrator(){ ; };

    virtual string Name () const { return "SDXIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }

    virtual void
    AddPointContribution (const MappedIntegrationPoint<D,D> & mip,
                          FlatVector<> & shape,
                          FlatMatrixFixWidth<D> & dshape,
                          FlatMatrixFixWidth<D*D> & ddshape_ref_h1,
                          FlatVector<> & lapshape,
                          FlatVector<> & dudwshape,
                          FlatVector<> & diffopshape,
                          FlatMatrix<double> & elmat, 
                          int ndof_h1,
                          int ndof_x,
                          int order,
                          DOMAIN_TYPE dt,
                          double convmax,
                          const ScalarFiniteElement<D>* scafe,
                          const XFiniteElement * xfe = NULL) const;

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

  };

/*
  template <int D>
  class XNitscheConvScaledIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * conv_neg;
    CoefficientFunction * conv_pos;
    CoefficientFunction * lambda;
  public:
    XNitscheConvScaledIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        conv_neg(coeffs[4]),conv_pos(coeffs[5]), 
        lambda(coeffs[6])
      { ; }
    virtual ~XNitscheConvScaledIntegrator(){ ; };

    virtual string Name () const { return "XNitsche(Convection scaled) Integrator"; }

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
  class SDXSourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * conv_neg;
    CoefficientFunction * conv_pos;
    CoefficientFunction * mass_neg;
    CoefficientFunction * mass_pos;
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    SDXSourceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : beta_neg(coeffs[0]),
        beta_pos(coeffs[1]),
        alpha_neg(coeffs[2]),
        alpha_pos(coeffs[3]),
        conv_neg(coeffs[4]),
        conv_pos(coeffs[5]),
        mass_neg(coeffs[6]),
        mass_pos(coeffs[7]),
        coef_neg(coeffs[8]),
        coef_pos(coeffs[9]) { ; }
    virtual ~SDXSourceIntegrator(){ ; };

    virtual string Name () const { return "SDXSourceIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    virtual void
    AddPointContribution (const MappedIntegrationPoint<D,D> & mip,
                          FlatMatrixFixWidth<D> & dshape,
                          FlatVector<> & dudwshape,
                          FlatVector<double> & elvec, 
                          int ndof_h1, int ndof_x,
                          int order,
                          DOMAIN_TYPE dt,
                          double convmax,
                          const ScalarFiniteElement<D>* scafe,
                          const XFiniteElement * xfe = NULL) const;

    // Calculates the element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> elvec,
                       LocalHeap & lh) const;

  };

}

#endif

