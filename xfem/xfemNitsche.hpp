#ifndef FILE_XFEMNITSCHE_HPP
#define FILE_XFEMNITSCHE_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"

namespace ngfem
{

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


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  class SpaceTimeXNitscheIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * lambda;
    CoefficientFunction * t_old;
    CoefficientFunction * t_new;
  public:
    SpaceTimeXNitscheIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        lambda(coeffs[4]),
        t_old(coeffs[5]), t_new(coeffs[6])
      { ; }
    virtual ~SpaceTimeXNitscheIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXNitscheIntegrator"; }

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

