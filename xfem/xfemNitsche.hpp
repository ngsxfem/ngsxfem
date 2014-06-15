#ifndef FILE_XFEMNITSCHE_HPP
#define FILE_XFEMNITSCHE_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"
#include "xfemIntegrators.hpp"
#include "stxfemIntegrators.hpp"

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
    CoefficientFunction * ab_neg;
    CoefficientFunction * ab_pos;
    XLaplaceIntegrator<D> * ablockintegrator;
    bool minimal_stabilization;
  public:
    XNitscheIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        lambda(coeffs.Size() > 4 ? coeffs[4] : NULL)
    { 
       minimal_stabilization = coeffs.Size() <= 4;

      const double abn = alpha_neg->EvaluateConst() * beta_neg->EvaluateConst(); 
      const double abp = alpha_pos->EvaluateConst() * beta_pos->EvaluateConst(); 
      ab_neg = new ConstantCoefficientFunction(abn);
      ab_pos = new ConstantCoefficientFunction(abp);
      Array<CoefficientFunction * > lapcoeffs(2);
      lapcoeffs[0] = ab_neg;
      lapcoeffs[1] = ab_pos;
      ablockintegrator = new XLaplaceIntegrator<D>(lapcoeffs);
    }

    virtual ~XNitscheIntegrator()
    { 
      if (ab_neg) delete ab_neg;
      if (ab_pos) delete ab_pos;
      if (ablockintegrator) delete ablockintegrator; 
    }

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
  class XNitscheRhsJumpIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * lambda;
    CoefficientFunction * ab_neg;
    CoefficientFunction * ab_pos;
    CoefficientFunction * coef_rhs;
    XLaplaceIntegrator<D> * ablockintegrator;
    bool minimal_stabilization;
  public:
    XNitscheRhsJumpIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        coef_rhs(coeffs[4]),
        lambda(coeffs.Size() > 5 ? coeffs[5] : NULL)
    { 
       minimal_stabilization = coeffs.Size() <= 4;

      const double abn = alpha_neg->EvaluateConst() * beta_neg->EvaluateConst(); 
      const double abp = alpha_pos->EvaluateConst() * beta_pos->EvaluateConst(); 
      ab_neg = new ConstantCoefficientFunction(abn);
      ab_pos = new ConstantCoefficientFunction(abp);
      Array<CoefficientFunction * > lapcoeffs(2);
      lapcoeffs[0] = ab_neg;
      lapcoeffs[1] = ab_pos;
      ablockintegrator = new XLaplaceIntegrator<D>(lapcoeffs);
    }

    virtual ~XNitscheRhsJumpIntegrator()
    { 
      if (ab_neg) delete ab_neg;
      if (ab_pos) delete ab_pos;
      if (ablockintegrator) delete ablockintegrator; 
    }

    virtual string Name () const { return "XNitscheRhsJumpIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> & elvec,
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

    double t1;
    double t0;
    double tau;

    CoefficientFunction * ab_neg;
    CoefficientFunction * ab_pos;
    CoefficientFunction * t_old;
    CoefficientFunction * t_new;

    SpaceTimeXLaplaceIntegrator<D> * ablockintegrator;
    bool minimal_stabilization;

  public:
    SpaceTimeXNitscheIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        lambda(coeffs.Size() > 6 ? coeffs[4] : NULL),
        t0(coeffs.Size() > 6 ? coeffs[5]->EvaluateConst() : coeffs[4]->EvaluateConst()),
        t1(coeffs.Size() > 6 ? coeffs[6]->EvaluateConst() : coeffs[5]->EvaluateConst())
    { 

      minimal_stabilization = coeffs.Size() <= 6;

      const double abn = alpha_neg->EvaluateConst() * beta_neg->EvaluateConst(); 
      const double abp = alpha_pos->EvaluateConst() * beta_pos->EvaluateConst(); 

      ab_neg = new ConstantCoefficientFunction(abn);
      ab_pos = new ConstantCoefficientFunction(abp);
      t_old = new ConstantCoefficientFunction(t0);
      t_new = new ConstantCoefficientFunction(t1);

      Array<CoefficientFunction * > lapcoeffs(4);

      lapcoeffs[0] = ab_neg;
      lapcoeffs[1] = ab_pos;
      lapcoeffs[2] = t_old;
      lapcoeffs[3] = t_new;
      ablockintegrator = new SpaceTimeXLaplaceIntegrator<D>(lapcoeffs);

      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
    }
    virtual ~SpaceTimeXNitscheIntegrator()
    { 
      if (ab_neg) delete ab_neg;
      if (ab_pos) delete ab_pos;
      if (ablockintegrator) delete ablockintegrator; 
    }

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

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }
  };


  template <int D, NITSCHE_VARIANTS::KAPPA_CHOICE kappa_choice>
  class FictXNitscheIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * alpha_neg;
    CoefficientFunction * alpha_pos;
    CoefficientFunction * beta_neg;
    CoefficientFunction * beta_pos;
    CoefficientFunction * lambda;
    CoefficientFunction * ab_neg;
    CoefficientFunction * ab_pos;
    XLaplaceIntegrator<D> * ablockintegrator;
    bool minimal_stabilization;
  public:
    FictXNitscheIntegrator (const Array<CoefficientFunction*> & coeffs)
      : alpha_neg(coeffs[0]),alpha_pos(coeffs[1]), 
        beta_neg(coeffs[2]),beta_pos(coeffs[3]), 
        lambda(coeffs.Size() > 4 ? coeffs[4] : NULL)
    { 
       minimal_stabilization = coeffs.Size() <= 4;

      const double abn = alpha_neg->EvaluateConst() * beta_neg->EvaluateConst(); 
      const double abp = alpha_pos->EvaluateConst() * beta_pos->EvaluateConst(); 
      ab_neg = new ConstantCoefficientFunction(abn);
      ab_pos = new ConstantCoefficientFunction(abp);
      Array<CoefficientFunction * > lapcoeffs(2);
      lapcoeffs[0] = ab_neg;
      lapcoeffs[1] = ab_pos;
      ablockintegrator = new XLaplaceIntegrator<D>(lapcoeffs);
    }

    virtual ~FictXNitscheIntegrator()
    { 
      if (ab_neg) delete ab_neg;
      if (ab_pos) delete ab_pos;
      if (ablockintegrator) delete ablockintegrator; 
    }

    virtual string Name () const { return "FictXNitscheIntegrator"; }

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

