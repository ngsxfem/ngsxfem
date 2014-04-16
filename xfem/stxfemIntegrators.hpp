#ifndef FILE_SPACETIMEXFEMINTEGRATORS_HPP
#define FILE_SPACETIMEXFEMINTEGRATORS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"
#include "../utils/stcoeff.hpp"

namespace ngfem
{
  template <int D>
  class SpaceTimeXMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    double t1;
    double t0;
    double tau;

  public:
    SpaceTimeXMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
    }
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

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

  };

  template <int D>
  class SpaceTimeXLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    double t1;
    double t0;
    double tau;
  public:
    SpaceTimeXLaplaceIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    { 
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
    }
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

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

  };

  template <int D, bool adjoint = false>
  class SpaceTimeXConvectionIntegrator : public BilinearFormIntegrator
  {
    const CoefficientFunction * coef_neg;
    const CoefficientFunction * coef_pos;
    double t1;
    double t0;
    double tau;
    bool made_neg;
    bool made_pos;
  public:
    SpaceTimeXConvectionIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
      made_neg = MakeHigherDimensionCoefficientFunction<D>(coeffs[0], coef_neg);
      made_pos = MakeHigherDimensionCoefficientFunction<D>(coeffs[1], coef_pos);
    }
    virtual ~SpaceTimeXConvectionIntegrator(){ 
      if (made_neg) delete coef_neg;
      if (made_pos) delete coef_pos;
    };

    virtual string Name () const { return "SpaceTimeXConvectionIntegrator"; }

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


  template <int D, bool adjoint = false>
  class SpaceTimeXTimeDerivativeIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    SpaceTimeXTimeDerivativeIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~SpaceTimeXTimeDerivativeIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXTimeDerivativeIntegrator"; }

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

  template <int D, TIME t>
  class SpaceTimeXTraceMassIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
  public:
    SpaceTimeXTraceMassIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; }
    virtual ~SpaceTimeXTraceMassIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXTraceMassIntegrator"; }

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
    const CoefficientFunction * coef_neg;
    const CoefficientFunction * coef_pos;
    double t1;
    double t0;
    double tau;
    bool made_neg;
    bool made_pos;
  public:
    SpaceTimeXSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
      made_neg = MakeHigherDimensionCoefficientFunction<D>(coeffs[0], coef_neg);
      made_pos = MakeHigherDimensionCoefficientFunction<D>(coeffs[1], coef_pos);
    }
    virtual ~SpaceTimeXSourceIntegrator()
    { 
      if (made_neg) delete coef_neg;
      if (made_pos) delete coef_pos;
    };

    virtual string Name () const { return "SpaceTimeXSourceIntegrator"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

    // Calculates the element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> & elvec,
                       LocalHeap & lh) const;
  };

  template <int D, TIME t>
  class SpaceTimeXTraceSourceIntegrator : public LinearFormIntegrator
  {
    const CoefficientFunction * coef_neg;
    const CoefficientFunction * coef_pos;
    bool made_neg;
    bool made_pos;
    double scale_pos = 1.0;
    double scale_neg = 1.0;
    double time = 0.0;
  public:
    SpaceTimeXTraceSourceIntegrator (const Array<CoefficientFunction*> & coeffs)
    { 
      made_neg = MakeHigherDimensionCoefficientFunction<D>(coeffs[0], coef_neg);
      made_pos = MakeHigherDimensionCoefficientFunction<D>(coeffs[1], coef_pos);
      if (coeffs.Size() > 2)
        time = coeffs[2]->EvaluateConst();
    }
    virtual ~SpaceTimeXTraceSourceIntegrator(){
      ;
    }

    virtual string Name () const { return "SpaceTimeXTraceSourceIntegrator"; }

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

    void SetTime(double t_) { time = t_; }

    virtual void ChangeNegPosCoefficient(CoefficientFunction * neg, CoefficientFunction * pos, double dneg = 0.0, double dpos = 0.0)
    {
      if (made_neg) {delete coef_neg; made_neg = false;}
      if (made_pos) {delete coef_pos; made_pos = false;}

      made_neg = MakeHigherDimensionCoefficientFunction<D>(neg, coef_neg);
      made_pos = MakeHigherDimensionCoefficientFunction<D>(pos, coef_pos);
      // coef_neg = neg;
      // coef_pos = pos;
      scale_neg = dneg;
      scale_pos = dpos;
    }

  };


  template <int D>
  class SpaceTimeXRobinIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    double t1;
    double t0;
    double tau;

  public:
    SpaceTimeXRobinIntegrator (const Array<CoefficientFunction*> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) 
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;
    }
    virtual ~SpaceTimeXRobinIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXRobinIntegrator"; }

    virtual int DimElement () const { return D-1; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return true; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> & elmat,
                       LocalHeap & lh) const;

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

  };

  template <int D>
  class SpaceTimeXNeumannIntegrator : public LinearFormIntegrator
  {
    const CoefficientFunction * coef_neg;
    const CoefficientFunction * coef_pos;
    
    double t1;
    double t0;
    double tau;

  public:
    SpaceTimeXNeumannIntegrator (const Array<CoefficientFunction*> & coeffs)
    {
      t0 = coeffs[2]->EvaluateConst(); 
      t1 = coeffs[3]->EvaluateConst();
      tau = t1 - t0;

      DomainVariableCoefficientFunction<D> * coefneg = dynamic_cast<DomainVariableCoefficientFunction<D> * > (coeffs[0]);
      if (coefneg != NULL)
      {
        int numreg = coefneg->NumRegions();
        if (numreg == INT_MAX) numreg = 1;
        Array< EvalFunction* > evals;
        evals.SetSize(numreg);
        for (int i = 0; i < numreg; ++i)
        {
          evals[i] = &coefneg->GetEvalFunction(i);
        }
        coef_neg = new DomainVariableCoefficientFunction<D+1>(evals); 
      }
      else
        coef_neg = coeffs[0];

      DomainVariableCoefficientFunction<D> * coefpos = dynamic_cast<DomainVariableCoefficientFunction<D> * > (coeffs[1]);
      if (coefpos != NULL)
      {
        int numreg = coefpos->NumRegions();
        if (numreg == INT_MAX) numreg = 1;
        Array< EvalFunction* > evals(numreg);
        evals.SetSize(numreg);
        for (int i = 0; i < numreg; ++i)
        {
          evals[i] = &coefpos->GetEvalFunction(i);
        }
        coef_pos = new DomainVariableCoefficientFunction<D+1>(evals); 
      }
      else
        coef_pos = coeffs[1];
    }
    virtual ~SpaceTimeXNeumannIntegrator(){ ; };

    virtual string Name () const { return "SpaceTimeXNeumannIntegrator"; }

    virtual int DimElement () const { return D-1; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return true; }

    virtual void SetTimeInterval (const TimeInterval & ti)
    { t0 = ti.first; t1=ti.second; tau = t1-t0; }

    // Calculates the element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> & elvec,
                       LocalHeap & lh) const;
  };

}

#endif

