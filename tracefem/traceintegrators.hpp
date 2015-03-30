#ifndef FILE_TRACEINTEGRATORS_HPP
#define FILE_TRACEINTEGRATORS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"
#include "../utils/stcoeff.hpp"

namespace ngfem
{
  template <int D>
  class TraceMassIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    TraceMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0]){;}
    virtual ~TraceMassIntegrator(){ ; };
    virtual string Name () const { return "TraceMassIntegrator"; }
    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    virtual bool BoundaryForm () const { return false; }
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                         const ElementTransformation & eltrans,
                         FlatMatrix<double> elmat,
                         LocalHeap & lh) const;
  };

  template <int D>
  class TraceSourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    TraceSourceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef(coeffs[0]){;}
    virtual ~TraceSourceIntegrator(){ ; };
    virtual string Name () const { return "TraceSourceIntegrator"; }
    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    virtual bool BoundaryForm () const { return false; }
    virtual void
    CalcElementVector (const FiniteElement & fel,
                         const ElementTransformation & eltrans,
                         FlatVector<double> elvec,
                         LocalHeap & lh) const;
  };

  template <int D>
  class DiffOpEvalExtTrace : public DiffOp<DiffOpEvalExtTrace<D> >
  {

  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // 2D space
    enum { DIM_ELEMENT = D };  // 2D elements (in contrast to 1D boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix is 2x2
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

  ///
  template <int D >
  class ExtTraceIntegrator
    : public T_BDBIntegrator<DiffOpEvalExtTrace<D>, DiagDMat<1>, CompoundFiniteElement >
  {
  public:
    ///
    ExtTraceIntegrator (shared_ptr<CoefficientFunction> coeff);
    ///
    ExtTraceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs);
    ///
    virtual ~ExtTraceIntegrator ();
    ///
    virtual string Name () const { return "ExtTraceIntegrator"; }
  };

}

#endif
