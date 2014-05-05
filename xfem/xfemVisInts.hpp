#ifndef FILE_XFEMVISINTS_HPP
#define FILE_XFEMVISINTS_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

// #include "../cuttriang/geom.hpp"
#include "xfiniteelement.hpp"
#include "../spacetime/spacetimeintegrators.hpp"

namespace ngfem
{

  template <int D>
  class DiffOpEvalX : public DiffOp<DiffOpEvalX<D> >
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
  class XVisIntegrator 
    : public T_BDBIntegrator<DiffOpEvalX<D>, DiagDMat<1>, CompoundFiniteElement >
  {
  public:
    ///
    XVisIntegrator (CoefficientFunction * coeff);
    ///
    XVisIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~XVisIntegrator ();
    ///
    virtual string Name () const { return "XVisIntegrator"; }
  };


  template <int D, TIME t>
  class DiffOpEvalSTX : public DiffOp<DiffOpEvalSTX<D,t> >
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
  template <int D, TIME t>
  class STXVisIntegrator 
    : public T_BDBIntegrator<DiffOpEvalSTX<D,t>, DiagDMat<1>, CompoundFiniteElement >
  {
  public:
    ///
    STXVisIntegrator (CoefficientFunction * coeff);
    ///
    STXVisIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~STXVisIntegrator ();
    ///
    virtual string Name () const { return "STXVisIntegrator"; }
  };




  template <int D, TIME t>
  class DiffOpEvalSTNegPos : public DiffOp<DiffOpEvalSTNegPos<D,t> >
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
  template <int D, TIME t>
  class STNegPosVisIntegrator 
    : public T_BDBIntegrator<DiffOpEvalSTNegPos<D,t>, DiagDMat<1>, CompoundFiniteElement >
  {
  public:
    ///
    STNegPosVisIntegrator (CoefficientFunction * coeff);
    ///
    STNegPosVisIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~STNegPosVisIntegrator ();
    ///
    virtual string Name () const { return "STNegPosVisIntegrator"; }
  };





  template <int D>
  class DiffOpEvalSigned : public DiffOp<DiffOpEvalSigned<D> >
  {
    
  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // 2D space
    enum { DIM_ELEMENT = D };  // 2D elements (in contrast to 1D boundary elements)
    enum { DIM_DMAT = 2 };     // D-matrix is 2x2
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };

  template <int D, TIME t>
  class DiffOpEvalSpaceTimeSigned : public DiffOp<DiffOpEvalSpaceTimeSigned<D,t> >
  {
    
  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // 2D space
    enum { DIM_ELEMENT = D };  // 2D elements (in contrast to 1D boundary elements)
    enum { DIM_DMAT = 2 };     // D-matrix is 2x2
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)
    
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);
  };


  /// diagonal tensor, one on first diag entry if lvlset-value is positive, 
  /// one on second diag entry otherwise
  class XHeavisideDMat : public DMatOp<XHeavisideDMat,2>
  {
    CoefficientFunction * lvlset;
  public:
    // typedef SCAL TSCAL;
    enum { DIM_DMAT = 2 };
    XHeavisideDMat (CoefficientFunction * acoef) : lvlset(acoef) { ; }
    
    XHeavisideDMat (Array<CoefficientFunction*> & acoefs) : lvlset(acoefs[0]) { ; }
    
    template <typename FEL, typename MIP, typename MAT>
    void GenerateMatrix (const FEL & fel, const MIP & mip,
                         MAT & mat, LocalHeap & lh) const;
  };

  ///
  template <int D >
  class SignedXMassIntegrator 
    : public T_BDBIntegrator<DiffOpEvalSigned<D>, XHeavisideDMat, CompoundFiniteElement >
  {
  public:
    ///
    SignedXMassIntegrator (CoefficientFunction * coeff);
    ///
    SignedXMassIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~SignedXMassIntegrator ();
    ///
    virtual string Name () const { return "SignedXMassIntegrator"; }
  };

  ///
  template <int D, TIME t >
  class SignedSpaceTimeXMassIntegrator 
    : public T_BDBIntegrator<DiffOpEvalSpaceTimeSigned<D,t>, XHeavisideDMat, CompoundFiniteElement >
  {
  public:
    ///
    SignedSpaceTimeXMassIntegrator (CoefficientFunction * coeff);
    ///
    SignedSpaceTimeXMassIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~SignedSpaceTimeXMassIntegrator ();
    ///
    virtual string Name () const { return "SignedSpaceTimeXMassIntegrator"; }
  };



  template <int D>
  class DiffOpEvalFict : public DiffOp<DiffOpEvalFict<D> >
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
  class FictVisIntegrator 
    : public T_BDBIntegrator<DiffOpEvalFict<D>, DiagDMat<1>, FiniteElement >
  {
  public:
    ///
    FictVisIntegrator (CoefficientFunction * coeff);
    ///
    FictVisIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    virtual ~FictVisIntegrator ();
    ///
    virtual string Name () const { return "FictVisIntegrator"; }
  };




}

#endif

