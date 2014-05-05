/*********************************************************************/
/* File:   spacetimeintegrators.hpp                                  */
/* Author: C. Lehrenfeld, RWTH                                       */
/* Date:   05. Feb. 2014                                             */
/*********************************************************************/
  
/*  
    Finite Element Integrators for space time elements
*/

#ifndef FILE_SPACETIMEINTEGRATORS_HPP
#define FILE_SPACETIMEINTEGRATORS_HPP
  
#include <fem.hpp>
#include "spacetimefe.hpp" 

namespace ngfem
{

  template <int D>
  class ST_MassIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_mass;
    double told;
    double tnew;
    double dt;
  public:
    ST_MassIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      told = coeffs[0] -> EvaluateConst();
      tnew = coeffs[1] -> EvaluateConst();
      coef_mass = coeffs[2];
      dt = tnew - told;
    }

    virtual ~ST_MassIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const;
  };


  template <int D>
  class ST_TimeDerivativeIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_;
    // double told;
    // double tnew;
    // double dt;
  public:
    ST_TimeDerivativeIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      // told = coeffs[0] -> EvaluateConst();
      // tnew = coeffs[1] -> EvaluateConst();
      coef_ = coeffs[0];
      // dt = tnew - told;
    }

    ST_TimeDerivativeIntegrator (CoefficientFunction* coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_ = coeffs;
    }

    virtual ~ST_TimeDerivativeIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const;

  };

  template <int D>
  class ST_LaplaceIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_;
    double told;
    double tnew;
    double dt;
  public:
    ST_LaplaceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      told = coeffs[0] -> EvaluateConst();
      tnew = coeffs[1] -> EvaluateConst();
      coef_ = coeffs[2];
      dt = tnew - told;
    }

    virtual ~ST_LaplaceIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const;
  };


  template <int D>
  class ST_SourceIntegrator : public LinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_source;
    double told;
    double tnew;
    double dt;
  public:
    ST_SourceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : LinearFormIntegrator()
    { 
      told = coeffs[0] -> EvaluateConst();
      tnew = coeffs[1] -> EvaluateConst();
      dt = tnew - told;
      coef_source = coeffs[2];
    }

    virtual ~ST_SourceIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> & elvec,
                                    LocalHeap & lh) const;
  };



  template <int D>
  class ST_TimeTraceMassIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_mass;
    double time;
  public:
    ST_TimeTraceMassIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_mass  = coeffs[1];
      time = coeffs[0] -> EvaluateConst();
    }

    virtual ~ST_TimeTraceMassIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const;
  };



  template <int D>
  class ST_TimeTraceSourceIntegrator : public LinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_source;
    double time;
  public:
    ST_TimeTraceSourceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : LinearFormIntegrator()
    { 
      coef_source = coeffs[1];
      time = coeffs[0] -> EvaluateConst();
    }

    virtual ~ST_TimeTraceSourceIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> & elvec,
                                    LocalHeap & lh) const;
  };


  enum TIME
  { PAST = 1, FUTURE = 2};

  template <int D, TIME t>  
  class DiffOpTimeTrace : public DiffOp<DiffOpTimeTrace<D,t> >
  {
    
  public:
    enum { DIM = 1 };          // just one copy of the spaces
    enum { DIM_SPACE = D };    // 2D space
    enum { DIM_ELEMENT = D };  // 2D elements (in contrast to 1D boundary elements)
    enum { DIM_DMAT = 1 };     // D-matrix is 1x1
    enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)
    
    static const ScalarSpaceTimeFiniteElement<D> & Cast (const FiniteElement & fel) 
    { return static_cast<const ScalarSpaceTimeFiniteElement<D>&> (fel); }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & bfel, const MIP & sip,
                                MAT & mat, LocalHeap & lh);

    static void GenerateMatrix (const FiniteElement & bfel,
                                const MappedIntegrationPoint<D,D> & mip,
                                FlatMatrix<> & mat, LocalHeap & lh);

    static void GenerateMatrix (const FiniteElement & bfel,
                                const MappedIntegrationPoint<D,D> & mip,
                                FlatMatrixFixHeight<1> & mat, LocalHeap & lh);

    // static void GenerateMatrix (const FiniteElement & bfel,
    //                             const MappedIntegrationPoint<D,D> & mip,
    //                             SliceMatrixColMajor<> & mat, LocalHeap & lh);
  };


  template <int D, TIME t>
  template <typename MIP, typename MAT>
  void DiffOpTimeTrace<D,t>::GenerateMatrix (const FiniteElement & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {
    const ScalarSpaceTimeFiniteElement<D> & cstfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (bfel);

    const int nd = cstfel.GetNDof();
    FlatVector<> shape_st (nd,lh);
    if (t == PAST)
      cstfel.CalcShapeSpaceTime(mip.IP(),0.0,shape_st,lh);
    else if (t==FUTURE)
      cstfel.CalcShapeSpaceTime(mip.IP(),1.0,shape_st,lh);
    mat = Trans(shape_st);
  }

  template <int D, TIME t>
  void DiffOpTimeTrace<D,t>::GenerateMatrix (const FiniteElement & bfel,
                                             const MappedIntegrationPoint<D,D> & mip,
                                             FlatMatrixFixHeight<1> & mat, LocalHeap & lh)
  {
    const ScalarSpaceTimeFiniteElement<D> & cstfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (bfel);

    const int nd = cstfel.GetNDof();
    FlatVector<> shape_st (nd,lh);
    if (t == PAST)
      cstfel.CalcShapeSpaceTime(mip.IP(),0.0,shape_st,lh);
    else if (t==FUTURE)
      cstfel.CalcShapeSpaceTime(mip.IP(),1.0,shape_st,lh);
    mat = Trans(shape_st);
  }

  /*
  template <int D, TIME t>
  void DiffOpTimeTrace<D,t>::GenerateMatrix (const FiniteElement & bfel,
                                             const MappedIntegrationPoint<D,D> & mip,
                                             SliceMatrixColMajor<> & mat, LocalHeap & lh)
  {
    const ScalarSpaceTimeFiniteElement<D> & cstfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (bfel);

    const int nd = cstfel.GetNDof();
    FlatVector<> shape_st (nd,lh);
    if (t == PAST)
      cstfel.CalcShapeSpaceTime(mip.IP(),0.0,shape_st,lh);
    else if (t==FUTURE)
      cstfel.CalcShapeSpaceTime(mip.IP(),1.0,shape_st,lh);
    mat = Trans(shape_st);
  }
  */

  template <int D, TIME t>
  void DiffOpTimeTrace<D,t>::GenerateMatrix (const FiniteElement & bfel,
                                             const MappedIntegrationPoint<D,D> & mip,
                                             FlatMatrix<> & mat, LocalHeap & lh)
  {
    const ScalarSpaceTimeFiniteElement<D> & cstfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (bfel);

    const int nd = cstfel.GetNDof();
    FlatVector<> shape_st (nd,lh);
    if (t == PAST)
      cstfel.CalcShapeSpaceTime(mip.IP(),0.0,shape_st,lh);
    else if (t==FUTURE)
      cstfel.CalcShapeSpaceTime(mip.IP(),1.0,shape_st,lh);
    mat = Trans(shape_st);
  }

  ///
  template <int D, TIME t>
  class NGS_DLL_HEADER SpaceTimeTimeTraceIntegrator 
    : public T_BDBIntegrator<DiffOpTimeTrace<D,t>, DiagDMat<1>, ScalarSpaceTimeFiniteElement<D> >
  {
  public:
    ///
    SpaceTimeTimeTraceIntegrator (Array<CoefficientFunction*> & coeffs);
    ///
    SpaceTimeTimeTraceIntegrator (CoefficientFunction* coeffs);
    ///
    virtual ~SpaceTimeTimeTraceIntegrator (){};
    ///
    virtual string Name () const { return "SpaceTimeTimeTraceIntegratorMass"; }
  };


#ifdef FILE_SPACETIMEINT_CPP
#define SPACETIMEINT_EXTERN
#else
#define SPACETIMEINT_EXTERN extern
  SPACETIMEINT_EXTERN template class DiffOpTimeTrace<2,PAST>;
  SPACETIMEINT_EXTERN template class DiffOpTimeTrace<3,PAST>;
  SPACETIMEINT_EXTERN template class DiffOpTimeTrace<2,FUTURE>;
  SPACETIMEINT_EXTERN template class DiffOpTimeTrace<3,FUTURE>;
  SPACETIMEINT_EXTERN template class T_DifferentialOperator<DiffOpTimeTrace<2,PAST> >;
  SPACETIMEINT_EXTERN template class T_DifferentialOperator<DiffOpTimeTrace<3,PAST> >;
  SPACETIMEINT_EXTERN template class T_DifferentialOperator<DiffOpTimeTrace<2,FUTURE> >;
  SPACETIMEINT_EXTERN template class T_DifferentialOperator<DiffOpTimeTrace<3,FUTURE> >;
#endif


}

#endif
