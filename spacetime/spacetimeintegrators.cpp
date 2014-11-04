/*********************************************************************/
/* File:   spacetimeintegrators.cpp                                  */
/* Author: C. Lehrenfeld, RWTH                                       */
/* Date:   03. Feb. 2014                                             */
/*********************************************************************/
  
/*  
    Finite Element Integrators for space time elements
*/

#define FILE_SPACETIMEINT_CPP  
#include "spacetimeintegrators.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{
  template <int D>
  void ST_MassIntegrator<D> :: CalcElementMatrix (const FiniteElement & fel,
                                                  const ElementTransformation & eltrans, 
                                                  FlatMatrix<double> elmat,
                                                  LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elmat = 0.0;

    if (false && stfel.IsDGFiniteElement()) //assume constant coefficients and DGFiniteElements
    {
      double m = coef_mass -> EvaluateConst();
      FlatVector<> diagmass(nd,lh);
      stfel.GetDiagMassMatrix(diagmass,lh);
      for (int i = 0; i < nd; ++i)
        elmat(i,i) = m * dt * diagmass(i);

    }
    else
    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
      const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

      FlatVector<> shape_st (nd,lh);
      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        // double current = told + ir_time[l](0) * dt;
        // *time = current;

        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          stfel.CalcShapeSpaceTime(ir_space[k],ir_time[l](0),shape_st,lh);
          MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);
          double m = coef_mass -> Evaluate(mip);
          const double fac = m * dt * ir_time[l].Weight() * mip.GetWeight();
          elmat += fac * shape_st * Trans(shape_st);
        }
      }
    }
  }

  template <int D>
  void ST_TimeDerivativeIntegrator<D> :: CalcElementMatrix (const FiniteElement & fel,
                                                            const ElementTransformation & eltrans, 
                                                            FlatMatrix<double> elmat,
                                                            LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elmat = 0.0;

    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
      const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

      FlatVector<> shape_st (nd,lh);
      FlatVector<> dtshape_st (nd,lh);
      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        // double current = told + ir_time[l](0) * dt;
        // *time = current;

        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          stfel.CalcDtShapeSpaceTime(ir_space[k],ir_time[l](0),dtshape_st,lh);
          stfel.CalcShapeSpaceTime(ir_space[k],ir_time[l](0),shape_st,lh);
          MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);
          double m = coef_ -> Evaluate(mip);
          const double fac = m * ir_time[l].Weight() * mip.GetWeight();
          elmat += fac * shape_st * Trans(dtshape_st);
        }
      }
    }
  }

  template <int D>
  void ST_LaplaceIntegrator<D> :: CalcElementMatrix (const FiniteElement & fel,
                                                     const ElementTransformation & eltrans, 
                                                     FlatMatrix<double> elmat,
                                                     LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elmat = 0.0;

    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
      const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

      FlatMatrixFixWidth<D> dshape_st (nd,lh);
      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        // double current = told + ir_time[l](0) * dt;
        // *time = current;

        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);
          stfel.CalcMappedDxShapeSpaceTime(mip,ir_time[l](0),dshape_st,lh);
          double m = coef_ -> Evaluate(mip);
          const double fac = m * dt * ir_time[l].Weight() * mip.GetWeight();
          elmat += fac * dshape_st * Trans(dshape_st);
        }
      }
    }
  }
    
  template <int D>
  void ST_SourceIntegrator<D> :: CalcElementVector (const FiniteElement & fel,
                                                    const ElementTransformation & eltrans, 
                                                    FlatVector<double> & elvec,
                                                    LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elvec = 0.0;

    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
      const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

      FlatVector<> shape_st (nd,lh);
      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        // double current = told + ir_time[l](0) * dt;
        // *time = current;
        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          stfel.CalcShapeSpaceTime(ir_space[k],ir_time[l](0),shape_st,lh);
            
          MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);
          double m = coef_source -> Evaluate(mip);
          const double fac = m * dt * ir_time[l].Weight() * mip.GetWeight();
          elvec += fac * shape_st;
        }
      }
    }
  }

  template <int D>
  void ST_TimeTraceMassIntegrator<D> :: CalcElementMatrix (const FiniteElement & fel,
                                                           const ElementTransformation & eltrans, 
                                                           FlatMatrix<double> elmat,
                                                           LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    // const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elmat = 0.0;

    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);

      FlatVector<> shape_st (nd,lh);

      for (int k = 0; k < ir_space.GetNIP(); k++)
      {
        stfel.CalcShapeSpaceTime(ir_space[k],time,shape_st,lh);
        MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);  
        double m = coef_mass -> Evaluate(mip);
        const double fac = m * mip.GetWeight();
        elmat += fac * shape_st * Trans(shape_st);
      }
    }
  }

  template <int D>
  void ST_TimeTraceSourceIntegrator<D> :: CalcElementVector (const FiniteElement & fel,
                                                             const ElementTransformation & eltrans, 
                                                             FlatVector<double> & elvec,
                                                             LocalHeap & lh) const
  {

    const ScalarSpaceTimeFiniteElement<D> & stfel = 
      dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
    ELEMENT_TYPE eltype = stfel.ElementType();
      
    const int order_space = stfel.OrderSpace();
    // const int order_time = stfel.OrderTime();

    // const int nd_space = stfel.GetNDofSpace();
    // const int nd_time = stfel.GetNDofTime();

    const int nd = stfel.GetNDof();

    elvec = 0.0;

    {
      const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);

      FlatVector<> shape_st (nd,lh);

      for (int k = 0; k < ir_space.GetNIP(); k++)
      {
        stfel.CalcShapeSpaceTime(ir_space[k],time,shape_st,lh);
        MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);  
        double m = coef_source -> Evaluate(mip);
        const double fac = m * mip.GetWeight();
        elvec += fac * shape_st;
      }
    }
  }

  static RegisterBilinearFormIntegrator<ST_MassIntegrator<2> > initstmass21 ("STmass", 2, 3);
  static RegisterBilinearFormIntegrator<ST_MassIntegrator<3> > initstmass31 ("STmass", 3, 3);

  static RegisterBilinearFormIntegrator<ST_TimeDerivativeIntegrator<2> > initsttimeder21 ("STtimeder", 2, 1);
  static RegisterBilinearFormIntegrator<ST_TimeDerivativeIntegrator<3> > initsttimeder31 ("STtimeder", 3, 1);

  static RegisterBilinearFormIntegrator<ST_LaplaceIntegrator<2> > initstlap21 ("STlaplace", 2, 3);
  static RegisterBilinearFormIntegrator<ST_LaplaceIntegrator<3> > initstlap31 ("STlaplace", 3, 3);

  static RegisterBilinearFormIntegrator<ST_TimeTraceMassIntegrator<2> > initsttrmass21 ("STtracemass", 2, 2);
  static RegisterBilinearFormIntegrator<ST_TimeTraceMassIntegrator<3> > initsttrmass31 ("STtracemass", 3, 2);

  static RegisterLinearFormIntegrator<ST_SourceIntegrator<2> > initstsource21 ("STsource", 2, 3);
  static RegisterLinearFormIntegrator<ST_SourceIntegrator<3> > initstsource31 ("STsource", 3, 3);

  static RegisterLinearFormIntegrator<ST_TimeTraceSourceIntegrator<2> > initsttrsrc21 ("STtracesource", 2, 2);
  static RegisterLinearFormIntegrator<ST_TimeTraceSourceIntegrator<3> > initsttrsrc31 ("STtracesource", 3, 2);

  template <int D, TIME t>  
  SpaceTimeTimeTraceIntegrator<D,t> :: SpaceTimeTimeTraceIntegrator (const Array<shared_ptr<CoefficientFunction> > & coeffs)
    : T_BDBIntegrator<DiffOpTimeTrace<D,t>, DiagDMat<1>, ScalarSpaceTimeFiniteElement<D> > (coeffs)
  { ; }

  template <int D, TIME t>  
  SpaceTimeTimeTraceIntegrator<D,t> :: SpaceTimeTimeTraceIntegrator (shared_ptr<CoefficientFunction> coeffs)
    : T_BDBIntegrator<DiffOpTimeTrace<D,t>, DiagDMat<1>, ScalarSpaceTimeFiniteElement<D> > (coeffs)
  { ; }

  template class SpaceTimeTimeTraceIntegrator<2, PAST>;
  template class SpaceTimeTimeTraceIntegrator<3, PAST>;
  template class SpaceTimeTimeTraceIntegrator<2, FUTURE>;
  template class SpaceTimeTimeTraceIntegrator<3, FUTURE>;

  static RegisterBilinearFormIntegrator<SpaceTimeTimeTraceIntegrator<2, PAST> > initsttrmasspast2 ("STtracepast", 2, 1);
  static RegisterBilinearFormIntegrator<SpaceTimeTimeTraceIntegrator<3, PAST> > initsttrmasspast3 ("STtracepast", 3, 1);
  static RegisterBilinearFormIntegrator<SpaceTimeTimeTraceIntegrator<2, FUTURE> > initsttrmassfut2 ("STtracefuture", 2, 1);
  static RegisterBilinearFormIntegrator<SpaceTimeTimeTraceIntegrator<3, FUTURE> > initsttrmassfut3 ("STtracefuture", 3, 1);



}

