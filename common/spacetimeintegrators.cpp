/*********************************************************************/
/* File:   spacetimeintegrators.cpp                                  */
/* Author: C. Lehrenfeld, RWTH                                       */
/* Date:   03. Feb. 2014                                             */
/*********************************************************************/
  
/*  
    Finite Element Integrators for space time elements
*/
  
#include <fem.hpp>
#include "spacetimefe.hpp" 

namespace ngfem
{


  template <int D>
  class ST_MassIntegrator : public BilinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_mass;
  public:
    ST_MassIntegrator (Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator()
    { 
      coef_mass  = coeffs[0];
      // alpha = coeffs[1] -> EvaluateConst();
    }

    virtual ~ST_MassIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> & elmat,
                                    LocalHeap & lh) const
    {

      const ScalarSpaceTimeFiniteElement<D> & stfel = 
        dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
      ELEMENT_TYPE eltype = stfel.ElementType();
      
      const int order_space = stfel.OrderSpace();
      const int order_time = stfel.OrderTime();

      const int nd_space = stfel.GetNDofSpace();
      const int nd_time = stfel.GetNDofTime();

      const int nd = stfel.GetNDof();

      elmat = 0.0;

      // scaling with time and space is missing
      // time interval is missing
      // how to set time?
      if (stfel.IsDGFiniteElement()) //assume constant coefficients and DGFiniteElements
      {
        double m = coef_mass -> EvaluateConst();
        FlatVector<> diagmass(nd,lh);
        stfel.GetDiagMassMatrix(diagmass,lh);
        for (int i = 0; i < nd; ++i)
          elmat(i,i) = diagmass(i);

      }
      else
      {
        const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
        const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

        FlatVector<> shape_st (nd,lh);
        for (int l = 0; l < ir_time.GetNIP(); l++)
        {
          double current = ir_time[l](0);
          // *time = current;
          for (int k = 0; k < ir_space.GetNIP(); k++)
          {
            stfel.CalcShapeSpaceTime(ir_space[k],current,shape_st,lh);
            
            const double fac = ir_time[l].Weight() * ir_space[k].Weight();
            elmat += fac * shape_st * Trans(shape_st);
          }
        }
        // throw Exception(" mass matrix - not implemented");
      }
      // cout << " elmat = " << elmat << endl;
    }
  };


  template <int D>
  class ST_SourceIntegrator : public LinearFormIntegrator
  {
  protected:
    CoefficientFunction *coef_source;
  public:
    ST_SourceIntegrator (Array<CoefficientFunction*> & coeffs) 
      : LinearFormIntegrator()
    { 
      coef_source  = coeffs[0];
      // alpha = coeffs[1] -> EvaluateConst();
    }

    virtual ~ST_SourceIntegrator () { ; }

    virtual bool BoundaryForm () const 
    { return 0; }

    virtual void CalcElementVector (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatVector<double> & elvec,
                                    LocalHeap & lh) const
    {

      const ScalarSpaceTimeFiniteElement<D> & stfel = 
        dynamic_cast<const ScalarSpaceTimeFiniteElement<D> &> (fel);
  
      ELEMENT_TYPE eltype = stfel.ElementType();
      
      const int order_space = stfel.OrderSpace();
      const int order_time = stfel.OrderTime();

      const int nd_space = stfel.GetNDofSpace();
      const int nd_time = stfel.GetNDofTime();

      const int nd = stfel.GetNDof();

      elvec = 0.0;

      // scaling with time and space is missing
      // time interval is missing
      // how to set time?
      // if (stfel.IsDGFiniteElement()) //assume constant coefficients and DGFiniteElements
      // {
      //   double m = coef_source -> EvaluateConst();
      //   FlatVector<> diagsource(nd,lh);
      //   stfel.GetDiagSourceMatrix(diagsource,lh);
      //   for (int i = 0; i < nd; ++i)
      //     elmat(i,i) = diagsource(i);

      // }
      // else
      {
        const IntegrationRule & ir_space = SelectIntegrationRule (eltype, 2*order_space);
        const IntegrationRule & ir_time = SelectIntegrationRule (ET_SEGM, 2*order_time);

        FlatVector<> shape_st (nd,lh);
        for (int l = 0; l < ir_time.GetNIP(); l++)
        {
          double current = ir_time[l](0);
          // *time = current;
          for (int k = 0; k < ir_space.GetNIP(); k++)
          {
            stfel.CalcShapeSpaceTime(ir_space[k],current,shape_st,lh);
            
            MappedIntegrationPoint<D,D> mip(ir_space[k],eltrans);
            double m = coef_source -> Evaluate(mip);
            const double fac = ir_time[l].Weight() * ir_space[k].Weight();
            elvec += fac * shape_st;
          }
        }
        // throw Exception(" source matrix - not implemented");
      }
      // cout << " elmat = " << elmat << endl;
    }
  };


  static RegisterBilinearFormIntegrator<ST_MassIntegrator<2> > initstmass21 ("STmass", 2, 1);
  static RegisterBilinearFormIntegrator<ST_MassIntegrator<3> > initstmass31 ("STmass", 3, 1);

  static RegisterLinearFormIntegrator<ST_SourceIntegrator<2> > initstsource21 ("STsource", 2, 1);
  static RegisterLinearFormIntegrator<ST_SourceIntegrator<3> > initstsource31 ("STsource", 3, 1);

}

