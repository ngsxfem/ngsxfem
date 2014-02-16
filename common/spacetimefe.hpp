#ifndef FILE_SPACETIMEFE_HPP
#define FILE_SPACETIMEFE_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

using namespace ngsolve;

namespace ngfem
{

typedef std::pair<double,double> TimeInterval;

/**
     
 */
  class SpaceTimeFiniteElement : public FiniteElement
  {
  protected:
    const FiniteElement & base_space;
    const FiniteElement & base_time;
    int order_space;
    int order_time;
    int ndof_space;
    int ndof_time;
  public:
    SpaceTimeFiniteElement(const FiniteElement & base_space,
                           const FiniteElement & base_time);
    virtual ~SpaceTimeFiniteElement();

    /// the name
    virtual string ClassName(void) const;

    int GetNDofSpace() const { return ndof_space;}
    int GetNDofTime() const { return ndof_time;}
    int OrderSpace() const { return order_space;}
    int OrderTime() const { return order_time;}
    virtual int GetNDof() const { return ndof;}

  };

  template <int D>  
  class ScalarSpaceTimeFiniteElement : public SpaceTimeFiniteElement
  {
  protected:
    const ScalarFiniteElement<D> & scalar_space;
    const DGFiniteElement<1> & scalar_time;
  public:
    ScalarSpaceTimeFiniteElement(const ScalarFiniteElement<D> & base_space,
                                 const DGFiniteElement<1> & base_time);
    virtual ~ScalarSpaceTimeFiniteElement();

    virtual int GetNDof() const { return ndof;}

    /// the name
    virtual string ClassName(void) const;

    virtual bool IsDGFiniteElement() const;
    /// mass diag
    virtual void GetDiagMassMatrix(FlatVector<> diagmass,
                                   LocalHeap & lh) const;

    /// compute shape
    virtual void CalcShapeTime (double time,
                                FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpace (const IntegrationPoint & ip,
                                 FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpaceTime (const IntegrationPoint & ip, double time,
                                     FlatVector<> shape, LocalHeap & lh) const;

    /// compute dshape, vector: ndof 
    virtual void CalcDtShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatVector<> dshape, LocalHeap & lh) const;

    virtual ELEMENT_TYPE ElementType() const{ return base_space.ElementType(); }

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDxShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDxShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time,
                                             FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;
  };
    
    // /// compute dshape, vector: ndof 
    // virtual void CalcMappedDtShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time
    //                                          FlatVector<> dshape) const;

    // /// compute dshape, matrix: ndof x (spacedim spacedim)
    // // virtual void CalcDDShape (const IntegrationPoint & ip, 
    // //                           FlatMatrix<> ddshape) const;

  class ScalarFieldEvaluator
  {
  public:
    virtual double operator()(const Vec<1>& point) const
    {
      throw Exception(" nonono 1");
    }

    virtual double operator()(const Vec<2>& point) const
    {
      throw Exception(" nonono 2");
    }

    virtual double operator()(const Vec<3>& point) const
    {
      throw Exception(" nonono 3");
    }

    virtual double operator()(const Vec<4>& point) const
    {
      throw Exception(" nonono 4");
    }

    static ScalarFieldEvaluator* Create(int dim, const FiniteElement & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const EvalFunction & evalf, const ElementTransformation& eltrans, LocalHeap & a_lh);
    static ScalarFieldEvaluator* Create(int dim, const EvalFunction & evalf, const ElementTransformation& eltrans, const TimeInterval & ti, LocalHeap & a_lh);

  };

  template <int D>
  class ScalarFEEvaluator : public ScalarFieldEvaluator
  {
  protected:
    const ScalarSpaceTimeFiniteElement<D> * st_fe;
    const ScalarFiniteElement<D> * s_fe;
    FlatVector<> linvec;
    mutable IntegrationPoint ip;
    LocalHeap & lh;
    mutable double fixedtime = 0;
    mutable bool timefixed = false;
  public:
    ScalarFEEvaluator(const FiniteElement & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh)
      : linvec(a_linvec), 
        lh(a_lh)
    {
      st_fe = dynamic_cast< const ScalarSpaceTimeFiniteElement<D> * >(&a_fe);
      s_fe = dynamic_cast< const ScalarFiniteElement<D> * >(&a_fe);

      if (st_fe == NULL && s_fe == NULL)
      {
        cout << " D = " << D << endl;
        throw Exception("ScalarFEEvaluator - constructor: cast failed...");
      }
    }
    
    void FixTime(double a_fixedtime ) const 
    {
      timefixed = true;
      fixedtime = a_fixedtime;
    }

    void UnFixTime() const
    {
      timefixed = false;
    }
      
    virtual double operator()(const Vec<D>& point) const;
    virtual double operator()(const Vec<D+1>& point) const;
  };


  template <int D> // D : resulting space dimension..
  class EvalFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    EvalFunction eval;
    const ElementTransformation & eltrans;
  public:
    EvalFunctionEvaluator( const string & str, const ElementTransformation & a_eltrans) 
      : eval(str), eltrans(a_eltrans) { ; }

    EvalFunctionEvaluator( const EvalFunction & str, const ElementTransformation & a_eltrans) 
      : eval(str), eltrans(a_eltrans) { ; }

    virtual double operator()(const Vec<1>& point) const
    {
      IntegrationPoint ip(point(0));
      MappedIntegrationPoint<1,D> mip(ip, eltrans);
      return eval.Eval(& mip.GetPoint()(0));
    }

    virtual double operator()(const Vec<2>& point) const
    {
      IntegrationPoint ip(point(0),point(1));
      MappedIntegrationPoint<2,D> mip(ip, eltrans);
      return eval.Eval(& mip.GetPoint()(0));
    }

    virtual double operator()(const Vec<3>& point) const
    {
      IntegrationPoint ip(point(0),point(1),point(2));
      MappedIntegrationPoint<3,D> mip(ip, eltrans);
      return eval.Eval(& mip.GetPoint()(0));
    }

  };

  template <int D> // D : resulting space dimension..
  class SpaceTimeEvalFunctionEvaluator : public ScalarFieldEvaluator
  {
  protected:
    EvalFunction eval;
    const ElementTransformation & eltrans;
    TimeInterval ti;
  public:
    SpaceTimeEvalFunctionEvaluator( const string & str, 
                                    const ElementTransformation & a_eltrans, 
                                    const TimeInterval & a_ti) 
      : eval(str), eltrans(a_eltrans), ti(a_ti) { ; }

    SpaceTimeEvalFunctionEvaluator( const EvalFunction & str, 
                                    const ElementTransformation & a_eltrans, 
                                    const TimeInterval & a_ti) 
      : eval(str), eltrans(a_eltrans), ti(a_ti) { ; }

    virtual double operator()(const Vec<1>& point) const
    {
      Vec<3> p(0.0); 
      p(0) = ( (1.0-point(0)) * ti.first + point(0) * ti.second );
      return eval.Eval(& p(0));
    }

    virtual double operator()(const Vec<2>& point) const
    {
      IntegrationPoint ip(point(0));
      MappedIntegrationPoint<1,D> mip(ip, eltrans);
      Vec<3> p = mip.GetPoint();
      p(D) = ( (1.0-point(1)) * ti.first + point(1) * ti.second );
      return eval.Eval(& p(0));
    }

    virtual double operator()(const Vec<3>& point) const
    {
      IntegrationPoint ip(point(0),point(1));
      MappedIntegrationPoint<2,D> mip(ip, eltrans);
      Vec<3> p = mip.GetPoint();
      p(D) = ( (1.0-point(2)) * ti.first + point(2) * ti.second );
      return eval.Eval(& p(0));
    }

  };



} // end of namespace


#endif
