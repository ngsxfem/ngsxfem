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

  class GR2PFiniteElement : public DGFiniteElement<1>
  {

  protected:
    int vnums[1<<1];  

    using ScalarFiniteElement<1>::ndof;
    using ScalarFiniteElement<1>::order;

  public:

    virtual string ClassName(void) const {return "GaussRadau2Point orthogonal Finite Element";}

    virtual void SetOrder (INT<1> p){ order = 1;};
    virtual void ComputeNDof(){};


    virtual void PrecomputeTrace (){};
    virtual void PrecomputeGrad (){};

    virtual void GetDiagMassMatrix (FlatVector<> mass) const { mass(0) = 0.75; mass(1) = 0.25; };

    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceVector<> shape) const;
    
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const;

    virtual ELEMENT_TYPE ElementType() const{ return ET_SEGM; }

    GR2PFiniteElement()  
    {
      ndof=2;
      order=1;
    }

  };
  



 // GR1FiniteElement:


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

} // end of namespace


#endif
