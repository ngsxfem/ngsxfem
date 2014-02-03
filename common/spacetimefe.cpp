
/// from ngxfem
#include "spacetimefe.hpp"

namespace ngfem
{

  SpaceTimeFiniteElement :: SpaceTimeFiniteElement(const FiniteElement & a_base_space,
                                                 const FiniteElement & a_base_time)
    : base_space(a_base_space),
      base_time(a_base_time)
  { 
    order_space = base_space.Order();
    order_time = base_time.Order();
    ndof_space = base_space.GetNDof();
    ndof_time = base_time.GetNDof();
    ndof = ndof_space * ndof_time;
  };

  SpaceTimeFiniteElement :: ~SpaceTimeFiniteElement() { ; };

  /// the name
  string SpaceTimeFiniteElement :: ClassName(void) const {return "SpaceTimeFiniteElement";};


  template <int D>
  ScalarSpaceTimeFiniteElement<D> :: ScalarSpaceTimeFiniteElement(const ScalarFiniteElement<D> & a_base_space,
                                                                  const DGFiniteElement<1> & a_base_time)
    :  scalar_space(a_base_space),
       scalar_time(a_base_time),
       SpaceTimeFiniteElement(a_base_space, a_base_time)
  { 
  };

  template <int D>
  ScalarSpaceTimeFiniteElement<D> :: ~ScalarSpaceTimeFiniteElement() { ; };

  /// the name
  template <int D>
  string ScalarSpaceTimeFiniteElement<D> :: ClassName(void) const {return "ScalarSpaceTimeFiniteElement";};

  
  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcShapeTime (double time,
                                                   FlatVector<> shape) const
  {
    IntegrationPoint ip(time,0.0,0.0);
    scalar_time.CalcShape(ip,shape);
  };

  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcShapeSpace (const IntegrationPoint & ip, 
                                                    FlatVector<> shape) const
  {
    scalar_space.CalcShape(ip,shape);
  };

  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcShapeSpaceTime (const IntegrationPoint & ip, double time,
                                                              FlatVector<> shape, LocalHeap & lh) const
  {
    FlatVector<> shapex(ndof_space,lh);
    FlatVector<> shapet(ndof_time,lh);
    scalar_space.CalcShape(ip,shapex);
    IntegrationPoint ipt(time,0.0,0.0);
    scalar_time.CalcShape(ipt,shapet);
    for (int j = 0; j < ndof_time; ++j)
      for (int i = 0; i < ndof_space; ++i)
        shape(ndof_space*j+i) = shapet(j) * shapex(i);
  };
  
/*
  /// compute dshape, matrix: ndof x spacedim
  template <int D>
  void SpaceTimeFiniteElement<D> :: CalcDShape (const IntegrationPoint & ip, 
                                       FlatMatrixFixWidth<D> dshape) const
  {
    base.CalcDShape(ip,dshape);
  };


  /// compute dshape, matrix: ndof x spacedim
  template <int D>
  void SpaceTimeFiniteElement<D> :: CalcMappedDShape (const SpecificIntegrationPoint<D,D> & sip, 
                                             FlatMatrixFixWidth<D> dshape) const
  {
    base.CalcMappedDShape(sip,dshape);
  };
    

  /// compute dshape, matrix: ndof x (spacedim spacedim)
  template <int D>
  void SpaceTimeFiniteElement<D> :: CalcDDShape (const IntegrationPoint & ip, 
                                        FlatMatrix<> ddshape) const
  {
    base.CalcDDShape(ip,ddshape);
  };
*/

  template class ScalarSpaceTimeFiniteElement<1>;
  template class ScalarSpaceTimeFiniteElement<2>;
  template class ScalarSpaceTimeFiniteElement<3>;

} // end of namespace
