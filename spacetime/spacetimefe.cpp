
/// from ngxfem
#include "spacetimefe.hpp"

namespace ngfem
{

  void GR2PFiniteElement::CalcShape (const IntegrationPoint & ip, 
                                    SliceVector<> shape) const
  {
    // std::cout << " GetNDof() = " << GetNDof() << std::endl;
    // std::cout << " ip = " << ip << std::flush << std::endl;
    // std::cout << " shape = " << shape << std::flush << std::endl;
    shape(0) = 1.5*(1-ip(0));
    shape(1) = 1.5*(ip(0)-1.0/3.0);
    // std::cout << " shape = " << shape << std::flush << std::endl;
    // getchar();
  };

  void GR2PFiniteElement::CalcDShape (const IntegrationPoint & ip, 
                                     SliceMatrix<> dshape) const 
  {
    // std::cout << " ip = " << ip << std::flush << std::endl;
    // std::cout << " dshape = " << dshape << std::flush << std::endl;
    dshape(0,0) = -1.5;
    dshape(1,0) = 1.5;
    // std::cout << " dshape = " << dshape << std::flush << std::endl;
    // getchar();
  };

  SpaceTimeFiniteElement :: SpaceTimeFiniteElement(const FiniteElement & a_base_space,
                                                 const FiniteElement & a_base_time)
    : base_space(a_base_space),
      base_time(a_base_time)
  { 
    order_space = base_space.Order();
    order_time = base_time.Order();
    order = order_space;
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
    :  SpaceTimeFiniteElement(a_base_space, a_base_time),
       scalar_space(a_base_space),
       scalar_time(a_base_time)
  { 
  };

  template <int D>
  ScalarSpaceTimeFiniteElement<D> :: ~ScalarSpaceTimeFiniteElement() { ; };

  /// the name
  template <int D>
  string ScalarSpaceTimeFiniteElement<D> :: ClassName(void) const {return "ScalarSpaceTimeFiniteElement";};


  template <int D>
  bool ScalarSpaceTimeFiniteElement<D> :: IsDGFiniteElement() const
  {
    const DGFiniteElement<D> * dg_space = dynamic_cast<const DGFiniteElement<D> *>(& scalar_space);
    return (dg_space != NULL);
  }

  /// diag mass
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: GetDiagMassMatrix(FlatVector<> diagmass, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    const DGFiniteElement<D> * dg_space = dynamic_cast<const DGFiniteElement<D> *>(& scalar_space);

    if (dg_space == NULL)
      throw Exception(" ScalarSpaceTimeFiniteElement<D> :: GetDiagMassMatrix - cast failed");
    
    FlatVector<> massdiagspace(ndof_space,lh);
    dg_space->GetDiagMassMatrix(massdiagspace);

    FlatVector<> massdiagtime(ndof_time,lh);
    scalar_time.GetDiagMassMatrix(massdiagtime);

    for (int m = 0; m < ndof_time; m++)
      for (int n = 0; n < ndof_space; n++)
        diagmass(ndof_space*m+n) = (massdiagtime(m) * massdiagspace(n));
  }


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

  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcDtShapeSpaceTime (const IntegrationPoint & ip, double time,
                                                                FlatVector<> dshape, LocalHeap & lh) const
  {
    FlatVector<> shapex(ndof_space,lh);
    FlatMatrixFixWidth<1> shapet(ndof_time,lh);
    scalar_space.CalcShape(ip,shapex);
    IntegrationPoint ipt(time,0.0,0.0);
    scalar_time.CalcDShape(ipt,shapet);
    for (int j = 0; j < ndof_time; ++j)
      for (int i = 0; i < ndof_space; ++i)
        dshape(ndof_space*j+i) = shapet(j,0) * shapex(i);
  };

  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcDxShapeSpaceTime (const IntegrationPoint & ip, 
                                                                double time,
                                                                FlatMatrixFixWidth<D> dshape, 
                                                                LocalHeap & lh) const
  {
    FlatMatrixFixWidth<D> dshapex(ndof_space,lh); 
    FlatVector<> shapet(ndof_time,lh);
    scalar_space.CalcDShape(ip,dshapex);
    IntegrationPoint ipt(time,0.0,0.0);
    scalar_time.CalcShape(ipt,shapet);
    for (int j = 0; j < ndof_time; ++j)
      for (int i = 0; i < ndof_space; ++i)
        dshape.Row(ndof_space*j+i) = shapet(j) * dshapex.Row(i);
  };

  /// compute shape
  template <int D>
  void ScalarSpaceTimeFiniteElement<D> :: CalcMappedDxShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, 
                                                                      double time,
                                                                      FlatMatrixFixWidth<D> dshape, 
                                                                      LocalHeap & lh) const
  {
    FlatMatrixFixWidth<D> dshapex(ndof_space,lh); 
    FlatVector<> shapet(ndof_time,lh);
    scalar_space.CalcMappedDShape(sip,dshapex);
    IntegrationPoint ipt(time,0.0,0.0);
    scalar_time.CalcShape(ipt,shapet);
    for (int j = 0; j < ndof_time; ++j)
      for (int i = 0; i < ndof_space; ++i)
          dshape.Row(ndof_space*j+i) = shapet(j) * dshapex.Row(i);
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
