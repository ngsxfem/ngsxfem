#pragma once

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>


using namespace ngsolve;
// using namespace cutinfo;

namespace ngcomp
{


  class ElementAggregation
  {
    
  protected:
    shared_ptr<MeshAccess> ma;
    
  public:
    ElementAggregation (shared_ptr<MeshAccess> ama);
    
    shared_ptr<MeshAccess> GetMesh () const { return ma; }
  
    
  
  };



}