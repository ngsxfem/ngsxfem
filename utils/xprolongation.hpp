#pragma once

/// from ngxfem
#include <multigrid.hpp>
//#include <comp.hpp>
// using namespace ngsolve;
// using namespace cutinfo;
// using namespace ngcomp;

namespace ngmg
{

  class P1Prolongation : public LinearProlongation
  {
    shared_ptr<MeshAccess> ma;
    Array<size_t> nvlevel;
    shared_ptr<FESpace> fes;
  public:
    P1Prolongation(shared_ptr<MeshAccess> ama)
      : LinearProlongation(ama) { ; }
    
    virtual ~P1Prolongation() { ; }

    virtual void Update (const FESpace & fes) override;

    virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("P1Prolongation::CreateProlongationMatrix not implemented!");
    }
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };

}
