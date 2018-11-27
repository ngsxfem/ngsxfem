#pragma once

/// from ngxfem
#include <multigrid.hpp>
//#include <comp.hpp>
// using namespace ngsolve;
// using namespace cutinfo;
// using namespace ngcomp;

namespace ngmg
{

  class P1Prolongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    Array<size_t> nvlevel;
    //Array<size_t> ndoflevel;
    Array<shared_ptr<BaseVector>> tmp_vecs;
    const FESpace* fes;
    Array<shared_ptr<Array<int>>> v2d_on_lvl;
  public:
    P1Prolongation(shared_ptr<MeshAccess> ama)
      : ma(ama), fes(nullptr) { ; }
    
    virtual ~P1Prolongation() { cout << "p1prolongation dying"; }

    virtual void Update (const FESpace & fes) override;

    virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("P1Prolongation::CreateProlongationMatrix not implemented!");
    }
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };

}
