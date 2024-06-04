#pragma once

/// from ngsxfem
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
    Array<shared_ptr<VVector<double> > > tmp_vecs;
    const FESpace* fes;
    Array<shared_ptr<Array<int>>> v2d_on_lvl;
  public:
    P1Prolongation(shared_ptr<MeshAccess> ama)
      : ma(ama), fes(nullptr) 
      { 
        tmp_vecs.SetSize(0); 
        nvlevel.SetSize(0);
        v2d_on_lvl.SetSize(0);
      }
    
    virtual ~P1Prolongation() { ; }

    virtual void Update (const FESpace & fes) override;

    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("P1Prolongation::CreateProlongationMatrix not implemented!");
    }
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };

  class P2Prolongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    Array<size_t> nVertLevel;
    Array<size_t> nEdgeLevel;
    //Array<size_t> ndoflevel;
    Array<shared_ptr<BaseVector>> tmp_vecs;
    const FESpace* fes;
    // Array<shared_ptr<Array<int>>> v2d_on_lvl;
  public:
    P2Prolongation(shared_ptr<MeshAccess> ama)
      : ma(ama), fes(nullptr) 
      { 
        // tmp_vecs.SetSize(0); 
        nVertLevel.SetSize(0);
        nEdgeLevel.SetSize(0);
      }
    
    virtual ~P2Prolongation() { cout << "p2prolongation dying"; }

    virtual void Update (const FESpace & fes) override;

    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("P2Prolongation::CreateProlongationMatrix not implemented!");
    }
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };


  class P2CutProlongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    Array<size_t> nVertLevel;
    Array<size_t> nEdgeLevel;
    //Array<size_t> ndoflevel;
    Array<shared_ptr<BaseVector>> tmp_vecs;
    const FESpace* fes;
    Array<shared_ptr<Array<int>>> v2d_on_lvl;
  public:
    P2CutProlongation(shared_ptr<MeshAccess> ama)
      : ma(ama), fes(nullptr) 
      { 
        // tmp_vecs.SetSize(0); 
        nVertLevel.SetSize(0);
        nEdgeLevel.SetSize(0);
      }
    
    virtual ~P2CutProlongation() { cout << "p2prolongation dying"; }

    virtual void Update (const FESpace & fes) override;

    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    {
      throw Exception("P2Prolongation::CreateProlongationMatrix not implemented!");
    }
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override;
    virtual void RestrictInline (int finelevel, BaseVector & v) const override;
  };

}
