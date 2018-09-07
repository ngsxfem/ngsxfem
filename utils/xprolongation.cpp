#include "../utils/xprolongation.hpp"

using namespace ngcomp;

namespace ngmg
{

  void P1Prolongation :: Update (const FESpace & afes)
  {
    fes = &afes;
    if (ma->GetNLevels() > nvlevel.Size())
    {
      nvlevel.Append (ma->GetNV());
      //ndoflevel.Append (fes->GetNDof());
    }
    else
      return;

    int nv = ma->GetNV();
    shared_ptr<Array<int>> vert2dof = make_shared<Array<int>>(nv);

    tmp_vecs.Append(make_shared<VVector<double>>(fes->GetNDof()));

    Array<int> dnums(1);
    for (auto i : Range(nv))
    {
      fes->GetDofNrs(NodeId(NT_VERTEX,i),dnums);    
      if (dnums.Size() > 0 && IsRegularDof(dnums[0]))
        (*vert2dof)[i] = dnums[0];
      else
        (*vert2dof)[i] = NO_DOF_NR;
    }
    v2d_on_lvl.Append(vert2dof);
    //cout << "vert2dof : " << *vert2dof << endl;
  }


  void P1Prolongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    if (fes == nullptr)
      throw Exception("call Update before prolongating");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");    
    static Timer t("Prolongate"); RegionTimer r(t);
    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);

    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel-1]->FV<double>();
    fw = fv;
    fv = 0.0;
    for (size_t i = 0; i < nc; i++)
      if (IsRegularDof(vert2dof_fine[i]))
        fv(vert2dof_fine[i]) = fw(vert2dof_coarse[i]);

    for (size_t i = nc; i < nf; i++)
      {
        if (!IsRegularDof(vert2dof_fine[i]))
          continue;

        auto parents = ma->GetParentNodes (i);
        for (auto j : Range(2))
          //if (IsRegularDof(vert2dof_coarse[parents[j]])) //not necessary
          fv(vert2dof_fine[i]) += 0.5 * fw(vert2dof_coarse[parents[j]]);
      }

  }


  void P1Prolongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    if (fes == nullptr)
      throw Exception("call Update before restricting");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");
    static Timer t("Restrict"); RegionTimer r(t);

    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);

    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel]->FV<double>();
    fw = fv;
    fv = 0;
    for (size_t i = 0; i < nc; i++)
    {
      if (!IsRegularDof(vert2dof_coarse[i]))
        continue;
      if (IsRegularDof(vert2dof_fine[i]))
        fv(vert2dof_coarse[i]) = fw(vert2dof_fine[i]);
    }

    for (size_t i = nf; i-- > nc; )
    {
      auto parents = ma->GetParentNodes (i);
      if (!IsRegularDof(vert2dof_fine[i]))
        continue;
      for (auto j : Range(2))
        //if (IsRegularDof(vert2dof_coarse[parents[j]])) // not necessary
        fv(vert2dof_coarse[parents[j]]) += 0.5 * fw(vert2dof_fine[i]);
    } 
  }

}

