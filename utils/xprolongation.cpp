#include "../utils/xprolongation.hpp"

using namespace ngcomp;

namespace ngmg
{

  void P1Prolongation :: Update (const FESpace & fes)
  {
    if (ma->GetNLevels() > nvlevel.Size())
      nvlevel.Append (ma->GetNV());
    cout << " nooooooo " << endl;
  }


  void P1Prolongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    static Timer t("Prolongate"); RegionTimer r(t);
    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];

    if (v.EntrySize() == 1)
      {
        FlatVector<> fv = v.FV<double>();
        fv.Range (nf, fv.Size()) = 0;
        for (size_t i = nc; i < nf; i++)
          {
            auto parents = ma->GetParentNodes (i);
            fv(i) = 0.5 * (fv(parents[0]) + fv(parents[1]));
          }
      }
    else
      {
        FlatSysVector<> sv = v.SV<double>();
        sv.Range (nf, sv.Size()) = 0;
        for (size_t i = nc; i < nf; i++)
          {
            auto parents = ma->GetParentNodes (i);
            sv(i) = 0.5 * (sv(parents[0]) + sv(parents[1]));
          }
      }
  }


  void P1Prolongation :: RestrictInline (int finelevel, BaseVector & v) const
    {
      static Timer t("Restrict"); RegionTimer r(t);

      size_t nc = nvlevel[finelevel-1];
      size_t nf = nvlevel[finelevel];

      FlatSysVector<> fv = v.SV<double>();

      for (size_t i = nf; i-- > nc; )
      {
        auto parents = ma->GetParentNodes (i);
        fv(parents[0]) += 0.5 * fv(i);
        fv(parents[1]) += 0.5 * fv(i);
      }

      for (size_t i = nf; i < fv.Size(); i++)
        fv(i) = 0;
    }


}

