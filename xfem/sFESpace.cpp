#include "../xfem/sFESpace.hpp"
#include "../xfem/xfemdiffops.hpp"

using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{
  void SFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    // cout << "GetDofNrs" << endl;
    // cout << firstdof_of_el << endl;
    if ((ei.VB() == VOL ) && (activeelem.Size() > 0 && activeelem.Test(ei.Nr())))
    {
      // cout << firstdof_of_el[elnr] << endl;
      dnums = IntRange(firstdof_of_el[ei.Nr()],firstdof_of_el[ei.Nr()+1]);
      // cout << dnums << endl;
      // getchar();
    }
    else
    {
      dnums.SetSize(0);
      // cout << "no cut" << endl;
    }
  }



  FiniteElement & SFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if (ei.VB() == VOL)
    {
      int elnr = ei.Nr();
      Ngs_Element ngel = ma->GetElement(ei);
      ELEMENT_TYPE eltype = ngel.GetType();
      if (eltype != ET_TRIG)
        throw Exception("can only work with trigs...");
      if (activeelem.Test(elnr))
      {
        return *(new (alloc) SFiniteElement(cuts_on_el[elnr],order,alloc));
      }
      else
      {
        return *dummy;
      }
    }
    else if (ei.VB() == BND)
    {
      return *dummy;
    }
    else
      throw Exception("only VB == VOL and VB == BND implemented");
  }

  SFESpace::SFESpace (shared_ptr<MeshAccess> ama,
                      shared_ptr<CoefficientFunction> a_coef_lset,
                      int aorder,
                      const Flags & flags) :
    FESpace(ama, flags), coef_lset(a_coef_lset), order(aorder)
  {
    type = "sfes";
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    dummy = new DummyFE<ET_TRIG>();
  }

  SFESpace::~SFESpace()
  {
    delete dummy;
  };

  
  void SFESpace::Update()
  {
    // throw Exception ("nothing done yet...");
    LocalHeapMem<100000> lh("SFESpace::Update");

    FESpace::Update();
    int ne=ma->GetNE();
    activeelem.SetSize(ne);

    cuts_on_el.SetSize(ne);

    activeelem.Clear();
    firstdof_of_el.SetSize(ne+1);
    ndof = 0;

    for (int elnr = 0; elnr < ne; ++elnr)
    {
      HeapReset hr(lh);
      Ngs_Element ngel = ma->GetElement(ElementId(VOL,elnr));
      ELEMENT_TYPE eltype = ngel.GetType();
      if (eltype != ET_TRIG)
        throw Exception("only trigs right now...");

      ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);


      Vec<2> ps [] = { Vec<2>(0.0,0.0),
                       Vec<2>(1.0,0.0),
                       Vec<2>(0.0,1.0)};

      IntegrationPoint ip1(ps[0](0),ps[0](1));
      IntegrationPoint ip2(ps[1](0),ps[1](1));
      IntegrationPoint ip3(ps[2](0),ps[2](1));
      MappedIntegrationPoint<2,2> mip1(ip1,eltrans);
      MappedIntegrationPoint<2,2> mip2(ip2,eltrans);
      MappedIntegrationPoint<2,2> mip3(ip3,eltrans);

      double lsets[] = { coef_lset->Evaluate(mip1),
                         coef_lset->Evaluate(mip2),
                         coef_lset->Evaluate(mip3)};
      // cout << lsets[0] << endl;
      // cout << lsets[1] << endl;
      // cout << lsets[2] << endl << endl;

      Array<Vec<2>> cuts(0);

      Array<IVec<2>> edges = { IVec<2>(0,1), IVec<2>(0,2), IVec<2>(1,2)};
      for (auto e: edges)
      {
        const double lset1 = lsets[e[0]];
        const double lset2 = lsets[e[1]];
        if ((lset1>0 && lset2>0) || (lset1<0 && lset2<0))
          continue;
        double cut_scalar = -lsets[e[0]]/(lsets[e[1]]-lsets[e[0]]);
        Vec<2> cut_pos = (1.0 - cut_scalar) * ps[e[0]] + cut_scalar * ps[e[1]];
        cuts.Append(cut_pos);
      }

      firstdof_of_el[elnr] = ndof;

      if (cuts.Size() > 0)
      {
        if (cuts.Size() == 1)
          throw Exception("error: only one cut?!");
        activeelem.SetBitAtomic(elnr);
        ndof += order + 1;
        cuts_on_el[elnr].Col(0) = cuts[0];
        cuts_on_el[elnr].Col(1) = cuts[1];
      }
      // getchar();
      // const double absdet = mip.GetJacobiDet();
      // const double h = sqrt(absdet);
      // ScalarFieldEvaluator * lset_eval_p = NULL;
    }
    firstdof_of_el[ne] = ndof;


  }

}

