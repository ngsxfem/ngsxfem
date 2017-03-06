#include "xFESpace.hpp"
#include "xfemdiffops.hpp"
using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{
  void XFESpace :: CleanUp ()
  {
    //empty
  }



  void XFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if (ei.VB() == VOL)
    {
      if (activeelem.Size() > 0 && activeelem.Test(ei.Nr()))
        dnums = (*el2dofs)[ei.Nr()];
      else
        dnums.SetSize(0);
    }
    else
    {
      if (activeselem.Size() > 0 && activeselem.Test(ei.Nr()))
        dnums = (*sel2dofs)[ei.Nr()];
      else
        dnums.SetSize(0);
    }
  }

  void XFESpace :: GetDomainNrs (ElementId ei, Array<DOMAIN_TYPE> & domnums) const
  {
    if (ei.VB() == VOL)
    {
      if (activeelem.Test(ei.Nr()))
      {
        FlatArray<int> dofs = (*el2dofs)[ei.Nr()];
        domnums.SetSize(dofs.Size());
        for (int i = 0; i < dofs.Size(); ++i)
        {
          domnums[i] = domofdof[dofs[i]];
        }
      }
      else
        domnums.SetSize(0);
    }
    else if (ei.VB() == BND)
    {
      if (activeselem.Test(ei.Nr()))
      {
        FlatArray<int> dofs = (*sel2dofs)[ei.Nr()];
        domnums.SetSize(dofs.Size());
        for (int i = 0; i < dofs.Size(); ++i)
        {
          domnums[i] = domofdof[dofs[i]];
        }
      }
      else
        domnums.SetSize(0);
    }
    else
      throw Exception("can only handle VOL and BND");
  }

  void XFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;

    for (int i = 0; i < basedof2xdof.Size(); ++i)
    {
      const int dof = basedof2xdof[i];
      if (dof != -1)
      {
        // if (trace)
        ctofdof[dof] = basefes->GetDofCouplingType(i); //INTERFACE_DOF; //
        // else
        //   ctofdof[dof] = INTERFACE_DOF; //
      }
    }

    if (trace && ma->GetDimension() == 3)
    // face bubbles on the outer part of the band will be local dofs... (for static cond.)
    {
      for (int facnr = 0; facnr < ma->GetNFaces(); ++facnr)
      {
        Array<int> elnums;
        ma->GetFaceElements (facnr, elnums);
        int cutels = 0;
        for (auto elnr : elnums)
        {
          if (activeelem.Test(elnr))
            cutels++;
        }
        if (cutels<2)
        {
          Array<int> facedofs;
          basefes->GetFaceDofNrs (facnr, facedofs);
          for (auto basedof : facedofs)
          {
            const int dof = basedof2xdof[basedof];
            if (dof != -1)
              ctofdof[dof] = LOCAL_DOF;
          }
        }
      }
    }
    *testout << "XFESpace, ctofdof = " << endl << ctofdof << endl;
    // cout << "XFESpace, ctofdof = " << endl << ctofdof << endl;
    // getchar();

  }


  void XFESpace::XToNegPos(shared_ptr<GridFunction> gf, shared_ptr<GridFunction> gf_neg_pos)
  {
    shared_ptr<GridFunction> gf_neg = gf_neg_pos->GetComponent(0);
    BaseVector & bv_neg = gf_neg->GetVector();
    FlatVector<> vneg = bv_neg.FVDouble();

    shared_ptr<GridFunction> gf_pos = gf_neg_pos->GetComponent(1);
    BaseVector & bv_pos = gf_pos->GetVector();
    FlatVector<> vpos = bv_pos.FVDouble();

    shared_ptr<GridFunction> gf_base = gf->GetComponent(0);
    BaseVector & bv_base = gf_base->GetVector();
    FlatVector<> vbase = bv_base.FVDouble();

    shared_ptr<GridFunction> gf_x = gf->GetComponent(1);
    BaseVector & bv_x = gf_x->GetVector();
    FlatVector<> vx = bv_x.FVDouble();

    auto xstdfes = dynamic_pointer_cast<CompoundFESpace>(gf->GetFESpace());
    if (!xstdfes)
      throw Exception("cast failed: not a CompoundFESpace");
    auto xfes = dynamic_pointer_cast<XFESpace>((*xstdfes)[1]);
    if (!xfes)
      throw Exception("cast failed: not an XFESpace");

    const int basendof = vneg.Size();
    for (int i = 0; i < basendof; ++i)
    {
      vneg(i) = vbase(i);
      vpos(i) = vbase(i);
      const int xdof = xfes->GetXDofOfBaseDof(i);
      if (xdof != -1)
      {
        if (xfes->GetDomOfDof(xdof) == POS)
          vpos(i) += vx(xdof);
        else
          vneg(i) += vx(xdof);
      }
    }
  }


  template <int D>
  T_XFESpace<D> :: T_XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
                               shared_ptr<CutInformation> cutinfo, const Flags & flags)
    : XFESpace (ama, basefes, cutinfo, flags)
  {
    // cout << "Constructor of XFESpace begin" << endl;
    vmax = flags.GetNumFlag("vmax",1e99);

    ref_lvl_space = (int) flags.GetNumFlag("ref_space",0);

    if (flags.GetDefineFlag("trace"))
      trace = true;
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND_GRAD>>>();

    private_cutinfo = false;
  }

  template <int D>
  T_XFESpace<D> :: T_XFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> basefes,
                               shared_ptr<CoefficientFunction> lset, const Flags & flags)
    : XFESpace (ama, basefes, lset, flags)
  {
    // cout << "Constructor of XFESpace begin" << endl;
    vmax = flags.GetNumFlag("vmax",1e99);

    ref_lvl_space = (int) flags.GetNumFlag("ref_space",0);

    if (flags.GetDefineFlag("trace"))
      trace = true;
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND_GRAD>>>();

    cutinfo = make_shared<CutInformation>(ma);
  }

  template <int D>
  SymbolTable<shared_ptr<DifferentialOperator>>
  T_XFESpace<D> :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    switch (ma->GetDimension())
    {
    case 1 :
      throw Exception("dim==1 not implemented"); break;
    case 2 :
      additional.Set ("extend", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND>>> ());
      additional.Set ("pos", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS>>> ());
      additional.Set ("neg", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG>>> ());
      additional.Set ("extendgrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND_GRAD>>> ());
      additional.Set ("posgrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS_GRAD>>> ());
      additional.Set ("neggrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG_GRAD>>> ()); break;
    case 3 :
      additional.Set ("extend", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND>>> ());
      additional.Set ("pos", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS>>> ());
      additional.Set ("neg", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG>>> ());
      additional.Set ("extendgrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND_GRAD>>> ());
      additional.Set ("posgrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS_GRAD>>> ());
      additional.Set ("neggrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG_GRAD>>> ()); break;
    default :
      ;
    }
    return additional;
  }


  template <int D>
  T_XFESpace<D> :: ~T_XFESpace ()
  {
    CleanUp();
    // if (eval_lset) delete eval_lset;
  }


  void IterateRange (int ne, LocalHeap & clh,
                     const function<void(int,LocalHeap&)> & func)
  {
    if (task_manager)
    {
      SharedLoop2 sl(ne);
      task_manager -> CreateJob
        ( [&] (const TaskInfo & ti)
      {
        LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
        for (int elnr : sl)
        {
          HeapReset hr(lh);
          func (elnr,lh);
        }

      } );
    }
    else
    {
      for (int elnr = 0; elnr < ne; elnr++)
      {
        HeapReset hr(clh);
        func (elnr,clh);
      }
    }
  }

  template <int D>
  void T_XFESpace<D> :: Update(LocalHeap & lh)
  {
    if ( basefes == NULL )
    {
      cout << " T_XFESpace<" <<  D << "> no basefes, Update postponed " << endl;
      ndof = 0.0;
      return;
    }
    if ( coef_lset == NULL )
    {
      cout << " T_XFESpace<" <<  D << "> no lset, Update postponed " << endl;
      ndof = 0.0;
      return;
    }
    CleanUp();

    if (private_cutinfo)
    {
      cout << IM(4) << " Calling cutinfo-Update from within XFESpace-Update " << endl;
      cutinfo->Update(coef_lset,lh);
    }

    static Timer timer ("XFESpace::Update");
    RegionTimer reg (timer);

    order_space = basefes->GetOrder();

    FESpace::Update(lh);

    int ne=ma->GetNE();
    int nedges=ma->GetNEdges();
    int nf=ma->GetNFaces();
    int nv=ma->GetNV();
    int nse=ma->GetNSE();

    activeelem.SetSize(ne);
    activeselem.SetSize(nse);
    activeelem.Clear();
    activeselem.Clear();

    BitArray activedofs(basefes->GetNDof());
    activedofs.Clear();

    domofel.SetSize(ne);
    domofsel.SetSize(nse);

    static int first = -1;
    first++;

    Array<double> kappa_pos(ne);
    BitArray element_most_pos(ne);
    element_most_pos.Clear();


    IterateRange
      (ne, lh,
      [&] (int elnr, LocalHeap & lh)
    {
      Ngs_Element ngel = ma->GetElement(elnr);
      ELEMENT_TYPE eltype = ngel.GetType();
      ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
      // IntegrationPoint ip(0.0);
      // MappedIntegrationPoint<D,D> mip(ip,eltrans);
      // const double absdet = mip.GetJacobiDet();
      // const double h = D==2 ? sqrt(absdet) : cbrt(absdet);
      ScalarFieldEvaluator * lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,lh);
      CompositeQuadratureRule<D> cquad;
      auto xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_p,
                                                     cquad, lh, 0, 0, ref_lvl_space, 0);
      // xgeom->SetDistanceThreshold(2.0*(h+1.0*vmax));
      DOMAIN_TYPE dt = xgeom->MakeQuadRule();

      QuadratureRule<D> & pquad =  cquad.GetRule(POS);
      QuadratureRule<D> & nquad =  cquad.GetRule(NEG);
      double pospart_vol = 0.0;
      double negpart_vol = 0.0;
      for (int i = 0; i < pquad.Size(); ++i)
        pospart_vol += pquad.weights[i];
      for (int i = 0; i < nquad.Size(); ++i)
        negpart_vol += nquad.weights[i];
      if (pospart_vol > negpart_vol)
        element_most_pos.Set(elnr);

      domofel[elnr] = dt;

      if (dt == IF)   // IsElementCut ?
      {
        activeelem.Set(elnr);
        Array<int> basednums;
        basefes->GetDofNrs(ElementId(VOL,elnr),basednums);
        for (int k = 0; k < basednums.Size(); ++k)
          activedofs.Set(basednums[k]);
      }
    });

    TableCreator<int> creator;
    for (; !creator.Done(); creator++)
    {
      for (int elnr = 0; elnr < ne; ++elnr)
      {
        if (! activeelem[elnr]) continue;

        Array<int> basednums;
        basefes->GetDofNrs(ElementId(VOL,elnr),basednums);
        for (int k = 0; k < basednums.Size(); ++k)
          creator.Add(elnr,basednums[k]);
      }
    }
    el2dofs = make_shared<Table<int>>(creator.MoveTable());

    IterateRange
      (nse, lh,
      [&] (int selnr, LocalHeap & lh)
    {
      ElementId ei(BND,selnr);
      Ngs_Element ngel = ma->GetElement(ei);
      ELEMENT_TYPE eltype = ngel.GetType();

      ElementTransformation & seltrans = ma->GetTrafo (ei, lh);

      ScalarFieldEvaluator * lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,seltrans,lh);

      CompositeQuadratureRule<D-1> cquad;
      auto xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_p,
                                                     cquad, lh, 0, 0, ref_lvl_space, 0);
      DOMAIN_TYPE dt = xgeom->MakeQuadRule();

      Array<int> fnums;
      ma->GetSElFacets(selnr,fnums);
      // int ed = fnums[0];

      domofsel[selnr] = dt;

      if (dt == IF)    // IsFacetCut(ed)
      {
        Array<int> basednums;
        basefes->GetDofNrs(ei,basednums);
        if (basednums.Size())
          activeselem.Set(selnr);
        for (int k = 0; k < basednums.Size(); ++k)
          activedofs.Set(basednums[k]);    // might be twice, but who cares..
      }
    });

    TableCreator<int> creator2;
    for (; !creator2.Done(); creator2++)
    {
      for (int selnr = 0; selnr < nse; ++selnr)
      {
        if (! activeselem[selnr]) continue;
        Array<int> basednums;
        basefes->GetDofNrs(ElementId(BND,selnr),basednums);
        for (int k = 0; k < basednums.Size(); ++k)
          creator2.Add(selnr,basednums[k]);
      }
    }
    sel2dofs = make_shared<Table<int>>(creator2.MoveTable());

    int nbdofs = basefes->GetNDof();
    basedof2xdof.SetSize(nbdofs);
    basedof2xdof = -1;
    ndof = 0;

    for (int i = 0; i < nbdofs; i++)
    {
      if (activedofs.Test(i))
        basedof2xdof[i] = ndof++;
    }

    xdof2basedof.SetSize(ndof);
    xdof2basedof = -1;

    ndof = 0;
    for (int i = 0; i < nbdofs; i++)
    {
      if (activedofs.Test(i))
        xdof2basedof[ndof++] = i;
    }

    for (int i = 0; i < ne; ++i)
    {
      if (activeelem.Test(i))
      {
        FlatArray<int> dofs = (*el2dofs)[i];
        for (int j = 0; j < (*el2dofs)[i].Size(); ++j)
          (*el2dofs)[i][j] = basedof2xdof[dofs[j] ];
      }
    }

    for (int i = 0; i < nse; ++i)
    {
      if (activeselem.Test(i))
      {
        FlatArray<int> dofs = (*sel2dofs)[i];
        for (int j = 0; j < (*sel2dofs)[i].Size(); ++j)
          (*sel2dofs)[i][j] = basedof2xdof[dofs[j] ];
      }
    }

    *testout << " x ndof : " << ndof << endl;

    // domain of dof
    domofdof.SetSize(ndof);
    domofdof = IF;

    if (D==3)
    {
      domofface.SetSize(nf);
      for (int facnr = 0; facnr < nf; ++facnr)
      {
        bool haspos = false;
        // bool hasneg = false;

        Array<int> elnums;
        ma->GetFaceElements (facnr, elnums);

        for (int k = 0; k < elnums.Size(); ++k)
        {
          DOMAIN_TYPE dt_cur = domofel[elnums[k]];
          // if (dt_cur == NEG)
          //   hasneg = true;

          if (dt_cur == POS)
            haspos = true;
        }

        if (haspos)
          domofface[facnr] = NEG;
        else
          domofface[facnr] = POS;

        Array<int> dnums;
        basefes->GetFaceDofNrs(facnr, dnums);
        for (int l = 0; l < dnums.Size(); ++l)
        {
          int xdof = basedof2xdof[dnums[l]];
          if ( xdof != -1)
            domofdof[xdof] = domofface[facnr];
        }
      }
    }

    domofedge.SetSize(nedges);
    for (int edgnr = 0; edgnr < nedges; ++edgnr)
    {
      bool haspos = false;
      // bool hasneg = false;

      Array<int> elnums;
      ma->GetEdgeElements (edgnr, elnums);

      for (int k = 0; k < elnums.Size(); ++k)
      {
        DOMAIN_TYPE dt_cur = domofel[elnums[k]];
        // if (dt_cur == NEG)
        //   hasneg = true;

        if (dt_cur == POS)
          haspos = true;
      }

      if (haspos)
        domofedge[edgnr] = NEG;
      else
        domofedge[edgnr] = POS;

      Array<int> dnums;
      basefes->GetEdgeDofNrs(edgnr, dnums);

      for (int l = 0; l < dnums.Size(); ++l)
      {
        int xdof = basedof2xdof[dnums[l]];
        if ( xdof != -1)
          domofdof[xdof] = domofedge[edgnr];
      }
    }

    nvertdofs = 0;
    domofvertex.SetSize(nv);
    for (int vnr = 0; vnr < nv; ++vnr)
    {
      bool haspos = false;
      // bool hasneg = false;

      Array<int> elnums;
      ma->GetVertexElements (vnr, elnums);

      for (int k = 0; k < elnums.Size(); ++k)
      {
        DOMAIN_TYPE dt_cur = domofel[elnums[k]];
        // if (dt_cur == NEG)
        //   hasneg = true;

        if (dt_cur == POS)
          haspos = true;
      }

      if (haspos)
        domofvertex[vnr] = NEG;
      else
        domofvertex[vnr] = POS;

      Array<int> dnums;
      basefes->GetVertexDofNrs(vnr, dnums);

      for (int l = 0; l < dnums.Size(); ++l)
      {
        int xdof = basedof2xdof[dnums[l]];
        if ( xdof != -1)
        {
          domofdof[xdof] = domofvertex[vnr];
          nvertdofs++;
        }
      }
    }

    domofinner.SetSize(ne);
    for (int elnr = 0; elnr < ne; ++elnr)
    {
      DOMAIN_TYPE dt_here = element_most_pos.Test(elnr) ? POS : NEG;
      domofinner[elnr] = dt_here;
      Array<int> dnums;
      basefes->GetInnerDofNrs(elnr, dnums);
      for (int l = 0; l < dnums.Size(); ++l)
      {
        int xdof = basedof2xdof[dnums[l]];
        if ( xdof != -1)
          domofdof[xdof] = dt_here;
      }
    }

    // domof dof on boundary
    for (int selnr = 0; selnr < nse; ++selnr)
    {
      ElementId ei(BND,selnr);
      DOMAIN_TYPE dt = domofsel[selnr];
      Array<int> dnums;
      basefes->GetDofNrs(ei, dnums);

      for (int i = 0; i < dnums.Size(); ++i)
      {
        const int xdof = basedof2xdof[dnums[i]];
        if (xdof != -1)
        {
          if (dt != IF)
            domofdof[xdof] = dt == POS ? NEG : POS;
        }
      }
    }

    BitArray dofs_with_cut_on_boundary(GetNDof());
    dofs_with_cut_on_boundary.Clear();

    for (int selnr = 0; selnr < nse; ++selnr)
    {
      ElementId ei(BND,selnr);
      DOMAIN_TYPE dt = domofsel[selnr];
      if (dt!=IF) continue;

      Array<int> dnums;
      GetDofNrs(ei, dnums);

      for (int i = 0; i < dnums.Size(); ++i)
      {
        const int xdof = dnums[i];
        dofs_with_cut_on_boundary.Set(xdof);
      }
    }

    UpdateCouplingDofArray();
    FinalizeUpdate (lh);

    dirichlet_dofs.SetSize (GetNDof());
    dirichlet_dofs.Clear();

    for (int i = 0; i < basedof2xdof.Size(); ++i)
    {
      const int dof = basedof2xdof[i];
      if (dof != -1 && basefes->IsDirichletDof(i))
        if (dofs_with_cut_on_boundary.Test(dof))
          dirichlet_dofs.Set (dof);
    }

    free_dofs->SetSize (GetNDof());
    *free_dofs = dirichlet_dofs;
    free_dofs->Invert();

    *testout << "ndof = " << ndof << endl;
    *testout << "basedof2xdof = " << basedof2xdof << endl;
    *testout << "el2dofs = " << *el2dofs << endl;
    *testout << "sel2dofs = " << *sel2dofs << endl;
    *testout << "domain of dofs = " << domofdof << endl;

    *testout << " basefes = " << basefes << endl;
    *testout << "basefes -> free_dofs = " << basefes->GetFreeDofs() << endl;
    if (basefes->GetFreeDofs())
      *testout << "*(basefes -> free_dofs) = " << *(basefes->GetFreeDofs()) << endl;
    *testout << "free_dofs = " << *free_dofs << endl;
  }

  template <int D>
  FiniteElement & T_XFESpace<D> :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if ( ( ei.VB() == BND && activeselem.Test(ei.Nr()) )
         || ( ei.VB() == VOL && activeelem.Test(ei.Nr()) ))
    {
      Array<DOMAIN_TYPE> domnrs;
      GetDomainNrs(ei,domnrs);
      return *(new (alloc) XFiniteElement(basefes->GetFE(ei,alloc),domnrs,alloc));
    }
    else
    {
      DOMAIN_TYPE dt = domofsel[ei.Nr()];
      Ngs_Element ngsel = ma->GetElement(ei);
      ELEMENT_TYPE eltype = ngsel.GetType();
      return *(new (alloc) XDummyFE(dt,eltype));
    }
  }

  template class T_XFESpace<2>;
  template class T_XFESpace<3>;

}

