#include "xFESpace.hpp"
#include "xfemdiffops.hpp"
using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{
  
  shared_ptr<FESpace> XFESpace::Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    if (ma->GetDimension() == 2)
      return make_shared<T_XFESpace<2> >(ma,flags);
    else
      return make_shared<T_XFESpace<3> >(ma,flags);
  }

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

    auto xstdfes = dynamic_pointer_cast<XStdFESpace>(gf->GetFESpace());
    if (!xstdfes)
      throw Exception("cast failed: not an XStdFESpace");
    auto xfes = xstdfes->GetXFESpace();
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
  }

  template <int D>
  T_XFESpace<D> :: T_XFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : XFESpace (ama, flags)
  {
    // cout << "Constructor of XFESpace begin" << endl;
    vmax = flags.GetNumFlag("vmax",1e99);

    ref_lvl_space = (int) flags.GetNumFlag("ref_space",0);

    if (flags.GetDefineFlag("trace"))
      trace = true;
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpX<D,DIFFOPX::EXTEND_GRAD>>>();
  }


  template <int D>
  SymbolTable<shared_ptr<DifferentialOperator>>
  T_XFESpace<D> :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    switch (ma->GetDimension())
    {
    case 1:
        throw Exception("dim==1 not implemented"); break;
    case 2:
      additional.Set ("extend", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND>>> ()); 
      additional.Set ("pos", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS>>> ()); 
      additional.Set ("neg", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG>>> ()); 
      additional.Set ("extendgrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::EXTEND_GRAD>>> ());
      additional.Set ("posgrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RPOS_GRAD>>> ()); 
      additional.Set ("neggrad", make_shared<T_DifferentialOperator<DiffOpX<2,DIFFOPX::RNEG_GRAD>>> ()); break;
    case 3:
      additional.Set ("extend", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND>>> ());
      additional.Set ("pos", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS>>> ());
      additional.Set ("neg", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG>>> ());
      additional.Set ("extendgrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::EXTEND_GRAD>>> ());
      additional.Set ("posgrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RPOS_GRAD>>> ());
      additional.Set ("neggrad", make_shared<T_DifferentialOperator<DiffOpX<3,DIFFOPX::RNEG_GRAD>>> ()); break;
    default:
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

    TableCreator<int> creator;
    for ( ; !creator.Done(); creator++)
    {
// #pragma omp parallel
      {
        LocalHeap llh(lh.Split());
// #pragma omp for schedule(static)
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(llh);

          Ngs_Element ngel = ma->GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();

          ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), llh);
        
          IntegrationPoint ip(0.0);
          MappedIntegrationPoint<D,D> mip(ip,eltrans);
          const double absdet = mip.GetJacobiDet();
          const double h = D==2 ? sqrt(absdet) : cbrt(absdet);
          ScalarFieldEvaluator * lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,llh);
          CompositeQuadratureRule<D> cquad;
          auto xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_p, 
                                                         cquad, llh, 0, 0, ref_lvl_space, 0);
          xgeom->SetDistanceThreshold(2.0*(h+1.0*vmax));
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

          if (dt == IF)// IsElementCut ?
          {
            activeelem.Set(elnr);
            Array<int> basednums;
            basefes->GetDofNrs(ElementId(VOL,elnr),basednums);
            for (int k = 0; k < basednums.Size(); ++k)
            {
              activedofs.Set(basednums[k]);
// #pragma omp critical(creatoraddel)
              creator.Add(elnr,basednums[k]);
            }
          }
        }
      }
    }
    el2dofs = make_shared<Table<int>>(creator.MoveTable());

    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
    {
      for (int selnr = 0; selnr < nse; ++selnr)
      {
        HeapReset hr(lh);
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

        if (dt == IF) // IsFacetCut(ed)
        {
          Array<int> basednums;
          basefes->GetDofNrs(ei,basednums);
          if (basednums.Size())
            activeselem.Set(selnr);
          for (int k = 0; k < basednums.Size(); ++k)
          {
            activedofs.Set(basednums[k]); // might be twice, but who cares..
            creator2.Add(selnr,basednums[k]);
          }
        }
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

  NumProcXToNegPos::NumProcXToNegPos (shared_ptr<PDE> apde, const Flags & flags)
  {
    gfxstd = apde->GetGridFunction (flags.GetStringFlag ("xstd_gridfunction"));
    gfnegpos = apde->GetGridFunction (flags.GetStringFlag ("negpos_gridfunction"));
  }
  void NumProcXToNegPos::Do (LocalHeap & lh)  {
    XFESpace::XToNegPos(gfxstd, gfnegpos);
  }

  static RegisterNumProc<NumProcXToNegPos> npxtonegpos("xtonegpos");




  XStdFESpace::XStdFESpace (shared_ptr<MeshAccess> ama, 		   
                          const Array<shared_ptr<FESpace>> & aspaces,
                          const Flags & flags)
    : CompoundFESpace(ama, aspaces, flags)
  {
    name="XStdFESpace";

    dynamic_pointer_cast<XFESpace>(spaces[1])->SetBaseFESpace(spaces[0]);

    shared_ptr<ConstantCoefficientFunction> one = make_shared<ConstantCoefficientFunction>(1);
    if (ma->GetDimension() == 2)
    {
      // integrator[VOL] = make_shared<XVisIntegrator<2> > (one);
      // evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpEvalX<2>>>();
      // flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradX<2>>>();
    }
    else
    {
      // integrator[VOL] = make_shared<XVisIntegrator<3> >(one);
      // evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpEvalX<3>>>();
      // flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradX<3>>>();
    }
  }


  ///SmoothingBlocks for good preconditioned iterative solvers
  shared_ptr<Table<int>> XStdFESpace::CreateSmoothingBlocks (const Flags & precflags) const
  {
    Table<int> * it;

    int nv = ma->GetNV();
    int ne = ma->GetNE();
    // int nf = ma->GetNFaces();
    int ned = ma->GetNEdges();      
    Array<int> dnums, dnums2, dnums3; //, verts, edges;

    // bool stdblock=false;   //put all std. dofs in one block + one block for each xfem dof

    // bool blocksystem=false;   //put all std. dofs in one block and all xfem dofs in another block (two large blocks)
    bool vertexpatch=precflags.GetDefineFlag("vertexpatches");   //put all dofs of elements around a vertex together
    bool elementpatch=precflags.GetDefineFlag("elementpatches");  //put all dofs of one element together
    bool edgepatch=precflags.GetDefineFlag("edgepatches");     //put all dofs of neighbouring elements together
    bool vertexdofs=precflags.GetDefineFlag("vertexdofs");    //put all dofs of a vertex together
    bool interfacepatch= precflags.GetDefineFlag("interfacepatch");//put all dofs located at (all) cut elements together (one large block)
    bool jacobi= precflags.GetDefineFlag("jacobi");        //each degree of freedom is one block

    bool dgjumps= precflags.GetDefineFlag("dgjumps");        //each degree of freedom is one block

    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    bool subassembled = precflags.GetDefineFlag("subassembled");
    COUPLING_TYPE dof_mode = eliminate_internal? (subassembled? WIREBASKET_DOF : EXTERNAL_DOF) : ANY_DOF;
    BitArray filter;
    GetFilteredDofs(dof_mode, filter, true);
    FilteredTableCreator creator(&filter);

    // FilteredTableCreator creator(GetFreeDofs());
    *testout << " dof_mode = " << dof_mode << endl;
    *testout << " filter = " << filter << endl;

    *testout << "*GetFreeDofs(): " << endl << *GetFreeDofs() << endl;

    if ( ma->GetDimension() == 2)
    {
      {
        for ( ; !creator.Done(); creator++)
        {      
          int basendof = spaces[0]->GetNDof();
          int xndof = spaces[1]->GetNDof();
          int offset = 0;

          if (jacobi)
          {
            for (int i = 0; i < basendof+xndof ; ++i)
            {
              creator.Add(i,i);
            }
            offset += basendof+xndof;
          }

          if (vertexdofs)
          {
            for (int i = 0; i < nv; i++)
            {
              GetVertexDofNrs(i,dnums);
              for (int j = 0; j < dnums.Size(); ++j)
                if (filter.Test(dnums[j]))
                  if (dnums[j] < basendof)
                    creator.Add(offset+dnums[j], dnums);
            }
            offset += nv;
          }

          if (vertexpatch)
          {

            for (int i = 0; i < nv; i++)
            {
              GetVertexDofNrs(i,dnums);
              for (int j = 0; j < dnums.Size(); ++j)
                for (int k = 0; k < dnums.Size(); ++k)
                  if (filter.Test(dnums[j]))
                    if (dnums[j] < basendof)
                      creator.Add(offset+i,dnums[k]);
            }
            
            for (int i = 0; i < ned; i++)
            {
              Ng_Node<1> edge = ma->GetNode<1> (i);

              GetVertexDofNrs(edge.vertices[0],dnums);
              GetVertexDofNrs(edge.vertices[1],dnums2);
              GetEdgeDofNrs(i,dnums3);
              
              for (int j = 0; j < dnums.Size(); ++j)
              {
                for (int k = 0; k < dnums2.Size(); ++k)
                {
                  if (filter.Test(dnums[j]))
                    if (dnums[j] < basendof)
                      creator.Add(offset+edge.vertices[0],dnums2[k]);
                  if (filter.Test(dnums2[k]))
                    if (dnums2[k] < basendof)
                      creator.Add(offset+edge.vertices[1],dnums[j]);
                }
              }

              
              for (int j = 0; j < dnums.Size(); ++j)
                if (filter.Test(dnums[j]))
                  if (dnums[j] < basendof)
                    creator.Add(offset+edge.vertices[0],dnums3);

              for (int j = 0; j < dnums2.Size(); ++j)
                if (filter.Test(dnums2[j]))
                  if (dnums2[j] < basendof)
                    creator.Add(offset+edge.vertices[1],dnums3);
            }

            for (int elnr = 0; elnr < ne; ++elnr)
            {
              Ng_Node<2> cell = ma->GetNode<2> (elnr);
              GetInnerDofNrs(elnr,dnums);
              for (int i = 0; i < dnums.Size(); ++i)
              {
                GetVertexDofNrs(cell.vertices[0],dnums);
                GetVertexDofNrs(cell.vertices[1],dnums);
                GetVertexDofNrs(cell.vertices[2],dnums);
              }
            }
              
            // }            
            offset += nv;
          }

          if (elementpatch)
          {
            for (int elnr = 0; elnr < ne; ++elnr)
            {
              GetDofNrs(elnr,dnums);
              for (int i = 0; i < dnums.Size(); ++i)
              {
                creator.Add(offset, dnums[i]);
              }
              offset++;
            }
          }


          if (edgepatch)
          {
            int nf = ma->GetNFacets();
            Array<int> elnums;
            Array<int> fnums;
            Array<int> ednums;
            Array<int> vnums;
            shared_ptr<T_XFESpace<2> > xfes22 = NULL;
            shared_ptr<T_XFESpace<3> > xfes33 = NULL;
            const int sD = ma->GetDimension();
            if (sD == 2)
              xfes22 = dynamic_pointer_cast<T_XFESpace<2> >(spaces[1]);
            else
              xfes33 = dynamic_pointer_cast<T_XFESpace<3> >(spaces[1]);
            // Array<int> elnums;
            for (int i = 0; i < nf; i++)
            {
              ma->GetFacetElements(i,elnums);
              /*
              bool cutedge = false;
              for (int k = 0; k < elnums.Size(); ++k)
              {
                int elnr = elnums[k];
                bool elcut = sD ? (spacetime ? xfes23->IsElementCut(elnr) 
                                   : xfes22->IsElementCut(elnr) )
                  :(spacetime ? xfes34->IsElementCut(elnr) 
                    : xfes33->IsElementCut(elnr) );

                if (elcut)
                  cutedge = true;
                else
                  break;
              }

              if (!cutedge)
                continue;
              */
              Ng_Node<1> edge = ma->GetNode<1> (i);
              int v1 = edge.vertices[0];
              int v2 = edge.vertices[1];

              GetEdgeDofNrs(nf,dnums);
              creator.Add(offset, dnums);

              GetVertexDofNrs(edge.vertices[0],dnums);
              creator.Add(offset, dnums);

              GetVertexDofNrs(edge.vertices[1],dnums);
              creator.Add(offset, dnums);
              
              for (int k = 0; k < elnums.Size(); ++k)
              {
                int elnr = elnums[k];
                ma->GetElEdges(elnr,ednums);
                for (int l = 0; l < ednums.Size(); ++l)
                  if (ednums[l] != nf)
                  {
                    GetEdgeDofNrs(ednums[l],dnums);
                    creator.Add(offset, dnums);
                  }

                ma->GetElVertices(elnr,vnums);
                for (int l = 0; l < vnums.Size(); ++l)
                  if (vnums[l] != v1 && vnums[l] != v2)
                  {
                    GetVertexDofNrs(vnums[l],dnums);
                    creator.Add(offset, dnums);
                  }
                    
                GetInnerDofNrs(elnr,dnums);
                creator.Add(offset, dnums);
              }
              offset++;
            }
          }

          if (interfacepatch)
          {
            BitArray mark(GetNDof());
            mark.Clear();
            int epcnt = 0;
            for (int elnr = 0; elnr < ne; ++elnr)
            {
              shared_ptr<T_XFESpace<2> > xfes22 = NULL;
              shared_ptr<T_XFESpace<3> > xfes33 = NULL;
              const int sD = ma->GetDimension();
              if (sD == 2)
                xfes22 = dynamic_pointer_cast<T_XFESpace<2> >(spaces[1]);
              else
                xfes33 = dynamic_pointer_cast<T_XFESpace<3> >(spaces[1]);
              bool elcut = sD == 2 ? xfes22->IsElementCut(elnr) : xfes33->IsElementCut(elnr);
              bool isnbelcut = sD == 2 ? xfes22->IsNeighborElementCut(elnr) : xfes33->IsNeighborElementCut(elnr);
              GetDofNrs(elnr,dnums);
              if (elcut || (dgjumps && isnbelcut))
                for (int i = 0; i < dnums.Size(); ++i)
                {
                  if (!mark.Test(dnums[i]))
                  {
                    mark.Set(dnums[i]);
                    creator.Add(offset, dnums[i]);
                    epcnt++;
                  }
                }
        
            }
            std::cout << " epcnt = " << epcnt << std::endl;
          }

        }
        it = creator.GetTable();
      }
    }
    else
    {
      throw Exception("nono not 3D precond. yet");
      // it = creator.GetTable();
    }
    *testout << "smoothingblocks: " << endl << *it << endl;
    return shared_ptr<Table<int>> (it);
  }


  Array<int> * XStdFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    bool bddc = false;

    if (true || flags.GetDefineFlag("subassembled"))
    {
      cout << "creating bddc-coarse grid(vertices)" << endl;
      Array<int> & clusters = *new Array<int> (GetNDof());
      clusters = 0;
      // int nv = ma->GetNV();
      // for (int i = 0; i < nv; i++)
      //   if (!IsDirichletVertex(i))
      //     clusters[i] = 1;		
      return &clusters;	
    }
    else
    {
      if(bddc)
      {
        Array<int> & clusters = *new Array<int> (GetNDof());
        clusters = 0;

        Array<int> dnums;
        // int nfa = ma->GetNFacets();
        int nv = ma->GetNV();

        for (int i = 0; i < nv; i++)
        {
          // basefes->GetVertexDofNrs (i, dnums);
          clusters[nv /*dnums[0]*/] = 1;
        }

        const BitArray & freedofs = *GetFreeDofs();
        for (int i = 0; i < freedofs.Size(); i++)
          if (!freedofs.Test(i)) clusters[i] = 0;
        *testout << "XStdFESpace, dsc = " << endl << clusters << endl;
        return &clusters;
      }
      else // put all degrees of freedoms at the interface on the coarse grid
      {

        Array<int> dnums;
        BitArray mark(GetNDof());
        mark.Clear();
        int epcnt = 0;
        for (int elnr = 0; elnr < ma->GetNE(); ++elnr)
        {
          shared_ptr<T_XFESpace<2> > xfes22 = NULL;
          shared_ptr<T_XFESpace<3> > xfes33 = NULL;
          const int sD = ma->GetDimension();
          if (sD == 2)
            xfes22 = dynamic_pointer_cast<T_XFESpace<2> >(spaces[1]);
          else
            xfes33 = dynamic_pointer_cast<T_XFESpace<3> >(spaces[1]);
          bool elcut = sD == 2 ? xfes22->IsElementCut(elnr) : xfes33->IsElementCut(elnr);

          GetDofNrs(elnr,dnums);
          if (elcut)
            for (int i = 0; i < dnums.Size(); ++i)
            {
              if (!mark.Test(dnums[i]))
              {
                mark.Set(dnums[i]);
                epcnt++;
              }
            }
        }
        std::cout << " epcnt = " << epcnt << std::endl;

        Array<int> & clusters = *new Array<int> (GetNDof());
        clusters = 0;

        for (int i = 0; i < GetNDof(); i++)
        {
          if (mark.Test(i)) clusters[i] = 1;
        }

        const BitArray & freedofs = *GetFreeDofs();
        for (int i = 0; i < freedofs.Size(); i++)
          if (!freedofs.Test(i)) clusters[i] = 0;
        *testout << "XStdFESpace, dsc = " << endl << clusters << endl;
        // cout << "XStdFESpace, dsc = " << endl << clusters << endl;
        // getchar();
        return &clusters;
      }
    }
  }


  namespace xfespace_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
  
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("xfespace", XFESpace::Create);
      GetFESpaceClasses().AddFESpace ("xstdfespace", XStdFESpace::Create);
      GetFESpaceClasses().AddFESpace ("xh1fespace", XStdFESpace::Create); //backward compatibility
    }
  
    Init init;

  }



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
    LocalHeap lh(10000000,"SFESpace::GetFE");
    if (ei.VB() == VOL)
    {
      int elnr = ei.Nr();
      Ngs_Element ngel = ma->GetElement(elnr);
      ELEMENT_TYPE eltype = ngel.GetType();
      if (eltype != ET_TRIG)
        throw Exception("can only work with trigs...");
      if (activeelem.Test(elnr))
      {
        return *(new (lh) SFiniteElement(cuts_on_el[elnr],order,lh));
      }
      else
        return *(new (lh) DummyFE<ET_TRIG>());
    }
    else if (ei.VB() == BND)
    {
      return *(new (lh) DummyFE<ET_SEGM>());
    }
    else
      throw Exception("only VB == VOL and VB == BND implemented");
  }

  SFESpace::SFESpace (shared_ptr<MeshAccess> ama,
                      shared_ptr<CoefficientFunction> a_coef_lset,
                      int aorder,
                      const Flags & flags):
    FESpace(ama, flags), coef_lset(a_coef_lset), order(aorder)
  {
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
  }

  void SFESpace::Update(LocalHeap & lh)
  {
    // throw Exception ("nothing done yet...");

    FESpace::Update(lh);
    int ne=ma->GetNE();
    activeelem.SetSize(ne);

    cuts_on_el.SetSize(ne);

    activeelem.Clear();
    firstdof_of_el.SetSize(ne+1);
    ndof = 0;

    for (int elnr = 0; elnr < ne; ++elnr)
    {
      HeapReset hr(lh);
      Ngs_Element ngel = ma->GetElement(elnr);
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

      Array<INT<2>> edges = { INT<2>(0,1), INT<2>(0,2), INT<2>(1,2)};
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
        activeelem.Set(elnr);
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
