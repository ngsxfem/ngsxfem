
#include "xFESpace.hpp"
#include "xfemVisInts.hpp"
using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{
  
  template <int D, int SD>
  XFESpace<D,SD> :: XFESpace (const MeshAccess & ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of XFESpace begin" << endl;
    spacetime = flags.GetDefineFlag("spacetime");
    empty = flags.GetDefineFlag("empty");
    if (empty)
        cout << " EMPTY XFESPACE active..." << endl;
    ti.first = flags.GetNumFlag("t0",0.0);
    ti.second = flags.GetNumFlag("t1",1.0);

    vmax = flags.GetNumFlag("vmax",1e99);

    ref_lvl_space = (int) flags.GetNumFlag("ref_space",0);
    ref_lvl_time = (int) flags.GetNumFlag("ref_time",0);

    std::cout << " ref_lvl_space = " << ref_lvl_space << std::endl;
    std::cout << " ref_lvl_time = " << ref_lvl_time << std::endl;

    string eval_lset_str(flags.GetStringFlag ("levelset","lset"));
    eval_lset = new EvalFunction(eval_lset_str);

    /* 
    static ConstantCoefficientFunction one(1);
    if (spacetime)
      integrator = new STXVisIntegrator<D,FUTURE>(&one) ;
    else
      integrator = new XVisIntegrator<D>(&one) ;
    // boundary_integrator = new RobinIntegrator<2> (&one);
    */

    cout << "Constructor of XFESpace end" << endl;
    // static ConstantCoefficientFunction one(1);
    // integrator = new MassIntegrator<D> (&one);
  }
    
  
  template <int D, int SD>
  void XFESpace<D,SD> :: CleanUp ()
  {

    if (el2dofs) delete el2dofs; 
    if (sel2dofs) delete sel2dofs; 
  }

  template <int D, int SD>
  XFESpace<D,SD> :: ~XFESpace ()
  {
    CleanUp();
    if (eval_lset) delete eval_lset;
  }
  
  template <int D, int SD>
  void XFESpace<D,SD> :: Update(LocalHeap & lh)
  {
    if ( basefes == NULL )
    {
      cout << " no basefes, Update postponed " << endl;
      return;
    }
    CleanUp();

    static Timer timer ("XFESpace::Update");
    RegionTimer reg (timer);

    order_space = basefes->GetOrder();

    if (spacetime)
      order_time = dynamic_cast<const SpaceTimeFESpace *>(basefes)->OrderTime();
    
    FESpace::Update(lh);

    int ne=ma.GetNE();
    int nedges=ma.GetNEdges();
    int nf=ma.GetNFaces();
    int nv=ma.GetNV();
    int nse=ma.GetNSE();

    activeelem.SetSize(ne);    
    activeselem.SetSize(nse);    
    activeelem.Clear();
    activeselem.Clear();

    BitArray activedofs(basefes->GetNDof());
    activedofs.Clear();

    domofel.SetSize(ne);
    domofsel.SetSize(nse);

    ELEMENT_TYPE et_time = D == SD ? ET_POINT : ET_SEGM;

    static int first = -1;
    first++;

    Array<double> kappa_pos(ne);
    BitArray element_most_pos(ne);
    element_most_pos.Clear();

    TableCreator<int> creator;
    for ( ; !creator.Done(); creator++)
    {
#pragma omp parallel
      {
        LocalHeap llh(lh.Split());
#pragma omp for schedule(static)
        for (int elnr = 0; elnr < ne; ++elnr)
        {
          HeapReset hr(llh);

          Ngs_Element ngel = ma.GetElement(elnr);
          ELEMENT_TYPE eltype = ngel.GetType();

          ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), llh);
        
          IntegrationPoint ip(0.0);
          MappedIntegrationPoint<D,D> mip(ip,eltrans);
          const double absdet = mip.GetJacobiDet();
          const double h = D==2 ? sqrt(absdet) : cbrt(absdet);
          ScalarFieldEvaluator * lset_eval_p = NULL;
          if (spacetime)
            lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,ti,llh);
          else
            lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,llh);
          CompositeQuadratureRule<SD> cquad;
          XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                                cquad, llh, 2*order_space, 2*order_time, ref_lvl_space, ref_lvl_time);
          xgeom->SetDistanceThreshold(2.0*(h+(ti.second-ti.first)*vmax));
          DOMAIN_TYPE dt = xgeom->MakeQuadRule();

          QuadratureRule<SD> & pquad =  cquad.GetRule(POS);
          QuadratureRule<SD> & nquad =  cquad.GetRule(NEG);
          double pospart_vol = 0.0;
          double negpart_vol = 0.0;
          for (int i = 0; i < pquad.Size(); ++i)
            pospart_vol += pquad.weights[i];
          for (int i = 0; i < nquad.Size(); ++i)
            negpart_vol += nquad.weights[i];
          if (pospart_vol > negpart_vol)
            element_most_pos.Set(elnr);

          delete xgeom;

          domofel[elnr] = dt;

          if (dt == IF)// IsElementCut ?
          {
            activeelem.Set(elnr);
            Array<int> basednums;
            basefes->GetDofNrs(elnr,basednums);
            for (int k = 0; k < basednums.Size(); ++k)
            {
              activedofs.Set(basednums[k]);
#pragma omp critical(creatoraddel)
              creator.Add(elnr,basednums[k]);
            }
          }
        }
      }
    }
    el2dofs = creator.GetTable();

    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
    {
      for (int selnr = 0; selnr < nse; ++selnr)
      {
        HeapReset hr(lh);

        Ngs_Element ngel = ma.GetSElement(selnr);
        ELEMENT_TYPE eltype = ngel.GetType();

        ElementTransformation & seltrans = ma.GetTrafo (selnr, BND, lh);

        ScalarFieldEvaluator * lset_eval_p = NULL;
        if (spacetime)
          lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,seltrans,ti,lh);
        else
          lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,seltrans,lh);

        CompositeQuadratureRule<SD-1> cquad;
        XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                              cquad, lh, 2*order_space, 2*order_time, ref_lvl_space, ref_lvl_time);
        DOMAIN_TYPE dt = xgeom->MakeQuadRule();

        delete xgeom;

        Array<int> fnums;
        ma.GetSElFacets(selnr,fnums);
        int ed = fnums[0];

        domofsel[selnr] = dt;

        if (dt == IF) // IsFacetCut(ed)
        {
          activeselem.Set(selnr);
          Array<int> basednums;
          basefes->GetSDofNrs(selnr,basednums);
          for (int k = 0; k < basednums.Size(); ++k)
          {
            activedofs.Set(basednums[k]); // might be twice, but who cares..
            creator2.Add(selnr,basednums[k]);
          }
        }
      }
    }
    sel2dofs = creator2.GetTable();


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
        bool hasneg = false;

        Array<int> elnums;
        ma.GetFaceElements (facnr, elnums);

        for (int k = 0; k < elnums.Size(); ++k)
        {
          DOMAIN_TYPE dt_cur = domofel[elnums[k]];
          if (dt_cur == NEG)
            hasneg = true;

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
      bool hasneg = false;

      Array<int> elnums;
      ma.GetEdgeElements (edgnr, elnums);

      for (int k = 0; k < elnums.Size(); ++k)
      {
        DOMAIN_TYPE dt_cur = domofel[elnums[k]];
        if (dt_cur == NEG)
          hasneg = true;

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

    domofvertex.SetSize(nv);
    for (int vnr = 0; vnr < nv; ++vnr)
    {
      bool haspos = false;
      bool hasneg = false;

      Array<int> elnums;
      ma.GetVertexElements (vnr, elnums);

      for (int k = 0; k < elnums.Size(); ++k)
      {
        DOMAIN_TYPE dt_cur = domofel[elnums[k]];
        if (dt_cur == NEG)
          hasneg = true;

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
          domofdof[xdof] = domofvertex[vnr];
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
        DOMAIN_TYPE dt = domofsel[selnr];
        Array<int> dnums;
        basefes->GetSDofNrs(selnr, dnums);

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
        DOMAIN_TYPE dt = domofsel[selnr];
        if (dt!=IF) continue;

        Array<int> dnums;
        GetSDofNrs(selnr, dnums);

        for (int i = 0; i < dnums.Size(); ++i)
        {
            const int xdof = dnums[i];
            dofs_with_cut_on_boundary.Set(xdof);
        }
    }

    if (empty)
        ndof = 0;

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
    
    free_dofs.SetSize (GetNDof());
    free_dofs = dirichlet_dofs;
    free_dofs.Invert();

    *testout << "ndof = " << ndof << endl;
    *testout << "basedof2xdof = " << basedof2xdof << endl;
    *testout << "el2dofs = " << *el2dofs << endl;
    *testout << "sel2dofs = " << *sel2dofs << endl;
    *testout << "domain of dofs = " << domofdof << endl;

    *testout << " basefes = " << basefes << endl;
    *testout << "basefes -> free_dofs = " << *(basefes->GetFreeDofs()) << endl;

    *testout << "free_dofs = " << free_dofs << endl;
  }
  

  template <int D, int SD>
  void XFESpace<D,SD> :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    if (activeelem.Test(elnr) && !empty)
      dnums = (*el2dofs)[elnr];
    else
      dnums.SetSize(0);
  }
  
  template <int D, int SD>
  void XFESpace<D,SD> :: GetDomainNrs (int elnr, Array<DOMAIN_TYPE> & domnums) const
  {
    if (activeelem.Test(elnr) && !empty)
    {
      FlatArray<int> dofs = (*el2dofs)[elnr];
      domnums.SetSize(dofs.Size());
      for (int i = 0; i < dofs.Size(); ++i)
      {
        domnums[i] = domofdof[dofs[i]];
      }
    }
    else
      domnums.SetSize(0);
  }

  template <int D, int SD>
  void XFESpace<D,SD> :: GetSurfaceDomainNrs (int selnr, Array<DOMAIN_TYPE> & domnums) const
  {
    if (activeselem.Test(selnr) && !empty)
    {
      FlatArray<int> dofs = (*sel2dofs)[selnr];
      domnums.SetSize(dofs.Size());
      for (int i = 0; i < dofs.Size(); ++i)
      {
        domnums[i] = domofdof[dofs[i]];
      }
    }
    else
      domnums.SetSize(0);
  }

  template <int D, int SD>
  void XFESpace<D,SD> :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    ctofdof = WIREBASKET_DOF;

    if (!empty)
        for (int i = 0; i < basedof2xdof.Size(); ++i)
        {
            const int dof = basedof2xdof[i];
            if (dof != -1)
                ctofdof[dof] = INTERFACE_DOF; //basefes->GetDofCouplingType(i);
        }
    *testout << "XFESpace, ctofdof = " << endl << ctofdof << endl;
  }


  template <int D, int SD>
  void XFESpace<D,SD> :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    if (activeselem.Test(selnr) && !empty)
      dnums = (*sel2dofs)[selnr];
    else
      dnums.SetSize(0);
  }
  template <int D, int SD>
  const FiniteElement & XFESpace<D,SD> :: GetFE (int elnr, LocalHeap & lh) const
  {
    static Timer timer ("XFESpace::GetFE");
    RegionTimer reg (timer);

    Ngs_Element ngel = ma.GetElement(elnr);
    ELEMENT_TYPE eltype = ngel.GetType();
    if (!activeelem.Test(elnr))
    {
      DOMAIN_TYPE dt = domofel[elnr];
      return *(new (lh) XDummyFE(dt,eltype));
    }
    else
    {
      Array<DOMAIN_TYPE> domnrs;
      GetDomainNrs(elnr,domnrs);  

      Ngs_Element ngel = ma.GetElement(elnr);
      ELEMENT_TYPE eltype = ngel.GetType();

      ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), lh);
        
      ScalarFieldEvaluator * lset_eval_p = NULL;
      if (spacetime)
        lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,ti,lh);
      else
        lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,lh);

      CompositeQuadratureRule<SD> * cquad = new CompositeQuadratureRule<SD>() ;

      ELEMENT_TYPE et_time = spacetime ? ET_SEGM : ET_POINT;

      XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                            *cquad, lh, 
                                                                            2*order_space, 2*order_time, 
                                                                            ref_lvl_space, ref_lvl_time);
      DOMAIN_TYPE dt;
      {
        static Timer timer ("XFESpace::GetFE::MakeQuadRule");
        RegionTimer regq (timer);
        dt = xgeom->MakeQuadRule();
      }
      XFiniteElement * retfel = NULL;

      if (spacetime)
      {
        bool also_future_trace = true;

        ScalarFieldEvaluator * lset_eval_past_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,ti.first,lh);
        CompositeQuadratureRule<D> * cquadp = new CompositeQuadratureRule<D>() ;
        XLocalGeometryInformation * xgeom_past = 
          XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_past_p, 
                                            *cquadp, lh, 
                                            2*order_space, 2*order_time, 
                                            ref_lvl_space, ref_lvl_time);
        {
          static Timer timer ("XFESpace::GetFE::PastMakeQuadRule");
          RegionTimer regq (timer);
          xgeom_past->MakeQuadRule();
        }

        if (also_future_trace)
        {
          XLocalGeometryInformation * xgeom_future = NULL;
          CompositeQuadratureRule<D> * cquadf = NULL;
          ScalarFieldEvaluator * lset_eval_future_p = NULL;
          lset_eval_future_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,ti.second,lh);
          cquadf = new CompositeQuadratureRule<D>() ;
          xgeom_future = XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_future_p, 
                                                           *cquadf, lh, 
                                                           2*order_space, 2*order_time, 
                                                           ref_lvl_space, ref_lvl_time);
          {
            static Timer timer ("XFESpace::GetFE::FutureMakeQuadRule");
            RegionTimer regq (timer);
            xgeom_future->MakeQuadRule();
          }
          retfel = new (lh) XFiniteElement(basefes->GetFE(elnr,lh),domnrs,xgeom,xgeom_past,xgeom_future, lh);

          delete xgeom_future;
          delete cquadf;
        }
        else
          retfel = new (lh) XFiniteElement(basefes->GetFE(elnr,lh),domnrs,xgeom,xgeom_past, lh);

        delete xgeom_past;
        delete cquadp;

      }
      else
      {
        retfel = new (lh) XFiniteElement(basefes->GetFE(elnr,lh),domnrs,xgeom, lh);
      }
      delete xgeom;
      delete cquad;

      if (empty)
          retfel->SetEmpty();
      
      return *retfel;
    }
  }


  template <int D, int SD>
  const FiniteElement & XFESpace<D,SD> :: GetSFE (int selnr, LocalHeap & lh) const
  {
    static Timer timer ("XFESpace::GetSFE");
    RegionTimer reg (timer);

    Ngs_Element ngsel = ma.GetSElement(selnr);
    ELEMENT_TYPE eltype = ngsel.GetType();
    if (!activeselem.Test(selnr))
    {
      DOMAIN_TYPE dt = domofsel[selnr];
      return *(new (lh) XDummyFE(dt,eltype));
    }
    else
    {
      Array<DOMAIN_TYPE> domnrs;
      GetSurfaceDomainNrs(selnr,domnrs);  

      Ngs_Element ngel = ma.GetSElement(selnr);
      ELEMENT_TYPE eltype = ngel.GetType();

      ElementTransformation & eltrans = ma.GetTrafo (selnr, BND, lh);
        
      ScalarFieldEvaluator * lset_eval_p = NULL;
      if (spacetime)
        lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,ti,lh);
      else
        lset_eval_p = ScalarFieldEvaluator::Create(D,*coef_lset,eltrans,lh);

      CompositeQuadratureRule<SD-1> * cquad = new CompositeQuadratureRule<SD-1>() ;

      ELEMENT_TYPE et_time = spacetime ? ET_SEGM : ET_POINT;

      XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                            *cquad, lh, 
                                                                            2*order_space, 2*order_time, 
                                                                            ref_lvl_space, ref_lvl_time);
      XFiniteElement * retfel = NULL;

      {
        static Timer timer ("XFESpace::GetSFE::PastMakeQuadRule");
        RegionTimer regq (timer);
        xgeom->MakeQuadRule();
      }

      FlatXLocalGeometryInformation fxgeom(*xgeom, lh);

      retfel = new (lh) XFiniteElement(basefes->GetSFE(selnr,lh),domnrs,xgeom,lh);

      delete xgeom;
      delete cquad;

      if (empty)
          retfel->SetEmpty();

      return *retfel;
    }
  }

  template <int D, int SD>
  void XFESpace<D,SD>::XToNegPos(const GridFunction & gf, GridFunction & gf_neg_pos) const
  {
    GridFunction & gf_neg = *gf_neg_pos.GetComponent(0);
    BaseVector & bv_neg = gf_neg.GetVector();
    FlatVector<> vneg = bv_neg.FVDouble();

    GridFunction & gf_pos = *gf_neg_pos.GetComponent(1);
    BaseVector & bv_pos = gf_pos.GetVector();
    FlatVector<> vpos = bv_pos.FVDouble();

    GridFunction & gf_base = *gf.GetComponent(0);
    BaseVector & bv_base = gf_base.GetVector();
    FlatVector<> vbase = bv_base.FVDouble();

    GridFunction & gf_x = *gf.GetComponent(1);
    BaseVector & bv_x = gf_x.GetVector();
    FlatVector<> vx = bv_x.FVDouble();

    const int basendof = vneg.Size();
    for (int i = 0; i < basendof; ++i)
    {
      vneg(i) = vbase(i);
      vpos(i) = vbase(i);
      const int xdof = basedof2xdof[i];
      if (xdof != -1)
      {
        if (domofdof[xdof] == POS)
          vpos(i) += vx(xdof);
        else
          vneg(i) += vx(xdof);
      }
    }
  }

  template class XFESpace<2,2>;
  template class XFESpace<2,3>;
  template class XFESpace<3,3>;
  template class XFESpace<3,4>;


  LevelsetContainerFESpace::LevelsetContainerFESpace(const MeshAccess & ama, const Flags & flags)
    : FESpace(ama,flags)
  {
    ;
  }


  NumProcInformXFESpace::NumProcInformXFESpace (PDE & apde, const Flags & flags)
    : NumProc (apde)
  { 
    
    FESpace* xh1fes = pde.GetFESpace(flags.GetStringFlag("xh1fespace","v"), true);
    FESpace* xfes = NULL;
    FESpace* basefes = NULL;

    if (xh1fes)
    {
      basefes = (*(dynamic_cast<CompoundFESpace*>(xh1fes)))[0];
      xfes = (*(dynamic_cast<CompoundFESpace*>(xh1fes)))[1];
    }
    else
    {
      xfes = pde.GetFESpace(flags.GetStringFlag("xfespace","vx"));
      basefes = pde.GetFESpace(flags.GetStringFlag("fespace","v"));
    }

    FESpace* fescl = pde.GetFESpace(flags.GetStringFlag("lsetcontfespace","vlc"),true);
    CoefficientFunction * coef_lset_in = pde.GetCoefficientFunction(flags.GetStringFlag("coef_levelset","coef_lset"));

    int mD = pde.GetMeshAccess().GetDimension();

    SpaceTimeFESpace * fes_st = dynamic_cast<SpaceTimeFESpace *>(basefes);
    int mSD = fes_st == NULL ? mD : mD + 1;

    if (mD == 2)
    {
      DomainVariableCoefficientFunction<2> * coef_lset_in_2 = dynamic_cast<DomainVariableCoefficientFunction<2> * > (coef_lset_in);
      int numreg = coef_lset_in_2->NumRegions();
      if (numreg == INT_MAX) numreg = 1;
      Array< EvalFunction* > evals;
      evals.SetSize(numreg);
      for (int i = 0; i < numreg; ++i)
      {
        evals[i] = &coef_lset_in_2->GetEvalFunction(i);
      }
      if (mSD == 2)
      {
        dynamic_cast<XFESpace<2,2>* >(xfes) -> SetBaseFESpace (basefes);
        CoefficientFunction * coef_lset = new DomainVariableCoefficientFunction<2>(evals); 
        dynamic_cast<XFESpace<2,2>* >(xfes) -> SetLevelSetCoefficient (coef_lset);
        if (fescl)
          dynamic_cast<LevelsetContainerFESpace* >(fescl) -> SetLevelSetCoefficient (coef_lset);
      }
      else
      {
        dynamic_cast<XFESpace<2,3>* >(xfes) -> SetBaseFESpace (basefes);
        CoefficientFunction * coef_lset = new DomainVariableCoefficientFunction<3>(evals); 
        dynamic_cast<XFESpace<2,3>* >(xfes) -> SetLevelSetCoefficient (coef_lset);
        if (fescl)
          dynamic_cast<LevelsetContainerFESpace* >(fescl) -> SetLevelSetCoefficient (coef_lset);
      }
    }
    else
    {
      DomainVariableCoefficientFunction<3> * coef_lset_in_3 = dynamic_cast<DomainVariableCoefficientFunction<3> * > (coef_lset_in);
      int numreg = coef_lset_in_3->NumRegions();
      if (numreg == INT_MAX) numreg = 1;
      Array< EvalFunction* > evals;
      evals.SetSize(numreg);
      for (int i = 0; i < numreg; ++i)
      {
        evals[i] = &coef_lset_in_3->GetEvalFunction(i);
      }
      if (mSD == 3)
      {
        dynamic_cast<XFESpace<3,3>* >(xfes) -> SetBaseFESpace (basefes);
        CoefficientFunction * coef_lset = new DomainVariableCoefficientFunction<3>(evals); 
        dynamic_cast<XFESpace<3,3>* >(xfes) -> SetLevelSetCoefficient (coef_lset);
        if (fescl)
          dynamic_cast<LevelsetContainerFESpace* >(fescl) -> SetLevelSetCoefficient (coef_lset);
      }
      else
      {
        dynamic_cast<XFESpace<3,4>* >(xfes) -> SetBaseFESpace (basefes);
        CoefficientFunction * coef_lset = new DomainVariableCoefficientFunction<4>(evals); 
        dynamic_cast<XFESpace<3,4>* >(xfes) -> SetLevelSetCoefficient (coef_lset);
        if (fescl)
          dynamic_cast<LevelsetContainerFESpace* >(fescl) -> SetLevelSetCoefficient (coef_lset);
      }
    }
    // if (xfes_ != NULL)
    // {
    //   xfes_->SetBaseFESpace(basefes);
    // }
    // else
    // {
    //   CompoundFESpace* cfes = dynamic_cast<CompoundFESpace*>(xfes);
    //   if (cfes != NULL)
    //   {
    //     for (int i=0; i < cfes->GetNSpaces(); ++i)
    //     {
    //       XFESpace* xfes_2 = dynamic_cast<XFESpace*>((*cfes)[i]);
    //       if (xfes_2 != NULL)
    //       {
    //         xfes_2->SetBaseFESpace(basefes);
    //       }
    //     }
    //   }
    //   else
    //     throw Exception("FESpace is not compatible to x-fespaces!");
    // }
  }

  NumProcInformXFESpace::~NumProcInformXFESpace(){ ; }
  string NumProcInformXFESpace::GetClassName () const {return "InformXFESpace";  }
  void NumProcInformXFESpace::Do (LocalHeap & lh)  { ; }
  
  static RegisterNumProc<NumProcInformXFESpace> npinfoxfe("informxfem");



  XH1FESpace::XH1FESpace (const MeshAccess & ama, 		   
                          const Array<FESpace*> & aspaces,
                          const Flags & flags)
    : CompoundFESpace(ama, aspaces, flags)
  {
    name="XH1FESpace";

    const int sD = ma.GetDimension();
    if (sD == 2)
      if (dynamic_cast<const XFESpace<2,3> * >(spaces[1]) != NULL)
        spacetime = true;
      else
        spacetime = false;
    else
      if (dynamic_cast<const XFESpace<3,4> * >(spaces[1]) != NULL)
        spacetime = true;
      else
        spacetime = false;


    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
    {
      integrator = new XVisIntegrator<2> (&one);
      // boundary_integrator = new RobinIntegrator<2> (&one);
    }
    else
    {
      integrator = new XVisIntegrator<3> (&one);
      // evaluator = new T_DifferentialOperator<DiffOpVecIdHDG<3> >();
      // boundary_integrator = new RobinVecHDGIntegrator<3> (&one);
    }
  }


  ///SmoothingBlocks for good preconditioned iterative solvers
  Table<int> * XH1FESpace::CreateSmoothingBlocks (const Flags & precflags) const
  {
    Table<int> * it;

    int nv = ma.GetNV();
    int ne = ma.GetNE();
    int nf = ma.GetNFaces();
    int ned = ma.GetNEdges();      
    Array<int> dnums, dnums2, dnums3; //, verts, edges;

    bool stdblock=false;   //put all std. dofs in one block + one block for each xfem dof

    bool blocksystem=true;   //put all std. dofs in one block and all xfem dofs in another block (two large blocks)
    bool vertexpatch=false;   //put all dofs of elements around a vertex together
    bool elementpatch=false;  //put all dofs of one element together
    bool edgepatch=false;     //put all dofs of neighbouring elements together
    bool vertexdofs=false;    //put all dofs of a vertex together
    bool interfacepatch= false;//put all dofs located at (all) cut elements together (one large block)
    bool jacobi= false;        //each degree of freedom is one block

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

    if ( ma.GetDimension() == 2)
    {
      if (blocksystem || stdblock)
      {
        for ( ; !creator.Done(); creator++)
        {
          int basendof = spaces[0]->GetNDof();
          int xndof = spaces[1]->GetNDof();
          for (int i = 0; i < basendof ; ++i)
          {
            creator.Add(0,i);
          }

          if (blocksystem)
          {
            for (int i = 0; i < xndof ; ++i)
            {
              creator.Add(1,i+basendof);
            }
          }
          else
          {
            for (int i = 0; i < xndof ; ++i)
            {
              creator.Add(i+1,i+basendof);
            }
          }
        }
        it = creator.GetTable();
      }
      else
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

          if (vertexpatch || vertexdofs)
          {
            for (int i = 0; i < nv; i++)
            {
              GetVertexDofNrs(i,dnums);
              for (int j = 0; j < dnums.Size(); ++j)
                if (filter.Test(dnums[j]))
                  if (dnums[j] < basendof)
                    creator.Add(dnums[j], dnums);
            }
          }

          if (vertexpatch)
          {
            for (int i = 0; i < ned; i++)
            {
              Ng_Node<1> edge = ma.GetNode<1> (i);

              GetVertexDofNrs(edge.vertices[0],dnums);
              GetVertexDofNrs(edge.vertices[1],dnums2);
              GetEdgeDofNrs(i,dnums3);
            
              for (int j = 0; j < dnums.Size(); ++j)
                for (int k = 0; k < dnums2.Size(); ++k)
                {
                  if (filter.Test(dnums[j]))
                    if (dnums[j] < basendof)
                      creator.Add(dnums[j],dnums2[k]);
                  if (filter.Test(dnums2[k]))
                    if (dnums2[k] < basendof)
                      creator.Add(dnums2[k],dnums[j]);
                }

            
              for (int j = 0; j < dnums.Size(); ++j)
                if (filter.Test(dnums[j]))
                  if (dnums[j] < basendof)
                    creator.Add(dnums[j],dnums3);

              for (int j = 0; j < dnums2.Size(); ++j)
                if (filter.Test(dnums2[j]))
                  if (dnums2[j] < basendof)
                    creator.Add(dnums2[j],dnums3);
            }
          }

          if (vertexdofs || vertexpatch)
            offset += nv;

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
            int nf = ma.GetNFacets();
            Array<int> elnums;
            Array<int> fnums;
            Array<int> ednums;
            Array<int> vnums;
            const XFESpace<2,2>* xfes22 = NULL;
            const XFESpace<2,3>* xfes23 = NULL;
            const XFESpace<3,3>* xfes33 = NULL;
            const XFESpace<3,4>* xfes34 = NULL;
            const int sD = ma.GetDimension();
            if (sD == 2)
              if (spacetime)
                xfes23 = dynamic_cast<const XFESpace<2,3> * >(spaces[1]);
              else
                xfes22 = dynamic_cast<const XFESpace<2,2> * >(spaces[1]);
            else
              if (spacetime)
                xfes34 = dynamic_cast<const XFESpace<3,4> * >(spaces[1]);
              else
                xfes33 = dynamic_cast<const XFESpace<3,3> * >(spaces[1]);
            // Array<int> elnums;
            for (int i = 0; i < nf; i++)
            {
              ma.GetFacetElements(i,elnums);
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

              Ng_Node<1> edge = ma.GetNode<1> (i);
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
                ma.GetElEdges(elnr,ednums);
                for (int l = 0; l < ednums.Size(); ++l)
                  if (ednums[l] != nf)
                  {
                    GetEdgeDofNrs(ednums[l],dnums);
                    creator.Add(offset, dnums);
                  }

                ma.GetElVertices(elnr,vnums);
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
              const XFESpace<2,2>* xfes22 = NULL;
              const XFESpace<2,3>* xfes23 = NULL;
              const XFESpace<3,3>* xfes33 = NULL;
              const XFESpace<3,4>* xfes34 = NULL;
              const int sD = ma.GetDimension();
              if (sD == 2)
                if (spacetime)
                  xfes23 = dynamic_cast<const XFESpace<2,3> * >(spaces[1]);
                else
                  xfes22 = dynamic_cast<const XFESpace<2,2> * >(spaces[1]);
              else
                if (spacetime)
                  xfes34 = dynamic_cast<const XFESpace<3,4> * >(spaces[1]);
                else
                  xfes33 = dynamic_cast<const XFESpace<3,3> * >(spaces[1]);
              bool elcut = sD ? 
                (spacetime ? xfes23->IsElementCut(elnr) 
                 : xfes22->IsElementCut(elnr) )
                :(spacetime ? xfes34->IsElementCut(elnr) 
                  : xfes33->IsElementCut(elnr) );
            
              GetDofNrs(elnr,dnums);
              if (elcut)
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
    return it;
  }


  Array<int> * XH1FESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    bool bddc = true;

    if (true || flags.GetDefineFlag("subassembled"))
    {
      cout << "creating bddc-coarse grid(vertices)" << endl;
      Array<int> & clusters = *new Array<int> (GetNDof());
      clusters = 0;
      // int nv = ma.GetNV();
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
        int nfa = ma.GetNFacets();
        int nv = ma.GetNV();

        for (int i = 0; i < nv; i++)
        {
          // basefes->GetVertexDofNrs (i, dnums);
          clusters[nv /*dnums[0]*/] = 1;
        }

        const BitArray & freedofs = *GetFreeDofs();
        for (int i = 0; i < freedofs.Size(); i++)
          if (!freedofs.Test(i)) clusters[i] = 0;
        *testout << "XH1FESpace, dsc = " << endl << clusters << endl;
        return &clusters;
      }
      else // put all degrees of freedoms at the interface on the coarse grid
      {
        Array<int> dnums;
        BitArray mark(GetNDof());
        mark.Clear();
        int epcnt = 0;
        for (int elnr = 0; elnr < ma.GetNE(); ++elnr)
        {
          const XFESpace<2,2>* xfes22 = NULL;
          const XFESpace<2,3>* xfes23 = NULL;
          const XFESpace<3,3>* xfes33 = NULL;
          const XFESpace<3,4>* xfes34 = NULL;
          const int sD = ma.GetDimension();
          if (sD == 2)
            if (spacetime)
              xfes23 = dynamic_cast<const XFESpace<2,3> * >(spaces[1]);
            else
              xfes22 = dynamic_cast<const XFESpace<2,2> * >(spaces[1]);
          else
            if (spacetime)
              xfes34 = dynamic_cast<const XFESpace<3,4> * >(spaces[1]);
            else
              xfes33 = dynamic_cast<const XFESpace<3,3> * >(spaces[1]);
          bool elcut = sD ? 
            (spacetime ? xfes23->IsElementCut(elnr) 
             : xfes22->IsElementCut(elnr) )
            :(spacetime ? xfes34->IsElementCut(elnr) 
              : xfes33->IsElementCut(elnr) );
            
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
        *testout << "XH1FESpace, dsc = " << endl << clusters << endl;
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
      GetFESpaceClasses().AddFESpace ("xfespace", XFESpace<2,2>::Create);
      GetFESpaceClasses().AddFESpace ("lsetcontfespace", LevelsetContainerFESpace::Create);
      GetFESpaceClasses().AddFESpace ("xh1fespace", XH1FESpace::Create);
    }
  
    Init init;
  }
}
