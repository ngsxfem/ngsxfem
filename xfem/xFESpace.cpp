
#include "xFESpace.hpp"
//#include "xfemVisInts.hpp"
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

    ti.first = flags.GetNumFlag("t0",0.0);
    ti.second = flags.GetNumFlag("t1",1.0);

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
  XFESpace<D,SD> :: ~XFESpace ()
  {
    if (el2dofs) delete el2dofs; 
    if (sel2dofs) delete sel2dofs; 
    if (eval_lset) delete eval_lset;
    // nothing to do
  }
  
  template <int D, int SD>
  void XFESpace<D,SD> :: Update(LocalHeap & lh)
  {
    if ( basefes == NULL )
    {
      cout << " no basefes, Update postponed " << endl;
      return;
    }
    
    order_space = basefes->GetOrder();
    
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

    TableCreator<int> creator;
    for ( ; !creator.Done(); creator++)
    {
      for (int elnr = 0; elnr < ne; ++elnr)
      {
        HeapReset hr(lh);

        netgen::Ng_Element ngel = ma.GetElement(elnr);
        ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());

        ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), lh);
        
        ScalarFieldEvaluator * lset_eval_p = NULL;
        if (spacetime)
          lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,ti,lh);
        else
          lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,lh);

        CompositeQuadratureRule<SD> cquad;
        XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                              cquad, lh, 2*order_space, 1, 0, 0);
        DOMAIN_TYPE dt = xgeom->MakeQuadRule();

        domofel[elnr] = dt;

        if (dt == IF)// IsElementCut ?
        {
          activeelem.Set(elnr);
          Array<int> basednums;
          basefes->GetDofNrs(elnr,basednums);
          for (int k = 0; k < basednums.Size(); ++k)
          {
            activedofs.Set(basednums[k]);
            creator.Add(elnr,basednums[k]);
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

        netgen::Ng_Element ngel = ma.GetSElement(selnr);
        ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());

        ElementTransformation & seltrans = ma.GetTrafo (selnr, BND, lh);

        ScalarFieldEvaluator * lset_eval_p = NULL;
        if (spacetime)
          lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,seltrans,ti,lh);
        else
          lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,seltrans,lh);

        CompositeQuadratureRule<SD-1> cquad;
        XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                              cquad, lh, 2*order_space, 1, 0, 0);
        DOMAIN_TYPE dt = xgeom->MakeQuadRule();

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

    *testout << "basefes -> free_dofs = " << *basefes->GetFreeDofs() << endl;

    *testout << "free_dofs = " << free_dofs << endl;
  }
  

  template <int D, int SD>
  void XFESpace<D,SD> :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    if (activeelem.Test(elnr))
      dnums = (*el2dofs)[elnr];
    else
      dnums.SetSize(0);
  }
  
  template <int D, int SD>
  void XFESpace<D,SD> :: GetDomainNrs (int elnr, Array<DOMAIN_TYPE> & domnums) const
  {
    if (activeelem.Test(elnr))
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
    if (activeselem.Test(selnr))
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

    for (int i = 0; i < basedof2xdof.Size(); ++i)
    {
      const int dof = basedof2xdof[i];
      if (dof != -1)
        ctofdof[dof] = basefes->GetDofCouplingType(i);
    }
    *testout << "XFESpace, ctofdof = " << endl << ctofdof << endl;
  }


  template <int D, int SD>
  void XFESpace<D,SD> :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    if (activeselem.Test(selnr))
      dnums = (*sel2dofs)[selnr];
    else
      dnums.SetSize(0);
  }
  template <int D, int SD>
  const FiniteElement & XFESpace<D,SD> :: GetFE (int elnr, LocalHeap & lh) const
  {
    netgen::Ng_Element ngel = ma.GetElement(elnr);
    ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());
    if (!activeelem.Test(elnr))
    {
      DOMAIN_TYPE dt = domofel[elnr];
      return *(new (lh) XDummyFE(dt,eltype));
    }
    else
    {
      Array<DOMAIN_TYPE> domnrs;
      GetDomainNrs(elnr,domnrs);  

      netgen::Ng_Element ngel = ma.GetElement(elnr);
      ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());

      ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), lh);
        
      ScalarFieldEvaluator * lset_eval_p = NULL;
      if (spacetime)
        lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,ti,lh);
      else
        lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,lh);

      CompositeQuadratureRule<SD> * cquad = new (lh) CompositeQuadratureRule<SD>() ;

      ELEMENT_TYPE et_time = spacetime ? ET_SEGM : ET_POINT;

      XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                            *cquad, lh, 2*order_space, 1, 0, 0);
      DOMAIN_TYPE dt = xgeom->MakeQuadRule();

      if (spacetime)
      {
          ScalarFieldEvaluator * lset_eval_past_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,ti.first,lh);
          ScalarFieldEvaluator * lset_eval_future_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,ti.second,lh);
          CompositeQuadratureRule<D> * cquadp = new (lh) CompositeQuadratureRule<D>() ;
          CompositeQuadratureRule<D> * cquadf = new (lh) CompositeQuadratureRule<D>() ;
          XLocalGeometryInformation * xgeom_past = 
              XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_past_p, 
                                                *cquadp, lh, 2*order_space, 1, 0, 0);
          XLocalGeometryInformation * xgeom_future = 
              XLocalGeometryInformation::Create(eltype, ET_POINT, *lset_eval_future_p, 
                                                *cquadf, lh, 2*order_space, 1, 0, 0);
          xgeom->SetPastTrace(xgeom_past);
          xgeom->SetFutureTrace(xgeom_future);
      }


      return *(new (lh) XFiniteElement(basefes->GetFE(elnr,lh),domnrs,xgeom));
    }
  }


  template <int D, int SD>
  const FiniteElement & XFESpace<D,SD> :: GetSFE (int selnr, LocalHeap & lh) const
  {
    netgen::Ng_Element ngsel = ma.GetSElement(selnr);
    ELEMENT_TYPE eltype = ConvertElementType(ngsel.GetType());
    if (!activeselem.Test(selnr))
    {
      DOMAIN_TYPE dt = domofel[selnr];
      return *(new (lh) XDummyFE(dt,eltype));
    }
    else
    {
      Array<DOMAIN_TYPE> domnrs;
      GetSurfaceDomainNrs(selnr,domnrs);  

      netgen::Ng_Element ngel = ma.GetSElement(selnr);
      ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());

      ElementTransformation & eltrans = ma.GetTrafo (selnr, BND, lh);
        
      ScalarFieldEvaluator * lset_eval_p = NULL;
      if (spacetime)
        lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,ti,lh);
      else
        lset_eval_p = ScalarFieldEvaluator::Create(D,*eval_lset,eltrans,lh);

      CompositeQuadratureRule<SD-1> * cquad = new (lh) CompositeQuadratureRule<SD-1>() ;

      ELEMENT_TYPE et_time = spacetime ? ET_SEGM : ET_POINT;

      XLocalGeometryInformation * xgeom = XLocalGeometryInformation::Create(eltype, et_time, *lset_eval_p, 
                                                                            *cquad, lh, 2*order_space, 1, 0, 0);
      return *(new (lh) XFiniteElement(basefes->GetSFE(selnr,lh),domnrs,xgeom));
    }
  }

  /*
  XH1FESpace :: XH1FESpace (const MeshAccess & ama, 
                            const Array<FESpace*> & aspaces,
                            const Flags & flags)
    : CompoundFESpace (ama, aspaces, flags)
  { 
    static ConstantCoefficientFunction one(1);
    integrator = new SumMassIntegrator<D> (&one);
    boundary_integrator = new RobinIntegrator<D> (&one);
    boundary_integrator = new CompoundBilinearFormIntegrator (*boundary_integrator, 0);
  }

  */


    
  NumProcInformXFESpace::NumProcInformXFESpace (PDE & apde, const Flags & flags)
    : NumProc (apde)
  { 
    FESpace* xfes = pde.GetFESpace(flags.GetStringFlag("xfespace","vx"));
    FESpace* basefes = pde.GetFESpace(flags.GetStringFlag("fespace","v"));
    
    int mD = pde.GetMeshAccess().GetDimension();
    
    SpaceTimeFESpace * fes_st = dynamic_cast<SpaceTimeFESpace *>(basefes);
    int mSD = fes_st == NULL ? mD : mD + 1;

    if (mD == 2)
      if (mSD == 2)
        dynamic_cast<XFESpace<2,2>* >(xfes) -> SetBaseFESpace (basefes);
      else
        dynamic_cast<XFESpace<2,3>* >(xfes) -> SetBaseFESpace (basefes);
    else
      if (mSD == 2)
        dynamic_cast<XFESpace<3,2>* >(xfes) -> SetBaseFESpace (basefes);
      else
        dynamic_cast<XFESpace<3,3>* >(xfes) -> SetBaseFESpace (basefes);

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
      // GetFESpaceClasses().AddFESpace ("xh1fespace", XH1FESpace::Create);
    }
  
    Init init;
  }
}
