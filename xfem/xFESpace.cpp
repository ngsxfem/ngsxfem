
#include "xFESpace.hpp"
using namespace ngsolve;
using namespace ngfem;

namespace ngcomp
{
  
  template <int D, int SD>
  XFESpace<D,SD> :: XFESpace (const MeshAccess & ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of XFESpace" << endl;

    throw Exception(" wooohooo - not a good constructor yet ");
    basefes = NULL;
       
    static ConstantCoefficientFunction one(1);
    integrator = new MassIntegrator<D> (&one);
  }
    
  
  template <int D, int SD>
  XFESpace<D,SD> :: ~XFESpace ()
  {
    if (el2dofs) delete el2dofs; 
    if (sel2dofs) delete sel2dofs; 
    // nothing to do
  }
  
  template <int D, int SD>
  void XFESpace<D,SD> :: Update(LocalHeap & lh)
  {
    FESpace<D,SD>::Update(lh);

    int ne=ma.GetNE();
    int nedges=ma.GetNEdges();
    int nv=ma.GetNV();
    int nse=ma.GetNSE();

    activeelem.SetSize(ne);    
    activeselem.SetSize(nse);    
    activeelem.Clear();
    activeselem.Clear();

    BitArray activedofs(basefes->GetNDof());
    activedofs.Clear();

    TableCreator<int> creator;
    for ( ; !creator.Done(); creator++)
    {
      for (int i = 0; i < ne; ++i)
      {

        const FESpace & fes = gf_lset->GetFESpace();
        const FiniteElement & fel = fes.GetFE(elnr, lh);

        // if (fel.ElementType() / 10 + 1 == 2)
        // { 
        //   ScalarFEEvaluator<2> lset_eval(fel, linvec, lh);


// ----------- V ----------- CONTINUE HERE ----------- V -----------
        const ScalarSpaceTimeFiniteElement<2> *  fel_2d_st 
          = dynamic_cast<const ScalarSpaceTimeFiniteElement<2> * >(&fel);


        if (gci->IsElementCut(i))
        {
          activeelem.Set(i);
          Array<int> basednums;
          basefes->GetDofNrs(i,basednums);
          for (int k = 0; k < basednums.Size(); ++k)
          {
            activedofs.Set(basednums[k]);
            creator.Add(i,basednums[k]);
          }
        }
      }
    }
    el2dofs = creator.GetTable();

    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
    {
      for (int i = 0; i < nse; ++i)
      {
        Array<int> fnums;
        ma.GetSElFacets(i,fnums);
        int ed = fnums[0];
        if (gci->IsFacetCut(ed))
        {
          activeselem.Set(i);
          Array<int> basednums;
          basefes->GetSDofNrs(i,basednums);
          for (int k = 0; k < basednums.Size(); ++k)
          {
            activedofs.Set(basednums[k]); // might be twice, but who cares..
            creator2.Add(i,basednums[k]);
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
   
    // domain of dof
    domofdof.SetSize(ndof);


    int domain = 0;
    Array<int> vdofs;
    for (int i = 0; i < nv; i++)
    {
      basefes->GetVertexDofNrs(i,vdofs);
      domain = gci->DomainOfVertex(i);
      for (int j = 0; j < vdofs.Size(); j++)
      {
        if (activedofs.Test(vdofs[j]))
          domofdof[basedof2xdof[vdofs[j]]] = domain;
      }
    }
    Array<int> edofs;
    for (int i = 0; i < nedges; i++)
    {
      basefes->GetEdgeDofNrs(i,edofs);
      domain = gci->DomainOfEdge(i);
      for (int j = 0; j < edofs.Size(); j++)
      {
        if (activedofs.Test(edofs[j]))
          domofdof[basedof2xdof[edofs[j]]] = domain;
      }
    }


    Array<int> indofs;
    for (int i = 0; i < ne; i++)
    {
      basefes->GetInnerDofNrs(i,indofs);
      domain = gci->DomainOfElement(i);
      for (int j = 0; j < indofs.Size(); j++)
      {
        if (activedofs.Test(indofs[j]))
          domofdof[basedof2xdof[indofs[j]]] = domain;
      }
    }


    dirichlet_dofs.SetSize (GetNDof());
    dirichlet_dofs.Clear();
    Array<int> dnums;
    for (int i = 0; i < basedof2xdof.Size(); ++i)
    {
      const int dof = basedof2xdof[i];
      if (dof != -1 && basefes->IsDirichletDof(i))
        dirichlet_dofs.Set (dof);
    }
    
    UpdateCouplingDofArray();
    FinalizeUpdate (lh);

    *testout << "ndof = " << ndof << endl;
    *testout << "el2dofs = " << *el2dofs << endl;
    *testout << "domain of dofs = " << domofdof << endl;
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
  void XFESpace<D,SD> :: GetDomainNrs (int elnr, Array<int> & domnums) const
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
  void XFESpace<D,SD> :: GetSurfaceDomainNrs (int selnr, Array<int> & domnums) const
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
    Ng_Element ngel = ma.GetElement(elnr);
    ELEMENT_TYPE eltype = ConvertElementType(ngel.GetType());
    if (!activeelem.Test(elnr))
    {
      const int sign = gci->DomainOfElement(elnr);
      switch (eltype)
      {
        case ET_SEGM:    
          return *(new (lh) XDummyFE<1>(eltype, sign)); break;
        case ET_TRIG:    
        case ET_QUAD:    
          return *(new (lh) XDummyFE<2>(eltype,sign)); break;
        case ET_TET:
        case ET_PYRAMID:
        case ET_PRISM:
        case ET_HEX:
          return *(new (lh) XDummyFE<3>(eltype,sign)); break;
        default:
        {
          throw Exception ("GetFE not supported for element");
        }
      }

    }
    else
    {
      Array<int> domnrs;
      GetDomainNrs(elnr,domnrs);  
      const MasterElement& masterelement = *(gci->GetCutElement(elnr));
      
      const ScalarFiniteElement<1> * basefel1 = NULL;
      const ScalarFiniteElement<2> * basefel2 = NULL;
      const ScalarFiniteElement<3> * basefel3 = NULL;
      switch (eltype)
      {
        case ET_SEGM:    
          basefel1 = dynamic_cast<const ScalarFiniteElement<1>*>(&(basefes->GetFE(elnr,lh)));
          return *(new (lh) XFiniteElement<1>(*basefel1,domnrs,masterelement)); break;
        case ET_TRIG:    
        case ET_QUAD:    
          basefel2 = dynamic_cast<const ScalarFiniteElement<2>*>(&(basefes->GetFE(elnr,lh)));
          return *(new (lh) XFiniteElement<2>(*basefel2,domnrs,masterelement)); break;
        case ET_TET:     
        case ET_PYRAMID: 
        case ET_PRISM:   
        case ET_HEX:     
          basefel3 = dynamic_cast<const ScalarFiniteElement<3>*>(&(basefes->GetFE(elnr,lh)));
          return *(new (lh) XFiniteElement<3>(*basefel3,domnrs,masterelement)); break;
        default:
        {
          throw Exception ("GetFE not supported for element");
        }
      }
    }
  }

  /// the same for the surface elements
  template <int D, int SD>
  const FiniteElement & XFESpace<D,SD> :: GetSFE (int selnr, LocalHeap & lh) const
  {
    Array<int> fnums;
    ma.GetSElFacets(selnr,fnums);
    int ed = fnums[0];
 
    ELEMENT_TYPE eltype = ma.GetSElType(selnr);
    if (!activeselem.Test(selnr))
    {
      const int sign = gci->DomainOfEdge(ed);
      switch (eltype)
      {
        case ET_SEGM:   
          return *(new (lh) XDummyFE<1>(eltype, sign )); break;
        case ET_TRIG:    
        case ET_QUAD:    
          return *(new (lh) XDummyFE<2>(eltype, sign )); break;
        default:
        {
          throw Exception ("GetSFE not supported for element");
        }
      }

    }
    else
    {
      Array<int> sdomnrs;
      GetSurfaceDomainNrs(selnr,sdomnrs); 
      const Edge & edge = *gci->GetBaseEdge(ed);
      // adapt (see GetFE)
      const ScalarFiniteElement<1> * basefel1 = NULL;
      const ScalarFiniteElement<2> * basefel2 = NULL;
      // const ScalarFiniteElement<3> * basefel3 = NULL;

      switch (eltype)
      {
        case ET_SEGM:    
          basefel1 = dynamic_cast<const ScalarFiniteElement<1>*>(&(basefes->GetSFE(selnr,lh)));
          return *(new (lh) XFiniteElement<1>(*basefel1,sdomnrs,edge)); break;
        case ET_TRIG:
        case ET_QUAD: 
          basefel1 = dynamic_cast<const ScalarFiniteElement<1>*>(&(basefes->GetSFE(selnr,lh)));
          return *(new (lh) XFiniteElement<2>(*basefel2,sdomnrs,edge)); break;
        default:
        {
          throw Exception ("GetSFE not supported for element");
        }
      }
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

  /*
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
      GetFESpaceClasses().AddFESpace ("xh1fespace", XH1FESpace::Create);
    }
  
    Init init;
  }
  */
}
