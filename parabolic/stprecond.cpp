/*********************************************************************/
/* File:   stprecond.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   June 2014                                                 */
/*********************************************************************/

#include <solve.hpp>
#include "stprecond.hpp"
#include "../xfem/xFESpace.hpp"

extern ngsolve::PDE * pde;
namespace ngcomp
{


  SpaceTimePreconditioner :: SpaceTimePreconditioner (const PDE & pde, const Flags & flags, const string & name)  
    : Preconditioner (&pde, flags)
  {
    bfa = dynamic_pointer_cast<S_BilinearForm<double> >
      (pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL)));
  }

  SpaceTimePreconditioner :: SpaceTimePreconditioner (shared_ptr<BilinearForm> bfa_in, const Flags & flags, const string & name)  
    : Preconditioner (bfa_in, flags)
  {
    bfa = dynamic_pointer_cast<S_BilinearForm<double> >(bfa_in);
  }
  

  SpaceTimePreconditioner :: SpaceTimePreconditioner (const BaseMatrix & matrix, const BitArray * afreedofs)
    : Preconditioner (pde, Flags()), freedofs(afreedofs)
  {
    Setup (matrix);
  }
  

  SpaceTimePreconditioner :: ~SpaceTimePreconditioner ()
  {
    delete AssBlock;
    // delete InvAssBlock;
    // delete blockjacobixtable;
    delete blockjacobix;
    ;
  }
  
  
  
  void SpaceTimePreconditioner :: Update()
  {
    freedofs = bfa->GetFESpace()->GetFreeDofs (bfa->UsesEliminateInternal());
    Setup (bfa->GetMatrix());
    if (test) Test();

  }
  

  void SpaceTimePreconditioner :: Setup (const BaseMatrix & matrix)
  {
    cout << IM(1) << "Setup SpaceTime preconditioner" << endl;
    static Timer t("SpaceTime setup");
    RegionTimer reg(t);

    const SparseMatrix<double> & mat = dynamic_cast< const SparseMatrix<double> &>(matrix);

    if (dynamic_cast< const SparseMatrixSymmetric<double> *> (&mat))
      throw Exception ("Please use fully stored sparse matrix for SpaceTime (bf -nonsymmetric)");

    shared_ptr<XStdFESpace> fesh1x = dynamic_pointer_cast<XStdFESpace>(bfa->GetFESpace());
    bool spacetime = fesh1x->IsSpaceTime();
    // LocalHeap lh(6000000,"test");
    // const_cast<FESpace &>(fesh1x).Update(lh); 
    shared_ptr<FESpace> fesh1 = (*fesh1x)[0];
    shared_ptr<FESpace> fesx = (*fesh1x)[1];

    shared_ptr<MeshAccess> ma = fesh1->GetMeshAccess();

    int ndof_h1 = fesh1->GetNDof();
    int ndof_x = fesh1x->GetNDof() - fesh1->GetNDof();
    // int ndof = fesh1x->GetNDof();

    Array<int> dnums;

    bool dgcoupl = fesh1x->UsesDGCoupling();

    TableCreator<int> creator(dgcoupl ? ma->GetNE() + ma->GetNFacets() : ma->GetNE());
    for ( ; !creator.Done(); creator++)
    {    
      for (ElementId ei : ma->Elements(VOL))
      {	
        // if (!DefinedOn (ei)) continue;
        int i = ei.Nr();
        // ELEMENT_TYPE eltype = ma->GetElType(ei); 
      
        fesh1->GetDofNrs(i,dnums);
        for (int j = 0; j < dnums.Size(); ++j)
          creator.Add(i,dnums[j]);
      }


      if (dgcoupl)
      {
        Array<int> nbelems;
        Array<int> elnums;
        //add dofs of neighbour elements as well
        for (int i = 0; i < ma->GetNFacets(); i++)
        {
          nbelems.SetSize(0);
          ma->GetFacetElements(i,elnums);
          for (int k=0; k<elnums.Size(); k++)
            nbelems.Append(elnums[k]);
          
          for (int k=0;k<nbelems.Size();k++){
            int elnr=nbelems[k];
            // if (!fesh1.DefinedOn (ma->GetElIndex(elnr))) continue;
            fesh1->GetDofNrs (elnr, dnums);
            for (int j = 0; j < dnums.Size(); j++)
              if (dnums[j] != -1)
                creator.Add (ma->GetNE()+i, dnums[j]);
          }
        }
      }
    }

    Table<int> & element2dof = *(creator.GetTable());
    delete AssBlock;
    AssBlock = new SparseMatrix<double>(ndof_h1,
                                        element2dof,
                                        element2dof,
                                        false);
    

    Array<int> ai(1);
    Array<int> aj(1);
    Matrix<double> aval(1,1);
    AssBlock->AsVector() = 0.0;
    for( int i = 0; i < ndof_h1; i++)
    {
      int row = i;
      if (row == -1) continue;

      FlatArray<int> cols = mat.GetRowIndices(i);
      FlatVector<double> values = mat.GetRowValues(i);
      
      // Array<int> cols_global;
      // Array<double> values;

      for( int j = 0; j < cols.Size(); j++)
        if (cols[j] != -1 && cols[j] < ndof_h1)
	    {
          ai[0]=i;
          aj[0]=cols[j];
          aval(0,0)=values[j];
          // std::cout << " ai = " << ai << std::endl;
          // std::cout << " aj = " << aj << std::endl;
          // std::cout << " aval = " << aval << std::endl;
          AssBlock->AddElementMatrix(ai,aj,aval);
	      // row = i
          // col = cols[j]
          // entry = values[j]
	    }
    }

    // std::cout << " AssBlock->Print() = ";
    // AssBlock->Print(cout);
    // cout << std::endl;
    // std::cout << " mat.Print() = ";
    // mat.Print(cout);
    // cout << std::endl;
    // getchar();

    const BitArray * freedofs = fesh1x->GetFreeDofs();

    TableCreator<int> creator2( spacetime ? ndof_x/2 : ndof_x);
    for ( ; !creator2.Done(); creator2++)
    {    
      
      int nv = ma->GetNV();
      int cnt = 0;
      for (int i = 0; i < nv; i++)
      {
        fesx->GetVertexDofNrs(i,dnums);
        if (dnums.Size()==0) continue;
        for (int j = 0; j < dnums.Size(); ++j)
            if(freedofs->Test(ndof_h1 + dnums[j]))
                creator2.Add(cnt, ndof_h1 + dnums[j]);
        cnt++;
      }
    }

    delete blockjacobix; //also deletes the table
    blockjacobixtable= creator2.GetTable();
    // std::cout << " *blockjacobixtable = " << *blockjacobixtable << std::endl;
    blockjacobix = new BlockJacobiPrecond<double> (mat, *blockjacobixtable);

    // delete InvAssBlock;
    dynamic_cast<BaseSparseMatrix&> (*AssBlock) . SetInverseType (inversetype);
    InvAssBlock = AssBlock->InverseMatrix(fesh1->GetFreeDofs());
    // std::cout << " *fesh1.GetFreeDofs() = " << *fesh1.GetFreeDofs() << std::endl;

    delete & element2dof;

  }


  void SpaceTimePreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("SpaceTime mult");
    RegionTimer reg(t);

    shared_ptr<FESpace> fesh1x = bfa->GetFESpace();
    shared_ptr<FESpace> fesh1 = (*(dynamic_pointer_cast<CompoundFESpace>(fesh1x)))[0];

    // const BitArray * freedofs = fesh1x->GetFreeDofs();

    int ndof_h1 = fesh1->GetNDof();

    const FlatVector<double> fvf =f.FVDouble();
    FlatVector<double> fu =u.FVDouble();
    
    VFlatVector<double> fh1(ndof_h1,&fvf(0));
    VFlatVector<double> uh1(ndof_h1,&fu(0));
    
    // std::cout << " f = " << f << std::endl;
    // std::cout << " fh1 = " << fh1 << std::endl;
    // std::cout << " uh1 = " << uh1 << std::endl;
    // getchar();
    u = 0.0;

    // u = f;
    // fu = *blockjacobix * fvf;
    blockjacobix->MultAdd(1.0,f,u);

    // for (int i = 0; i < freedofs->Size(); ++i)
    //     if(!freedofs->Test(i))
    //         fu(i) = 0.0;

    uh1 = 0.0;
    // std::cout << " uh1 = " << uh1 << std::endl;
    // std::cout << " u = " << u << std::endl;
    // getchar();
    // std::cout << " InvAssBlock->Height() = " << InvAssBlock->Height() << std::endl;
    // std::cout << " ndof_h1 = " << ndof_h1 << std::endl;
    uh1 = *InvAssBlock * fh1;
    // // std::cout << " f = " << f << std::endl;
    // std::cout << " u = " << u << std::endl;
    // getchar();
    
    //u=f;

  }

  static RegisterPreconditioner<SpaceTimePreconditioner> init_SpaceTimepre ("spacetime");
}
