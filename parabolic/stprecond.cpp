/*********************************************************************/
/* File:   stprecond.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   June 2014                                                 */
/*********************************************************************/

#include <solve.hpp>
#include "stprecond.hpp"

extern ngsolve::PDE * pde;
namespace ngcomp
{


  SpaceTimePreconditioner :: SpaceTimePreconditioner (const PDE & pde, const Flags & flags, const string & name)  
    : Preconditioner (&pde, flags)
  {
    bfa = dynamic_cast<const S_BilinearForm<double>*>
      (pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL)));
  }
  

  SpaceTimePreconditioner :: SpaceTimePreconditioner (const BaseMatrix & matrix, const BitArray * afreedofs)
    : Preconditioner (pde, Flags()), freedofs(afreedofs)
  {
    Setup (matrix);
  }
  

  SpaceTimePreconditioner :: ~SpaceTimePreconditioner ()
  {
    delete AssBlock;
    delete InvAssBlock;
    // delete blockjacobixtable;
    delete blockjacobix;
    ;
  }
  
  
  
  void SpaceTimePreconditioner :: Update()
  {
    freedofs = bfa->GetFESpace().GetFreeDofs (bfa->UsesEliminateInternal());
    Setup (bfa->GetMatrix());
  }
  

  void SpaceTimePreconditioner :: Setup (const BaseMatrix & matrix)
  {
    cout << IM(1) << "Setup SpaceTime preconditioner" << endl;
    static Timer t("SpaceTime setup");
    RegionTimer reg(t);

    const SparseMatrix<double> & mat = dynamic_cast< const SparseMatrix<double> &>(matrix);

    if (dynamic_cast< const SparseMatrixSymmetric<double> *> (&mat))
      throw Exception ("Please use fully stored sparse matrix for SpaceTime (bf -nonsymmetric)");


    const FESpace & fesh1x = bfa->GetFESpace();
    const FESpace & fesh1 = *((dynamic_cast<const CompoundFESpace &>(fesh1x))[0]);
    const FESpace & fesx = *((dynamic_cast<const CompoundFESpace &>(fesh1x))[1]);

    const MeshAccess & ma = fesh1.GetMeshAccess();

    int ndof_h1 = fesh1.GetNDof();
    int ndof_x = fesh1x.GetNDof() - fesh1.GetNDof();
    int ndof = fesh1x.GetNDof();

    Array<int> dnums;

    TableCreator<int> creator(ma.GetNE());
    for ( ; !creator.Done(); creator++)
    {    
      for (ElementId ei : ma.Elements(VOL))
      {	
        // if (!DefinedOn (ei)) continue;
        int i = ei.Nr();
        ELEMENT_TYPE eltype = ma.GetElType(ei); 
      
        fesh1.GetDofNrs(i,dnums);
        for (int j = 0; j < dnums.Size(); ++j)
          creator.Add(i,dnums[j]);
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
          ai=i;
          aj=cols[j];
          aval=values[j];
          // std::cout << " ai = " << ai << std::endl;
          // std::cout << " aj = " << aj << std::endl;
          // std::cout << " aval = " << aval << std::endl;
          AssBlock->AddElementMatrix(ai,aj,aval);
	      // row = i
          // col = cols[j]
          // entry = values[j]
	    }
    }


    TableCreator<int> creator2(ndof_x/2);
    for ( ; !creator2.Done(); creator2++)
    {    
      
      int nv = ma.GetNV();
      int cnt = 0;
      for (int i = 0; i < nv; i++)
      {
        fesx.GetVertexDofNrs(i,dnums);
        // std::cout << " dnums = " << dnums << std::endl;
        if (dnums.Size()==0) continue;
        for (int j = 0; j < dnums.Size(); ++j)
          creator2.Add(cnt, ndof_h1 + dnums[j]);
        cnt++;
      }

      // for (int i = 0; i < ndof_x; ++i)
      // {	
      //   creator2.Add(i, ndof_h1+i);
      // }
    }

    delete blockjacobix; //also deletes the table
    blockjacobixtable= creator2.GetTable();
    // std::cout << " *blockjacobixtable = " << *blockjacobixtable << std::endl;
    blockjacobix = new BlockJacobiPrecond<double> (mat, *blockjacobixtable);

    delete InvAssBlock;
    dynamic_cast<BaseSparseMatrix&> (*AssBlock) . SetInverseType (inversetype);
    InvAssBlock = AssBlock->InverseMatrix(fesh1.GetFreeDofs());

    delete & element2dof;
  }


  void SpaceTimePreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("SpaceTime mult");
    RegionTimer reg(t);

    const FESpace & fesh1x = bfa->GetFESpace();
    const FESpace & fesh1 = *((dynamic_cast<const CompoundFESpace &>(fesh1x))[0]);

    int ndof_h1 = fesh1.GetNDof();

    const FlatVector<double> fvf =f.FVDouble();
    FlatVector<double> fu =u.FVDouble();
    
    VFlatVector<double> fh1(ndof_h1,&fvf(0));
    VFlatVector<double> uh1(ndof_h1,&fu(0));
    
    // std::cout << " fh1 = " << fh1 << std::endl;
    // std::cout << " uh1 = " << uh1 << std::endl;
    // getchar();
    u = 0.0;

    // u = f;
    u = *blockjacobix * f;
    // uh1 = 0.0;
    // std::cout << " f = " << f << std::endl;
    // std::cout << " u = " << u << std::endl;
    // getchar();
    uh1 += *InvAssBlock * fh1;
    // // std::cout << " f = " << f << std::endl;
    // std::cout << " u = " << u << std::endl;
    // getchar();
    
    u=f;

  }

  static RegisterPreconditioner<SpaceTimePreconditioner> init_SpaceTimepre ("spacetime");
}
