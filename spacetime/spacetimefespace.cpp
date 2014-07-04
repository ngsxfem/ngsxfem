
#include "spacetimefespace.hpp"
#include "spacetimefe.hpp"
#include "spacetimeintegrators.hpp"

// #include "../fem/spacetimefe.hpp"

namespace ngcomp
{

  SpaceTimeFESpace ::  SpaceTimeFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags)
    : FESpace(ama, flags)
  {
    name="SpaceTimeFESpace";
    // defined flags
    DefineNumFlag("order_time");
    DefineNumFlag("order_space");
    DefineStringFlag("type_space");
    DefineDefineFlag("gaussradau");
    // DefineDefineFlag("variableorder"); 

    if(checkflags) CheckFlags(flags);
    
    // ndlevel.SetSize(0);
    spacedim = ma.GetDimension();
    
    order_time = flags.GetNumFlag("order_time",1);
    gaussradau = flags.GetDefineFlag("gaussradau");
    order_space = flags.GetNumFlag("order_space",1);
    if (gaussradau)
        fel_time = new GR2PFiniteElement();
    else
        fel_time = new L2HighOrderFE<ET_SEGM> (order_time);

    static ConstantCoefficientFunction one(1);
    if (spacedim == 2)
    {
      integrator = new SpaceTimeTimeTraceIntegrator<2,FUTURE>(&one) ;
      boundary_integrator = new RobinIntegrator<2> (&one);
    }
    else
    {
      integrator = new SpaceTimeTimeTraceIntegrator<3,FUTURE>(&one) ;
      // integrator = new MassIntegrator<3> (&one);
      boundary_integrator = new RobinIntegrator<3> (&one);
    }

    if (dimension > 1)
    {
      integrator = new BlockBilinearFormIntegrator (*integrator, dimension);
      boundary_integrator =
        new BlockBilinearFormIntegrator (*boundary_integrator, dimension);
    }
    
    const FESpaceClasses::FESpaceInfo * info;

    string fet_space = flags.GetStringFlag("type_space","l2ho");

    info = GetFESpaceClasses().GetFESpace(fet_space);
    if (!info) throw Exception("SpaceTimeFESpace ::  SpaceTimeFESpace : fespace not given ");
    
    Flags fespaceflags(flags);
    fespaceflags.SetFlag("order",order_space);
    fes_space = info->creator(ma, fespaceflags);

  }

  

  SpaceTimeFESpace :: ~SpaceTimeFESpace ()
  {
    delete fel_time; 
  }




  void SpaceTimeFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);

    fes_space->Update(lh);

    if(print) 
      *testout << " SpaceTimeFESpace" << endl; // with order " << order << " rel_order " << rel_order << " var_order " << var_order << endl; 

    nel = ma.GetNE();
    nfa = ma.GetNFacets(); 
    
    ndof_space = fes_space->GetNDof();
    ndof_time  =  fel_time->GetNDof();
    ndof = ndof_space * ndof_time;
    
    // while (ma.GetNLevels() > ndlevel.Size())
    //   ndlevel.Append (ndof);
    // ndlevel.Last() = ndof;
      
    if(print)
    {
      *testout << "*** Update SpaceTImeFESpace: General Information" << endl;
      *testout << " order_time: " << order_time << endl; 
      *testout << " order_space: " << order_space << endl; 
    } 

    UpdateCouplingDofArray();

    if (timing) Timing();
  }


  void SpaceTimeFESpace :: UpdateCouplingDofArray()
  {
    ctofdof.SetSize(ndof);
    for (int i = 0; i < ndof_space; ++i)
    {
      COUPLING_TYPE ct = fes_space->GetDofCouplingType(i);
      for (int j = 0; j < ndof_time; ++j)
        ctofdof[j*ndof_space+i] = ct;
      
    }
    // take space-fespaces coupling dofs for all time levels
    if (print)
      *testout << "SpaceTimeFESpace, ctofdof = " << endl << ctofdof << endl;

    // throw Exception("SpaceTimeFESpace :: UpdateCouplingDofArray - not implemented yet");
  }


// ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetFE (int elnr, LocalHeap & lh) const
  { 

    const DGFiniteElement<1> * fel_time_radau = new (lh) GR2PFiniteElement();
    const DGFiniteElement<1> * fel_time_dg = new (lh) L2HighOrderFE<ET_SEGM> (order_time);
    const DGFiniteElement<1> * fel_time = gaussradau ? fel_time_radau : fel_time_dg;
    
    const FiniteElement * fel_space = &(fes_space->GetFE(elnr,lh));   
    if ( spacedim == 2 ) 
    {
      const ScalarFiniteElement<2> * fel_scalar = dynamic_cast<const ScalarFiniteElement<2> *>(fel_space);    
      if (fel_scalar == NULL)
        throw Exception("SpaceTimeFESpace :: GetFE Cast to ScalarFiniteElement failed!");
      return *(new (lh) ScalarSpaceTimeFiniteElement<2>(*fel_scalar,*fel_time));
    }
    else
    {
      const ScalarFiniteElement<3> * fel_scalar = dynamic_cast<const ScalarFiniteElement<3> *>(fel_space);    
      if (fel_scalar == NULL)
        throw Exception("SpaceTimeFESpace :: GetFE Cast to ScalarFiniteElement failed!");
      return *(new (lh) ScalarSpaceTimeFiniteElement<3>(*fel_scalar,*fel_time));
    }
  }

// ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetTimeFE (LocalHeap & lh) const
  {
    if (gaussradau)
      return *(new (lh) GR2PFiniteElement() );
    else
      return *(new (lh) L2HighOrderFE<ET_SEGM> (order_time));
  }

// ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetSpaceFE (int elnr, LocalHeap & lh) const
  {
    return fes_space->GetFE( elnr, lh);   
  }

// ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    const DGFiniteElement<1> * fel_time_radau = new (lh) GR2PFiniteElement();
    const DGFiniteElement<1> * fel_time_dg = new (lh) L2HighOrderFE<ET_SEGM> (order_time);
    const DGFiniteElement<1> * fel_time = gaussradau ? fel_time_radau : fel_time_dg;
        
    const FiniteElement * fel_space = &(fes_space->GetSFE(selnr,lh));   
    if ( spacedim == 2 ) 
    {
      const ScalarFiniteElement<1> * fel_scalar = dynamic_cast<const ScalarFiniteElement<1> *>(fel_space);    
      if (fel_scalar == NULL)
        throw Exception("SpaceTimeFESpace :: GetFE Cast to ScalarFiniteElement failed!");
      return *(new (lh) ScalarSpaceTimeFiniteElement<1>(*fel_scalar,*fel_time));
    }
    else
    {
      const ScalarFiniteElement<2> * fel_scalar = dynamic_cast<const ScalarFiniteElement<2> *>(fel_space);    
      if (fel_scalar == NULL)
        throw Exception("SpaceTimeFESpace :: GetFE Cast to ScalarFiniteElement failed!");
      return *(new (lh) ScalarSpaceTimeFiniteElement<2>(*fel_scalar,*fel_time));
    }

    // throw Exception("SpaceTimeFESpace :: GetSFE - not implemented yet");
    // // take space-fespaces GetSFE and tp-extend with fel_time
  }



// ------------------------------------------------------------------------
  int SpaceTimeFESpace :: GetNDof () const
  {
    return ndof;
  }

// ------------------------------------------------------------------------
  int SpaceTimeFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }



// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    Array<int> dnums_space;
    fes_space->GetDofNrs(elnr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }

// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    Array<int> dnums_space;
    fes_space->GetSDofNrs(selnr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }


// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetVertexDofNrs ( int nr, Array<int> & dnums ) const
  {
    Array<int> dnums_space;
    fes_space->GetVertexDofNrs(nr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }

// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetEdgeDofNrs ( int nr, Array<int> & dnums ) const
  {
    Array<int> dnums_space;
    fes_space->GetEdgeDofNrs(nr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }

// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetFaceDofNrs (int nr, Array<int> & dnums) const
  {
    Array<int> dnums_space;
    fes_space->GetFaceDofNrs(nr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }
  
// ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    Array<int> dnums_space;
    fes_space->GetInnerDofNrs(elnr, dnums_space);
    const int ndof_local_space = dnums_space.Size();
    dnums.SetSize(ndof_local_space*ndof_time);

    for (int i = 0; i < ndof_local_space; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        dnums[j*ndof_local_space+i] = dnums_space[i] + j * ndof_space;
      }
    }
  }



// ------------------------------------------------------------------------
  Table<int> * SpaceTimeFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    throw Exception("SpaceTimeFESpace :: CreateDirectSolverClusters - not implemented yet");
    Table<int> & table = *new Table<int> (0,0);
    return &table;

  }


  
  Array<int> * SpaceTimeFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    throw Exception("SpaceTimeFESpace :: CreateDirectSolverClusters - not implemented yet");
    Array<int> & clusters = *new Array<int> (0,0);
    return & clusters;
  }
  
// ------------------------------------------------------------------------

  
  static RegisterFESpace<SpaceTimeFESpace> init_spacetimefes ("spacetimefes");
}


