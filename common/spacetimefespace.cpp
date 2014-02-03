
#include <comp.hpp>
#include <fem.hpp> 

#include "spacetimefespace.hpp"

// #include "../fem/spacetimefe.hpp"

namespace ngcomp
{


  SpaceTimeFESpace ::  SpaceTimeFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags)
    : FESpace(ama, flags)
  {
    name="SpaceTimeFESpace(facet)";
    // defined flags
    // DefineNumFlag("relorder");
    // DefineDefineFlag("variableorder"); 

    // if(checkflags) CheckFlags(flags);
    
    
    // ndlevel.SetSize(0);
    
    // TODO: Evaluator for shape tester 
    static ConstantCoefficientFunction one(1);
    if (ma.GetDimension() == 2)
      {
        integrator = new MassIntegrator<2> (&one);
        boundary_integrator = new RobinIntegrator<2> (&one);
      }
    else
      {
        integrator = new MassIntegrator<3> (&one);
        boundary_integrator = new RobinIntegrator<3> (&one);
      }

    if (dimension > 1)
      {
        integrator = new BlockBilinearFormIntegrator (*integrator, dimension);
        boundary_integrator =
	  new BlockBilinearFormIntegrator (*boundary_integrator, dimension);
      }
  }

  

  SpaceTimeFESpace :: ~SpaceTimeFESpace ()
  { ; }




  void SpaceTimeFESpace :: Update(LocalHeap & lh)
  {
    FESpace :: Update (lh);

    fes_space->Update(lh);

    if(print) 
      *testout << " SpaceTimeFESpace with order " << order << " rel_order " << rel_order << " var_order " << var_order << endl; 

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
	*testout << " order_time: no info yet " <<  << endl; 
      } 

    UpdateCouplingDofArray();

    if (timing) Timing();
  }


  void SpaceTimeFESpace :: UpdateCouplingDofArray()
  {
    throw Exception("SpaceTimeFESpace :: UpdateCouplingDofArray - not implemented yet");
    ctofdof.SetSize(ndof);
    // take space-fespaces coupling dofs for all time levels
    if (print)
      *testout << "SpaceTimeFESpace, ctofdof = " << endl << ctofdof << endl;
  }


  // ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    throw Exception("SpaceTimeFESpace :: GetFE - not implemented yet");
    // take space-fespaces GetFE and tp-extend with fel_time
  }


  // ------------------------------------------------------------------------
  const FiniteElement & SpaceTimeFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    throw Exception("SpaceTimeFESpace :: GetSFE - not implemented yet");
    // take space-fespaces GetSFE and tp-extend with fel_time
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
    throw Exception("SpaceTimeFESpace :: GetDofNrs - not implemented yet");
    // take space-fespaces GetDofNrs for each time level
  }

  // ------------------------------------------------------------------------
  void SpaceTimeFESpace :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    throw Exception("SpaceTimeFESpace :: GetSDofNrs - not implemented yet");
    // take space-fespaces GetSDofNrs for each time level
  }

  // ------------------------------------------------------------------------
  Table<int> * SpaceTimeFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    throw Exception("SpaceTimeFESpace :: CreateDirectSolverClusters - not implemented yet");
    Table<int> & table = *new Table<int> (0);
    return &table;

  }


  
  Array<int> * SpaceTimeFESpace :: CreateDirectSolverClusters (const Flags & precflags) const
  {
    throw Exception("SpaceTimeFESpace :: CreateDirectSolverClusters - not implemented yet");
    Array<int> & clusters = *new Array<int> (0);
    return & clusters;
  }
  
  // ------------------------------------------------------------------------

  
  //  static RegisterFESpace<SpaceTimeFESpace> init_facet ("spacetimefes");
}


