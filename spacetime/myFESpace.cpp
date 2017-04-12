/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

My own FESpace for linear and quadratic triangular elements.  

A fe-space provides the connection between the local reference
element, and the global mesh.

*/


#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>


#include "myElement.hpp"
#include "myFESpace.hpp"

/*
#include <diffop_impl.hpp>
#ifdef WIN32
      template ngcomp::T_DifferentialOperator<ngcomp::DiffOpId<2> >;
#endif
*/

namespace ngcomp
{

  MyFESpace :: MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyFESpace" << endl;
    cout << "Flags = " << flags << endl;

    order = int(flags.GetNumFlag ("order", 2));

    cout << "Hello from myFESpace.cpp" << endl;
    cout << "Order: " << order << endl;
 
    if ( order > 2 ) {
       cout << "Sorry, only first and second order elements are implemented. Switch to second order." << endl;
       order = 2;
    }
   
    // needed to draw solution function
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(), 
                                                 make_shared<ConstantCoefficientFunction>(1));
  }
    
  
  MyFESpace :: ~MyFESpace ()
  {
    // nothing to do
  }

  
  void MyFESpace :: Update(LocalHeap & lh)
  {
    // some global update: 

    cout << "Update MyFESpace, #vert = " << ma->GetNV() 
         << ", #edge = " << ma->GetNEdges() << endl;

    int n_vert = ma->GetNV();
    int n_edge =  ma->GetNEdges();
    int n_cell =  ma->GetNE();

    // number of dofs:
    if (order == 1)
          ndof = n_vert;

    if( order == 2) {
       first_edge_dof.SetSize(n_edge+1);
       int ii = n_vert;
       for ( int i=0; i < n_edge; i++, ii+=1)
          first_edge_dof[i] = ii;
       first_edge_dof[n_edge] = ii;
    
       first_cell_dof.SetSize(n_cell+1);
       for ( int i = 0; i< n_cell; i++, ii+=1) 
         first_cell_dof[i] = ii;
       first_cell_dof[n_cell] = ii;
       cout << "first_cell_dof = " << endl << first_cell_dof << endl;
    
       ndof = ii;
    }
       
  }

  void MyFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment
    
    dnums.SetSize(0);

    Ngs_Element ngel = ma->GetElement (ei);


    Array<int> vert_nums;
    ma->GetElVertices (ei, vert_nums);  // global vertex numbers


        for (auto v : vert_nums) 
      dnums.Append (v);

    if (order == 2)
      { 
        for (auto e : ngel.Edges())
	  {
	    int first = first_edge_dof[e];
            int next  = first_edge_dof[e+1];
            for ( int j = first; j < next; j++) 
              dnums.Append(j);
           }
        if (ei.IsVolume())
          {
            int first = first_cell_dof[ei.Nr()];
            int next  = first_cell_dof[ei.Nr()+1];
            for ( int j = first; j < next; j++)
              dnums.Append(j);
          }
    }
  }
  
  FiniteElement & MyFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  { 
    Ngs_Element ngel = ma->GetElement (ei);
    if (ei.IsVolume())
      {
          MyLinearQuad * quad = new (alloc) MyLinearQuad(order);
		
          for (int i = 0; i < 4; i++) 
            quad->SetVertexNumber (i, ngel.vertices[i]);
	
      	    return *quad;
            
      }              
    else
      {
        return * new (alloc) FE_Segm1;
      }
  }



  /*
    register fe-spaces
    Object of type MyFESpace can be defined in the pde-file via
    "define fespace v -type=myfespace"
  */

  static RegisterFESpace<MyFESpace> initifes ("myfespace");
}
