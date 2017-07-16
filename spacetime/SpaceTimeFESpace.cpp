
// SpaceTimeFESpace based on:

/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>

#include <fem.hpp>


#include "SpaceTimeFE.hpp"
#include "SpaceTimeFESpace.hpp"

/*
#include <diffop_impl.hpp>
#ifdef WIN32
      template ngcomp::T_DifferentialOperator<ngcomp::DiffOpId<2> >;
#endif
*/

namespace ngcomp
{

SpaceTimeFESpace :: SpaceTimeFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> aVh, shared_ptr<ScalarFiniteElement<1>> atfe, const Flags & flags)
  : FESpace (ama, flags), Vh_ptr(aVh)
  {
    cout << "Constructor of SpaceTimeFESpace" << endl;
    cout << "Flags = " << flags << endl;

    dimension = aVh->GetDimension ();

    int order_s = aVh->GetOrder();
    int order_t = atfe->Order();
    bool linear_time = order_t == 1;

    Vh = aVh.get();
    tfe = atfe.get();

    cout << "Hello from SpaceTimeFESpace.cpp" << endl;
    cout << "Order Space: " << order_s << endl;
    cout << "Order Time: " << order_t << endl;

    // needed to draw solution function
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
                                                 make_shared<ConstantCoefficientFunction>(1));

    if (dimension > 1)
    {
      evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
      flux_evaluator[VOL] = make_shared<BlockDifferentialOperator> (flux_evaluator[VOL], dimension);
      evaluator[BND] = 
        make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
      // flux_evaluator[BND] = 
      //   make_shared<BlockDifferentialOperator> (flux_evaluator[BND], dimension);
    }

    time=0;
  }


  SpaceTimeFESpace :: ~SpaceTimeFESpace ()
  {
    // nothing to do
  }


  void SpaceTimeFESpace :: Update(LocalHeap & lh)
  {
    // some global update:
    if(dirichlet_boundaries.Size() == 0) {
      dirichlet_boundaries.SetSize(ma->GetNBoundaries());
      dirichlet_boundaries.Clear();
      for(int i = 0; i < ma->GetNBoundaries();i++) {
          if(Vh->IsDirichletBoundary(i))
            dirichlet_boundaries.Set(i);
       }
    }
    FESpace::Update(lh);
    Vh->Update(lh);
    cout << "Dofs in base: " << Vh->GetNDof() << endl;

    // number of dofs:
    ndof = (Vh->GetNDof()) * tfe->GetNDof();
    cout << "Total number of Dofs: " << ndof << endl;


  }

  void SpaceTimeFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {

    dnums.SetSize(0);

    Array<int> Vh_dofs;
    Vh->GetDofNrs(ei,Vh_dofs);

    for( int i = 0; i < tfe->GetNDof(); i++) {
        for (auto v : Vh_dofs)
            dnums.Append (v+i*Vh->GetNDof());
     }

  }


  FiniteElement & SpaceTimeFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {

     ScalarFiniteElement<2>* s_FE = dynamic_cast<ScalarFiniteElement<2>*>(&(Vh->GetFE(ei,alloc)));

     ScalarFiniteElement<1>* t_FE = tfe;

     SpaceTimeFE * st_FE =  new (alloc) SpaceTimeFE(s_FE,t_FE,override_time,time);

     return *st_FE;

   }


  shared_ptr<GridFunction> SpaceTimeFESpace :: CreateRestrictedGF(shared_ptr<GridFunction> st_GF, double time)
  {

     //auto restricted_GF = make_shared < T_GridFunction < double > >( *Vh) ;
     auto restricted_GF =  make_shared < T_GridFunction < double > >( Vh_ptr);
     auto st_vec = st_GF->GetVectorPtr()->FV<double>();
     restricted_GF->Update();
     auto restricted_vec = restricted_GF->GetVectorPtr()->FV<double>();

     if(time == 0.0) {
         for(int i = 0; i < Vh->GetNDof(); i++)
             restricted_vec[i] = st_vec[i];
     } else {
         cout << "Only t = 0 implemented so far." << endl;
     }

     return restricted_GF;

  }


}
