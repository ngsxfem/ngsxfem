
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
#include "../utils/p1interpol.hpp"

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
    cout << IM(3) << "Constructor of SpaceTimeFESpace" << endl;
    cout << IM(3) <<"Flags = " << flags << endl;

    dimension = aVh->GetDimension ();

    int order_s = aVh->GetOrder();
    int order_t = atfe->Order();
    bool linear_time = order_t == 1;

    Vh = aVh.get();
    tfe = atfe.get();

    cout << IM(3) <<"Hello from SpaceTimeFESpace.cpp" << endl;
    cout << IM(3) <<"Order Space: " << order_s << endl;
    cout << IM(3) <<"Order Time: " << order_t << endl;

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
    cout << IM(3) << "Dofs in base: " << Vh->GetNDof() << endl;

    // number of dofs:
    ndof = (Vh->GetNDof()) * tfe->GetNDof();
    cout << IM(3) << "Total number of Dofs: " << ndof << endl;


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

  void SpaceTimeFESpace :: RestrictGFInTime(shared_ptr<GridFunction> st_GF, double time, shared_ptr<GridFunction> s_GF)
  {
     auto st_vec = st_GF->GetVectorPtr()->FV<double>();   
     auto restricted_vec = s_GF->GetVectorPtr()->FV<double>();
     restricted_vec = 0.0;

     Vector<double> nodes(order_time() +1 );
     TimeFE_nodes(nodes);

     //cout << "Vhdim = " << Vh->GetDimension() << endl;
     // Using nodal property for special case
     for(int i= 0; i < nodes.Size(); i++) {
         if(time == nodes[i]) {
             cout << IM(3) <<"Node case" << endl;
             for(int j = 0; j < Vh->GetNDof();j++)
                 restricted_vec[j] = st_vec[j+i*Vh->GetNDof()];
             return;
         }
     }
     
     cout << IM(3) <<"General case" << endl;
     // General case
     //cout << IM(3) <<"time fe:" << GetTimeFE() << endl;
     NodalTimeFE * time_FE = dynamic_cast<NodalTimeFE*>(tfe);
     const int dim = Vh->GetDimension();     
     for(int i= 0; i < nodes.Size(); i++) {
         double weight_time = time_FE->Lagrange_Pol(time,nodes,i);
         for(int j = 0; j < Vh->GetNDof();j++)
             for(int d = 0; d < dim;d++)
                 restricted_vec[dim*j+d] += weight_time * st_vec[dim*(j+i*Vh->GetNDof())+d];
     }
  }


  shared_ptr<GridFunction> SpaceTimeFESpace :: CreateRestrictedGF(shared_ptr<GridFunction> st_GF, double time)
  {
     shared_ptr<GridFunction> restricted_GF = nullptr;
     
     switch (Vh->GetDimension())
     {
       case 1:
         restricted_GF = make_shared < T_GridFunction < double > >( Vh_ptr);
         break;
       case 2:
         restricted_GF = make_shared < T_GridFunction < Vec<2> > >( Vh_ptr);
         break;
     
       default:
         throw Exception("cannot handle GridFunction type (dimension too large?).");
         break;
     }
     restricted_GF->Update();
     RestrictGFInTime(st_GF, time, restricted_GF);
     return restricted_GF;
  }

  void SpaceTimeFESpace ::InterpolateToP1(shared_ptr<CoefficientFunction> st_CF, shared_ptr<CoefficientFunction> ctref, double t, double dt, shared_ptr<GridFunction> st_GF)
  {
    LocalHeapMem<100000> lh("SpacetimeInterpolateToP1");
    auto node_gf = make_shared < T_GridFunction < double > >( Vh_ptr);
    node_gf->Update();
    auto gf_vec = st_GF->GetVectorPtr()->FV<double>();
    auto node_gf_vec = node_gf->GetVectorPtr()->FV<double>();
    shared_ptr<ParameterCoefficientFunction> coef_tref = dynamic_pointer_cast<ParameterCoefficientFunction>(ctref);
    if (!coef_tref)
      throw Exception("SpaceTimeFESpace ::InterpolateToP1 : tref is not a ParameterCF");
    const double backup_tref = coef_tref->GetValue();
    Vector<double> nodes(order_time() +1 );
    TimeFE_nodes(nodes);
    for(int i= 0; i < nodes.Size(); i++) {
      coef_tref->SetValue(t+nodes[i]*dt);
      InterpolateP1 iP1(st_CF, node_gf);
      iP1.Do(lh);
      for(int j = 0; j < Vh_ptr->GetNDof();j++)
      {
        gf_vec(i*Vh_ptr->GetNDof()+j) = node_gf_vec(j);
      }
    }        
    coef_tref->SetValue(backup_tref);
  }

}
