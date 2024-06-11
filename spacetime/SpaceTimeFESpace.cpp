
/*********************************************************************/
/* File:   SpaceTimeFESpace.cpp                                      */
/* Author: Janosch Preuss & Christoph Lehrenfeld                     */
/* Date:   2017                                                      */
/*********************************************************************/

#include "SpaceTimeFE.hpp"
#include <nginterface.h>
#include "SpaceTimeFESpace.hpp"

#include "../utils/p1interpol.hpp"
#include "../utils/ngsxstd.hpp"

#include "timecf.hpp"
#include "diffopDt.hpp"

/*
#include <diffop_impl.hpp>
#ifdef WIN32
      template ngcomp::T_DifferentialOperator<ngcomp::DiffOpId<2> >;
#endif
*/

namespace ngcomp
{

SpaceTimeFESpace :: SpaceTimeFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> aVh, shared_ptr<ScalarFiniteElement<1>> atfe, const Flags & flags)
  : FESpace (ama, flags), Vh(aVh), tfe(atfe)
  {
    *testout << "AMA DIM: " << ama->GetDimension() << endl;
    *testout << "Constructor of SpaceTimeFESpace" << endl;
    *testout <<"Flags = " << flags << endl;

    dimension = Vh->GetDimension ();

    int order_s = Vh->GetOrder();
    int order_t = atfe->Order();
    bool linear_time = order_t == 1;

    *testout <<"Hello from SpaceTimeFESpace.cpp" << endl;
    *testout <<"Order Space: " << order_s << endl;
    *testout <<"Order Time: " << order_t << endl;

    // needed to draw solution function
    Switch<3> (ma->GetDimension()-1, [&] (auto SDIM) {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<SDIM+1>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<SDIM+1>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<SDIM+1>>>();
    });
    //if (ma->GetDimension() < 2)
    //    throw Exception ("Unsupported spatial dimension in SpaceTimeFESpace :: SpaceTimeFESpace");

     integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
                                                 make_shared<ConstantCoefficientFunction>(1));

    if (dimension > 1)
    {
      evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
      flux_evaluator[VOL] = make_shared<BlockDifferentialOperator> (flux_evaluator[VOL], dimension);
      evaluator[BND] = 
        make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
      Switch<2> (ma->GetDimension()-2, [&] (auto SDIM) {
        enum {SDIM2 = SDIM()+2}; // work around MSVC compile issue
        Switch<3> (dimension-1, [&] (auto DIM) {
          additional_evaluators.Set ("dt", make_shared<T_DifferentialOperator<DiffOpDtVec<SDIM2,DIM()+1>>>());
          additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<DIM()+1>>> ());
        });
      });
    }
    else
      Switch<3> (ma->GetDimension()-1, [&] (auto DIM) {
        additional_evaluators.Set ("dt", make_shared<T_DifferentialOperator<DiffOpDt<DIM+1>>>());
        additional_evaluators.Set ("fix_tref_bottom", make_shared<T_DifferentialOperator<DiffOpFixt<DIM+1,0>>>());
        additional_evaluators.Set ("fix_tref_top", make_shared<T_DifferentialOperator<DiffOpFixt<DIM+1,1>>>());
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<DIM+1>>> ());
      });
    
        



    time=0;
  }


  SpaceTimeFESpace :: ~SpaceTimeFESpace ()
  {
    // nothing to do
  }


  void SpaceTimeFESpace :: Update()
  {
    // some global update:
    if(dirichlet_boundaries.Size() == 0) {
      dirichlet_boundaries.SetSize(ma->GetNBoundaries());
      dirichlet_boundaries.Clear();
      for(int i = 0; i < ma->GetNBoundaries();i++) {
          if(Vh->IsDirichletBoundary(i))
            dirichlet_boundaries.SetBit(i);
       }
    }
    FESpace::Update();
    Vh->Update();
    *testout << "Dofs in base: " << Vh->GetNDof() << endl;

    // number of dofs:
    ndof = (Vh->GetNDof()) * tfe->GetNDof();
    *testout << "Total number of Dofs: " << ndof << endl;


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
     ScalarFiniteElement<1>* t_FE = tfe.get();
     if(ma->GetDimension() == 2){
       ScalarFiniteElement<2>* s_FE2 = dynamic_cast<ScalarFiniteElement<2>*>(&(Vh->GetFE(ei,alloc)));
       //cout << "SpaceTimeFESpace :: GetFE for 2D called" << endl;
       SpaceTimeFE<2> * st_FE =  new (alloc) SpaceTimeFE<2>(s_FE2,t_FE,override_time,time);
       return *st_FE;
     }
     else if (ma->GetDimension() == 3){
       ScalarFiniteElement<3>* s_FE3 = dynamic_cast<ScalarFiniteElement<3>*>(&(Vh->GetFE(ei,alloc)));
       //cout << "SpaceTimeFESpace :: GetFE for 3D called" << endl;
       SpaceTimeFE<3> * st_FE =  new (alloc) SpaceTimeFE<3>(s_FE3,t_FE,override_time,time);
       return *st_FE;
     }
     else if (ma->GetDimension() == 1){
       ScalarFiniteElement<1>* s_FE1 = dynamic_cast<ScalarFiniteElement<1>*>(&(Vh->GetFE(ei,alloc)));
       //cout << "SpaceTimeFESpace :: GetFE for 3D called" << endl;
       SpaceTimeFE<1> * st_FE =  new (alloc) SpaceTimeFE<1>(s_FE1,t_FE,override_time,time);
       return *st_FE;
     }
     else
       throw Exception("SpaceTimeFESpace :: GetFE cannot help dimension != 1,2,3");
   }

  template<typename SCAL>
  void SpaceTimeFESpace :: RestrictGFInTime(shared_ptr<GridFunction> st_GF, double time, shared_ptr<GridFunction> s_GF)
  {
     auto st_vec = st_GF->GetVectorPtr()->FV<SCAL>();   
     auto restricted_vec = s_GF->GetVectorPtr()->FV<SCAL>();

     Array<double> & nodes = TimeFE_nodes();

     //cout << "Vhdim = " << Vh->GetDimension() << endl;
     // Using nodal property for special case
     int cnt = 0;
     for(int i= 0; i < nodes.Size(); i++) {
         if (!IsTimeNodeActive(i))
         {
           if (abs(time - nodes[i]) < globxvar.EPS_STFES_RESTRICT_GF)
           {
             restricted_vec = SCAL(0.0);
             return;
           }
           else
           {
             continue;
           }
         }
           
         if(abs(time - nodes[i]) < globxvar.EPS_STFES_RESTRICT_GF) {
             *testout <<"Node case" << endl;
             for(int j = 0; j < Vh->GetNDof();j++)
                 restricted_vec[j] = st_vec[j+cnt*Vh->GetNDof()];
             return;
         }
         cnt++;
     }
     *testout <<"General case" << endl;
     // General case
     //*testout <<"time fe:" << GetTimeFE() << endl;
     shared_ptr<NodalTimeFE> time_FE = dynamic_pointer_cast<NodalTimeFE>(tfe);

     const int dim = Vh->GetDimension();
     for(int j=0; j< Vh->GetNDof(); j++) restricted_vec[j] = 0.;

     for(int i= 0; i < nodes.Size(); i++) {
       if (IsTimeNodeActive(i))
       {
         double weight_time = time_FE->Lagrange_Pol(time,i);
         for(int j = 0; j < Vh->GetNDof();j++)
             //for(int d = 0; d < dim;d++)
                 //restricted_vec[dim*j+d] += weight_time * st_vec[dim*(j+i*Vh->GetNDof())+d];
             restricted_vec[j] += weight_time * st_vec[j+i*Vh->GetNDof()];
       }
     }
  }

  
  shared_ptr<GridFunction> SpaceTimeFESpace :: CreateRestrictedGF(shared_ptr<GridFunction> st_GF, double time)
  {
     shared_ptr<GridFunction> restricted_GF = nullptr;
     restricted_GF = make_shared < S_GridFunction < double > >( Vh);
     restricted_GF->Update();
     switch (Vh->GetDimension())
     {
       case 1: RestrictGFInTime<double>(st_GF, time, restricted_GF); break;
       case 2: RestrictGFInTime<Vec<2>>(st_GF, time, restricted_GF); break;
       case 3: RestrictGFInTime<Vec<3>>(st_GF, time, restricted_GF); break;
       default: throw Exception("cannot handle GridFunction type (dimension too large?)."); break;
     }

     return restricted_GF;
  }

  void SpaceTimeFESpace ::InterpolateToP1(shared_ptr<CoefficientFunction> st_CF, shared_ptr<CoefficientFunction> ctref, shared_ptr<GridFunction> st_GF)
  {
    LocalHeapMem<100000> lh("SpacetimeInterpolateToP1");
    auto node_gf = make_shared < S_GridFunction < double > >( Vh);
    node_gf->Update();
    auto gf_vec = st_GF->GetVectorPtr()->FV<double>();
    auto node_gf_vec = node_gf->GetVectorPtr()->FV<double>();
    shared_ptr<TimeVariableCoefficientFunction> coef_tref = dynamic_pointer_cast<TimeVariableCoefficientFunction>(ctref);
    if (!coef_tref)
      throw Exception("SpaceTimeFESpace ::InterpolateToP1 : tref is not a TimeVariableCoefficientFunction");
    //const double backup_tref = coef_tref->GetValue();
    Array<double> & nodes = TimeFE_nodes();
    for(int i= 0; i < nodes.Size(); i++) {
      if (IsTimeNodeActive(i))
      {
        //coef_tref->SetValue(t+nodes[i]*dt);
          coef_tref->FixTime(nodes[i]);

        InterpolateP1 iP1(st_CF, node_gf);
        iP1.Do(lh, globxvar.EPS_INTERPOLATE_TO_P1, nodes[i]);
        for(int j = 0; j < Vh->GetNDof();j++)
        {
          gf_vec(i*Vh->GetNDof()+j) = node_gf_vec(j);
        }
      }
    }        
    coef_tref->UnfixTime();
    //coef_tref->SetValue(backup_tref);
  }


  template void SpaceTimeFESpace :: RestrictGFInTime<double>(shared_ptr<GridFunction>, double, shared_ptr<GridFunction>);
  template void SpaceTimeFESpace :: RestrictGFInTime<Vec<2>>(shared_ptr<GridFunction>, double, shared_ptr<GridFunction>);
  template void SpaceTimeFESpace :: RestrictGFInTime<Vec<3>>(shared_ptr<GridFunction>, double, shared_ptr<GridFunction>);

}
