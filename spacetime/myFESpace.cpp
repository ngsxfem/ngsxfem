
/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

A fe-space provides the connection between the local reference
element, and the global mesh.

*/


#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>

#include <fem.hpp>


#include "myElement.hpp"
#include "myHOElement.hpp"
#include "myFESpace.hpp"

/*
#include <diffop_impl.hpp>
#ifdef WIN32
      template ngcomp::T_DifferentialOperator<ngcomp::DiffOpId<2> >;
#endif
*/

namespace ngcomp
{

  //SpaceTimeFESpace :: SpaceTimeFESpace (shared_ptr<MeshAccess> ama, FESpace& aVh, ScalarFiniteElement<1>& atfe, const Flags & flags)
  //  : FESpace (ama, flags)
SpaceTimeFESpace :: SpaceTimeFESpace (shared_ptr<MeshAccess> ama, shared_ptr<FESpace> aVh, shared_ptr<ScalarFiniteElement<1>> atfe, const Flags & flags)
  : FESpace (ama, flags)
  {
    cout << "Constructor of MyFESpace" << endl;
    cout << "Flags = " << flags << endl;

    order_s = int(flags.GetNumFlag ("order", 2));

    // How to get additional NumFlags?
    bool linear_time = flags.GetDefineFlag ("order_time");
    order_t = 1;

    Vh = aVh.get();
    tfe = atfe.get();

    cout << "Hello from myFESpace.cpp" << endl;
    cout << "Order Space: " << order_s << endl;
    cout << "Order Time: " << order_t << endl;

    // needed to draw solution function
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
                                                 make_shared<ConstantCoefficientFunction>(1));
    time=0;
  }


  SpaceTimeFESpace :: ~SpaceTimeFESpace ()
  {
    // nothing to do
  }


  void SpaceTimeFESpace :: Update(LocalHeap & lh)
  {
    // some global update:

    Vh->Update(lh);
    cout << "Dofs in space: " << Vh->GetNDof() << endl;

    // number of dofs:
    ndof = (Vh->GetNDof()) * tfe->GetNDof();


  }

  void SpaceTimeFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment

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

     // Ok to do this ?
     ScalarFiniteElement<2>* s_FE = dynamic_cast<ScalarFiniteElement<2>*>(&(Vh->GetFE(ei,alloc)));

     ScalarFiniteElement<1>* t_FE = tfe;

     SpaceTimeFE * st_FE =  new (alloc) SpaceTimeFE(order_s,s_FE,t_FE,time);

     return *st_FE;

   }




 // ************************************************************** //

  // Registering the SpaceTimeFESpace

  // From MylittleNGSolve
  // static RegisterFESpace<SpaceTimeFESpace> initifes ("SpaceTimeFESpace");


  /*
  template <typename FES>
  class MyRegisterFESpace
  {
  public:
    /// constructor registers fespace
    MyRegisterFESpace (string label)
    {
      GetFESpaceClasses().AddFESpace (label, Create);
      // cout << "register fespace '" << label << "'" << endl;
    }

    /// creates an fespace of type FES
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<FES> (ma, flags);
    }
  };

    string label = "SpaceTimeFESpace";
    static MyRegisterFESpace<SpaceTimeFESpace> initifes (label);

    //  error: no matching function for call to
    //   'ngcomp::SpaceTimeFESpace::SpaceTimeFESpace(std::shared_ptr<ngcomp::MeshAccess>&, const ngstd::Flags&)'

    */




 /*
  template <typename FES>
  class MyRegisterFESpace
  {
  public:
    /// constructor registers fespace
    MyRegisterFESpace (string label)
    {
      GetFESpaceClasses().AddFESpace (label, Create);
      // cout << "register fespace '" << label << "'" << endl;
    }

    /// creates an fespace of type FES
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma,FESpace& aVh, ScalarFiniteElement<1>& atfe, const Flags & flags)
    {
      return make_shared<FES> (ma, aVh, atfe, flags);
    }
  };

    string label = "SpaceTimeFESpace";
    static MyRegisterFESpace<SpaceTimeFESpace> initifes (label);

 // invalid conversion from 'std::shared_ptr<ngcomp::FESpace> (*)(std::shared_ptr<ngcomp::MeshAccess>, ngcomp::FESpace&, ngfem::ScalarFiniteElement<1>&, const ngstd::Flags&)'
 //    to 'std::shared_ptr<ngcomp::FESpace> (*)(std::shared_ptr<ngcomp::MeshAccess>, const ngstd::Flags&)' [-fpermissive]
 //      GetFESpaceClasses().AddFESpace (label, Create);


  */



}
