/*********************************************************************/
/* File:   setdeformation.cpp                                        */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   21. July 2015                                             */
/*********************************************************************/


#include <solve.hpp>

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

namespace ngcomp
{ 

  class NumProcUnsetDeformation : public NumProc
  {
  public:
    shared_ptr<MeshAccess> ma;
    NumProcUnsetDeformation (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      ma = apde->GetMeshAccess();
    }

    virtual ~NumProcUnsetDeformation()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcUnsetDeformation";
    }


    virtual void Do (LocalHeap & clh)
    {
      ma->SetDeformation(nullptr);
    }
  };

  class NumProcSetDeformation : public NumProc
  {
    shared_ptr<GridFunction> gf_deform;
  public:
    shared_ptr<MeshAccess> ma;
    NumProcSetDeformation (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      gf_deform = apde->GetGridFunction(flags.GetStringFlag("gridfunction","deform"));
      ma = apde->GetMeshAccess();
    }

    virtual ~NumProcSetDeformation()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcSetDeformation";
    }


    virtual void Do (LocalHeap & clh)
    {
      ma->SetDeformation(gf_deform);
    }
  };

  
}

static RegisterNumProc<NumProcUnsetDeformation > npxunsetdeform("unsetdeformation");
static RegisterNumProc<NumProcSetDeformation > npxsetdeform("setdeformation");
