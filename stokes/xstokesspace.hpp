#ifndef FILE_XSTOKESSPACE_HPP
#define FILE_XSTOKESSPACE_HPP

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement

#include "../xfem/xFESpace.hpp"

using namespace ngsolve;

namespace ngcomp
{

  class XStokesFESpace : public CompoundFESpace
  {
    shared_ptr<CoefficientFunction> lset = nullptr;
  public:
    XStokesFESpace (shared_ptr<MeshAccess> ama,
                    const Flags & flags);
    virtual ~XStokesFESpace () { ; }
    virtual string GetClassName () const { return "XStokesFESpace"; }

    void SetLevelSet(shared_ptr<CoefficientFunction> alset){
      lset = alset;
      for (int i = 0; i < spaces.Size(); ++i)
        dynamic_pointer_cast<XStdFESpace>(spaces[i]) -> SetLevelSet(alset);
    };

    void SetLevelSet(shared_ptr<GridFunction> lset){ 
      SetLevelSet(make_shared<GridFunctionCoefficientFunction>(lset));
    };

    virtual void Update(LocalHeap & lh)
    {
      if ( lset == nullptr )
        cout << "xstokesspace: no lset, Update postponed " << endl;
      else
        CompoundFESpace::Update(lh);
    }    
    
    virtual void FinalizeUpdate(LocalHeap & lh)
    {
      if ( lset == nullptr )
        cout << "xstokesspace: no lset, Update postponed " << endl;
      else
        CompoundFESpace::FinalizeUpdate(lh);
    }    
    
  };


  class NumProcInformXStokesFESpace : public NumProc
  {
    const CoefficientFunction * coef = NULL;
  public:
    NumProcInformXStokesFESpace (shared_ptr<PDE> apde, const Flags & flags);
    ~NumProcInformXStokesFESpace();
    virtual string GetClassName () const;
    virtual void Do (LocalHeap & lh);
  };


}

#endif

