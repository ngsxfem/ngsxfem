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
  public:
    XStokesFESpace (shared_ptr<MeshAccess> ama,
                    const Flags & flags);
    virtual ~XStokesFESpace () { ; }
    virtual string GetClassName () const { return "XStokesFESpace"; }

    void SetLevelSet(shared_ptr<CoefficientFunction> lset){ 
      for (int i = 0; i < spaces.Size(); ++i)
        dynamic_pointer_cast<XStdFESpace>(spaces[i]) -> SetLevelSet(lset);
    };

    void SetLevelSet(shared_ptr<GridFunction> lset){ 
      SetLevelSet(make_shared<GridFunctionCoefficientFunction>(lset));
    };

    
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

