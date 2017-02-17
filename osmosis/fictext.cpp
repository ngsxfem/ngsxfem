#include <solve.hpp>
#include "../xfem/xFESpace.hpp"
using namespace ngfem;
using namespace ngsolve;

class NumProcSolveInnerAndExtend : public NumProc
{
protected:
  shared_ptr<BilinearForm> bfinner;
  shared_ptr<BilinearForm> bfouter;
  shared_ptr<LinearForm> lff;
  shared_ptr<GridFunction> gfu;
public:
  NumProcSolveInnerAndExtend (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  { 
    bfinner = apde->GetBilinearForm (flags.GetStringFlag ("bilinearforminner", "a1"));
    bfouter = apde->GetBilinearForm (flags.GetStringFlag ("bilinearformouter", "a2"));
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
  }

  virtual ~NumProcSolveInnerAndExtend()
  { ; }


  virtual void Do(LocalHeap & lh)
  {
    cout << "Solve Inner and extend" << endl;

    shared_ptr<FESpace> fes = bfinner->GetFESpace();

    shared_ptr<XFESpace> xfes = dynamic_pointer_cast<XFESpace>((*dynamic_pointer_cast<CompoundFESpace>(fes))[1]);

    shared_ptr<BitArray> innerdofs = make_shared<BitArray>(fes->GetNDof());
    innerdofs->Clear();

    shared_ptr<MeshAccess> ma = bfinner->GetMeshAccess();
    int ne = ma->GetNE();

    Array<int> dnums;
    for (int i : Range(ne))
    {
      if (xfes->GetDomainOfElement(i) != POS)
      {
        fes->GetDofNrs(ElementId(VOL,i), dnums);
        for (int k : dnums)
          innerdofs->Set(k);
      }
    }

    shared_ptr<BitArray> outerdofs = make_shared<BitArray>(*innerdofs);
    outerdofs->Invert();
    if ( fes-> GetFreeDofs())
    {
      innerdofs->And(*fes->GetFreeDofs());
      outerdofs->And(*fes->GetFreeDofs());
    }
    
    cout << " innerdofs = \n" << *innerdofs << endl;
    cout << " outdofs = \n" << *outerdofs << endl;
    
    BaseVector & vecu = gfu->GetVector();
    const BaseVector & vecf = lff->GetVector();

    shared_ptr<BaseVector> w = vecu.CreateVector();

    BaseMatrix& matinner = bfinner->GetMatrix();
    BaseMatrix& matouter = bfouter->GetMatrix();

    shared_ptr<BaseMatrix> solveinner = matinner.InverseMatrix(innerdofs);

    vecu = *solveinner * vecf;




    shared_ptr<BaseMatrix> solveouter = matouter.InverseMatrix(outerdofs);

    *w = -1.0 *  matouter * vecu;
    
    //vecu += *solveouter * *w;
  }
};



static RegisterNumProc<NumProcSolveInnerAndExtend> npinit("solveinnerandextend");
