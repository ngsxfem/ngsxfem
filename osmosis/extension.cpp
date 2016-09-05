#include <solve.hpp>
#include "../xfem/xFESpace.hpp"
using namespace ngfem;
using namespace ngsolve;

/*
  first gf is gf coming from a fictitious domain approach with an "-empty" xfespace
  second gf lives on xfespace that is not "-empty"
  
  gridfunctions of nodes in the first layer of elements on the interior of gamma will not be visuali  zed correctly (no idea how to do this)
  
  xfes has to be the compound fespace with nonempy xfespace component
  
  is needed because xdof<->basedof mapping in XFESpace is protected (do not want to modify that file) - edit: i ended up having to modify it

  dt is the domain that is the "interior"
*/
void GetGFVisualization(shared_ptr<GridFunction> gf, shared_ptr<GridFunction> gf_vis, shared_ptr<FESpace> xfes, shared_ptr<MeshAccess> ma, DOMAIN_TYPE dt = NEG)
{
  if(dt==IF)
    {
      cout << "called GFVisualization with domain type IF !!" << endl;
      return;
    }

  shared_ptr<CompoundFESpace> comp_fes = dynamic_pointer_cast<CompoundFESpace>(xfes);
  shared_ptr<XFESpace> xfes_x = dynamic_pointer_cast<XFESpace>((*comp_fes)[1]);
  shared_ptr<FESpace> xfes_c = (*comp_fes)[0];
  
  int ne = ma->GetNE();
  int ndof_x = xfes_x->GetNDof();
  int ndof_c = xfes_c->GetNDof();

  BitArray isneg (comp_fes->GetNDof());
  isneg = 0;
  for(int i=0;i<ne;i++)
    if(xfes_x->GetDomainOfElement(i)==dt)
      {
	Array<int> dnums_c;
	xfes_c->GetDofNrs(i, dnums_c);
	for (int k : dnums_c)
	  isneg.Set(k);
      }
  for(int i=0;i<ne;i++)
    if(xfes_x->GetDomainOfElement(i)==IF)
      {
	Array<int> dnums_x;
	xfes_x->GetDofNrs(i, dnums_x);
	Array<int> dnums_c;
	xfes_c->GetDofNrs(i, dnums_c);
	for (int k : dnums_c)
	  isneg.Clear(k);
	for (int k : dnums_x)
	  if(xfes_x->GetDomainOfDof(k)==NEG)
	    isneg.Set(ndof_c+k);
      }

  //*testout << "isneg:" << endl << isneg << endl;

  int is_neg = 0;
  int is_neg2 = 0;
  for(int k=0;k<ndof_c;k++)
    if(isneg.Test(k))
      is_neg++;
  is_neg2 = is_neg;
  for(int k=0;k<ndof_x;k++)
    if(isneg.Test(ndof_c+k))
      is_neg2++;

  //*testout << "is_neg: " << is_neg << ", " << is_neg2 << endl;
  //*testout << "ndof_c: " << ndof_c << ", ndof_x: " << xfes_x->GetNDof() << endl;

  gf_vis->GetVector().FVDouble() = 0.0;

  for(int i=0;i<ne;i++)
    if(xfes_x->GetDomainOfElement(i)==dt)
      {
       
	Array<int> dnums_c;
	xfes_c->GetDofNrs(i, dnums_c);
	//continuous DOFs in the interior carry over!!
	
	//*testout << "element " << i << " with dofs " << dnums_c[0] << ", " << dnums_c[1]  << ", " << dnums_c[2] << "is in interior->write over values" << endl;
	for(int k=0;k<dnums_c.Size();k++)
	  if(isneg.Test(dnums_c[k]))
	    gf_vis->GetVector().FVDouble()[dnums_c[k]] = gf->GetVector().FVDouble()[dnums_c[k]];
      }
    else if(xfes_x->GetDomainOfElement(i)==IF)
      {
	Array<DOMAIN_TYPE> domnrs;
	xfes_x->GetDomainNrs(i,domnrs);
	Array<int> dnums_c;
	xfes_c->GetDofNrs(i, dnums_c);
	Array<int> dnums_x;
	xfes_x->GetDofNrs(i, dnums_x);

	//*testout << "element " << i << " with dofs " << dnums_c[0] << ", " << dnums_c[1]  << ", " << dnums_c[2] << "is on IF" << endl;

	for(int k=0;k<dnums_x.Size();k++)
	  if(isneg.Test(ndof_c+dnums_x[k]))//x-dof is interior->dof is exterior
	    {
	      //*testout << "dof nr " << xfes_x->x2base(dnums_x[k]) << " with x-dof " << dnums_x[k] << "/" << dnums_x[k]+ndof_c << " should be interior/exterior" << endl;
	      gf_vis->GetVector().FVDouble()[ndof_c+dnums_x[k]] = gf->GetVector().FVDouble()[xfes_x->x2base(dnums_x[k])];
	    }
	  else//x-dof is exterior->dof is interior
	    {
	      //*testout << "dof nr " << xfes_x->x2base(dnums_x[k]) << " with x-dof " << dnums_x[k] << "/" << dnums_x[k]+ndof_c << " should be exterior/interior" << endl;
	      
	      gf_vis->GetVector().FVDouble()[xfes_x->x2base(dnums_x[k])] = 1.0 * gf->GetVector().FVDouble()[xfes_x->x2base(dnums_x[k])];
	      gf_vis->GetVector().FVDouble()[ndof_c+dnums_x[k]] = -1.0 * gf->GetVector().FVDouble()[xfes_x->x2base(dnums_x[k])];
	      
	    }
	

      }
}


class NumProcTestExt : public NumProc
{
protected:
  
  shared_ptr<BilinearForm> bf_x;
  shared_ptr<BilinearForm> bf_c; 
  shared_ptr<BilinearForm> bf_extension; 
  shared_ptr<BilinearForm> bf_rob;
  

  shared_ptr<GridFunction> u_x;   //disc u
  shared_ptr<GridFunction> u_c;   //cont extension of u
  shared_ptr<GridFunction> u_vis; //"correct" visualization of u

  shared_ptr<LinearForm>   lf_x;

  shared_ptr<FESpace> xfes;
  shared_ptr<FESpace> xfes_dummy;
  shared_ptr<FESpace> cfes;

public:

  NumProcTestExt (shared_ptr<PDE> apde, const Flags & flags) : NumProc(apde)
  {
    *testout << endl << "start numproc test ext construct" << endl << endl;
    bf_x = apde->GetBilinearForm("bf_x");
    bf_c = apde->GetBilinearForm("bf_c");
    bf_extension= apde->GetBilinearForm("bf_extension");
    bf_rob = apde->GetBilinearForm("bf_rob");
    u_x = apde->GetGridFunction("u_d");
    u_c = apde->GetGridFunction("u_c");
    u_vis = apde->GetGridFunction("u_vis");
    xfes = apde->GetFESpace("xfes");
    cfes = apde->GetFESpace("cfes");
    lf_x = apde->GetLinearForm("lf_x");
    xfes_dummy = apde->GetFESpace("xfes_dummy");

    *testout << endl << "end numproc test ext construct" << endl << endl;
  }

  virtual ~NumProcTestExt(){ ; }


  virtual void Do(LocalHeap & lh)
  {
    shared_ptr<CompoundFESpace> comp_fes = dynamic_pointer_cast<CompoundFESpace>(xfes);
    shared_ptr<XFESpace> xfes_x = dynamic_pointer_cast<XFESpace>((*comp_fes)[1]);
    shared_ptr<FESpace> xfes_c = (*comp_fes)[0];

    int c_ndof = cfes->GetNDof();
    int x_ndof = xfes->GetNDof();
    int x_ndof_x = xfes_x->GetNDof();
    int x_ndof_c = xfes_c->GetNDof();

    shared_ptr<MeshAccess> ma = bf_x->GetMeshAccess();
    int ne = ma->GetNE();
    
    BitArray innerdofs(x_ndof);
    innerdofs=false;
    Array<int> dnums;
    for (int i : Range(ne))
      {
	if (xfes_x->GetDomainOfElement(i) != POS)
	  {
	    xfes_c->GetDofNrs(i, dnums);
	    for (int k : dnums)
	      innerdofs.Set(k);
	  }
      }
    *testout << "innerdofs:" << endl << innerdofs << endl;
    BitArray outerdofs(innerdofs);
    outerdofs.Invert();
    if ( xfes-> GetFreeDofs())
      {
	innerdofs.And(*xfes->GetFreeDofs());
	outerdofs.And(*xfes->GetFreeDofs());
      }

    BaseVector & vecu = u_x->GetVector();
    const BaseVector & vecf = lf_x->GetVector();

    BaseMatrix & matinner = bf_x->GetMatrix();
    BaseMatrix & matouter = bf_c->GetMatrix();
    
    shared_ptr<BaseVector> w = u_c->GetVector().CreateVector();
    
    shared_ptr<BaseMatrix> solveinner = matinner.InverseMatrix(&innerdofs);
    shared_ptr<BaseMatrix> solveouter = matouter.InverseMatrix(&outerdofs);
    vecu = *solveinner * vecf;
    GetGFVisualization(u_x,u_vis,xfes, ma, NEG);
    for(int k=0;k<x_ndof_c;k++)
     u_c->GetVector().FVDouble()[k] = vecu.FVDouble()[k];
    *w = -1.0 * matouter * u_c->GetVector();
    u_c->GetVector()+= *solveouter * *w;
    
    
  }    
};


static RegisterNumProc<NumProcTestExt> init_npte ("test_ext");


