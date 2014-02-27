/*
  Solver for a parabolic pde

  Solves
  M du/dt  +  A u = f
  by an implicite Euler method

  Please include this file to the src files given in netgen/ngsolve/Makefile
*/

#define FILE_SOLVEINSTATX_CPP
#define FILE_SPACETIMEINT_CPP  

#include <solve.hpp>
#include "../spacetime/spacetimeintegrators.hpp"
#include "../xfem/stxfemIntegrators.hpp"
#include "../xfem/xfemNitsche.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/setvaluesx.hpp"

#include <diffop_impl.hpp>

// #include "../utils/stcoeff.hpp"

using namespace ngsolve;



/*
  Every solver is a class derived from the class NumProc.
  It collects objects (such as bilinear-forms, gridfunctions) and parameters.
*/
template <int D>
class NumProcSolveInstatX : public NumProc
{
protected:
  // // bilinear-form for matrix
  //BilinearForm * bftau;
  // // bilinear-form for the mass-matrix
  // BilinearForm * bfm;
  // linear-form providing the right hand side
  // LinearForm * lff; //initial condition
  // solution vector
  
  FESpace * fes;
  FESpace * fesvis;

  GridFunction * gfu;

  GridFunction * gfu_vis = NULL;

  // time step
  double dt;
  // total time
  double tend;
  // direct solver type
  string inversetype;
  string fesstr;
	
  bool userstepping;

  double bneg = 1.0;
  double bpos = 1.0;
  double aneg = 1.0;
  double apos = 1.0;
  double lambda = 10.0;

  CoefficientFunction* coef_bconvneg = NULL;
  CoefficientFunction* coef_bconvpos = NULL;

  CoefficientFunction* coef_bneg = NULL;
  CoefficientFunction* coef_bpos = NULL;
  CoefficientFunction* coef_aneg = NULL;
  CoefficientFunction* coef_apos = NULL;
  CoefficientFunction* coef_abneg = NULL;
  CoefficientFunction* coef_abpos = NULL;

  CoefficientFunction* coef_lambda = NULL;

  CoefficientFunction* coef_zero = NULL;
  CoefficientFunction* coef_one = NULL;

  CoefficientFunction* coef_told = NULL;
  CoefficientFunction* coef_tnew = NULL;

  CoefficientFunction* coef_upwneg = NULL;
  CoefficientFunction* coef_upwpos = NULL;

  CoefficientFunction* coef_brhsneg = NULL;
  CoefficientFunction* coef_brhspos = NULL;

  CoefficientFunction* coef_binineg = NULL;
  CoefficientFunction* coef_binipos = NULL;

  CoefficientFunction* coef_bndneg = NULL;
  CoefficientFunction* coef_bndpos = NULL;

public:
  /*
	In the constructor, the solver class gets the flags from the pde - input file.
	the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcSolveInstatX (PDE & apde, const Flags & flags) : NumProc (apde)
  {
	// in the input-file, you specify the bilinear-forms for the stiffness and for the mass-term
	// like  "-bilinearforma=k". Default arguments are 'a' and 'm'

	// bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
	// bfm = pde.GetBilinearForm (flags.GetStringFlag ("bilinearformm", "m"));
	// lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));

	gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
	gfu_vis = pde.GetGridFunction (flags.GetStringFlag ("gf_vis", "u_vis"),true);

	inversetype = flags.GetStringFlag ("solver", "pardiso");
	fesstr = flags.GetStringFlag ("fespace", "fes_st");
	fes = pde.GetFESpace (fesstr.c_str());

	string fesvisstr = flags.GetStringFlag ("fespacevis", "fes_negpos");
	fesvis = pde.GetFESpace (fesvisstr.c_str());

	userstepping = flags.GetDefineFlag ("userstepping");

	dt = flags.GetNumFlag ("dt", 0.001);
	tend = flags.GetNumFlag ("tend", 1);

	aneg = flags.GetNumFlag ("aneg", 1.0);
	apos = flags.GetNumFlag ("apos", 1.0);
	bneg = flags.GetNumFlag ("bneg", 1.0);
	bpos = flags.GetNumFlag ("bpos", 1.0);
	lambda = flags.GetNumFlag ("lambda", 10.0);
	

	coef_bconvneg = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_conv_neg", "bconvneg"));
	coef_bconvpos = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_conv_pos", "bconvpos"));

	coef_brhsneg = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_rhs_neg", "brhsneg"));
	coef_brhspos = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_rhs_pos", "brhspos"));

	coef_binineg = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_ini_neg", "binineg"));
	coef_binipos = pde.GetCoefficientFunction (flags.GetStringFlag ("beta_ini_pos", "binipos"));

	coef_bndneg = pde.GetCoefficientFunction (flags.GetStringFlag ("boundary_neg", "bndneg"));
	coef_bndpos = pde.GetCoefficientFunction (flags.GetStringFlag ("boundary_pos", "bndpos"));

	coef_aneg = new ConstantCoefficientFunction(aneg);
	coef_apos = new ConstantCoefficientFunction(apos);
	coef_bneg = new ConstantCoefficientFunction(bneg);
	coef_bpos = new ConstantCoefficientFunction(bpos);
	coef_abneg = new ConstantCoefficientFunction(aneg*bneg);
	coef_abpos = new ConstantCoefficientFunction(apos*bpos);

	coef_lambda = new ConstantCoefficientFunction(lambda);

	coef_zero = new ConstantCoefficientFunction(0.0);
	coef_one = new ConstantCoefficientFunction(1.0);

	coef_told = new ConstantCoefficientFunction(10.0);
	coef_tnew = new ConstantCoefficientFunction(10.0+dt);

  }


  virtual ~NumProcSolveInstatX() 
  { 
	delete coef_lambda;
	delete coef_zero;
	delete coef_one;
	delete coef_told;
	delete coef_tnew;

	delete coef_aneg;
	delete coef_bneg;
	delete coef_abneg;
	delete coef_apos;
	delete coef_bpos;
	delete coef_abpos;
  }


  /*
	creates an solver object
  */
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
	return new NumProcSolveInstatX (pde, flags);
  }


  /*
	solve at one level
  */
  virtual void Do(LocalHeap & lh)
  {

	
	// cout << "TESTING" << endl;

	// Array <EvalFunction*> evals(1);
	
	// std::cout << " a " << std::endl;
	// evals[0] = new EvalFunction("x+y+z");
	// std::cout << " b " << std::endl;
	// // SpaceTimeDomainVariableCoefficientFunction<D> coef_test( evals );
	// DomainVariableCoefficientFunction<D+1> coef_test( evals );
	// std::cout << " c " << std::endl;
	// ElementTransformation & eltrans = pde.GetMeshAccess().GetTrafo (0, VOL, lh);
	// IntegrationPoint ip (0.0);
	// MappedIntegrationPoint<D,D> mip(ip,eltrans);
	// DimMappedIntegrationPoint<D+1> mip2(ip,eltrans);
	// mip2.Point().Range(0,D) = mip.GetPoint();
	// mip2.Point()[D] = 1.0;
	// Vec<1> test;
	// coef_test.Evaluate(mip2,test);	
	// std::cout << " mip2.GetPoint() = " << mip2.GetPoint() << std::endl;
	// cout << " test: " << test  << endl;
	
	cout << "solve solveinstatx pde" << endl;

	// PointContainer<2+1> pc;
	// Array<Simplex<2+1> *> ret(0);
	// DecomposeIntoSimplices<2,3>(ET_TRIG,ET_SEGM,ret,pc,lh);

	// PointContainer<3+1> pc2;
	// Array<Simplex<3+1> *> ret2(0);
	// DecomposeIntoSimplices<3,4>(ET_TET,ET_SEGM,ret2,pc2,lh);

	// PointContainer<2> pc3;
	// Array<Simplex<2> *> ret3(0);
	// DecomposeIntoSimplices<2,2>(ET_TRIG,ET_POINT,ret3,pc3,lh);

	// PointContainer<3> pc4;
	// Array<Simplex<3> *> ret4(0);
	// DecomposeIntoSimplices<3,3>(ET_TET,ET_POINT,ret4,pc4,lh);

	BilinearForm * bftau;
	Flags massflags;
	massflags.SetFlag ("fespace", fesstr.c_str());
	bftau = CreateBilinearForm (fes, "bftau", massflags);
	bftau -> SetUnusedDiag (0);

	Array<CoefficientFunction*> coefs_timeder(2);
	coefs_timeder[0] = coef_bneg;
	coefs_timeder[1] = coef_bpos;

	SpaceTimeXTimeDerivativeIntegrator<D> bfidt(coefs_timeder);
	bftau -> AddIntegrator (&bfidt);

	Array<CoefficientFunction*> coefs_xlaplace(4);
	coefs_xlaplace[0] = coef_abneg;
	coefs_xlaplace[1] = coef_abpos;
	coefs_xlaplace[2] = coef_told;
	coefs_xlaplace[3] = coef_tnew;

	SpaceTimeXLaplaceIntegrator<D> bfilap(coefs_xlaplace);
	bftau -> AddIntegrator (&bfilap);

	Array<CoefficientFunction *> coefs_xnitsche(7);
	coefs_xnitsche[0] = coef_aneg;
	coefs_xnitsche[1] = coef_apos;
	coefs_xnitsche[2] = coef_bneg;
	coefs_xnitsche[3] = coef_bpos;
	coefs_xnitsche[4] = coef_lambda;
	coefs_xnitsche[5] = coef_told;
	coefs_xnitsche[6] = coef_tnew;

	SpaceTimeXNitscheIntegrator<D,NITSCHE_VARIANTS::HANSBOBETA> bfixnit (coefs_xnitsche);
	bftau -> AddIntegrator (&bfixnit);

	Array<CoefficientFunction*> coefs_xconvection(4);
	coefs_xconvection[0] = coef_bconvneg;
	coefs_xconvection[1] = coef_bconvpos;
	coefs_xconvection[2] = coef_told;
	coefs_xconvection[3] = coef_tnew;

	SpaceTimeXConvectionIntegrator<D> bfixconv(coefs_xconvection);
	bftau -> AddIntegrator (&bfixconv);


	Array<CoefficientFunction*> coefs_timetr(2);
	coefs_timetr[0] = coef_bneg;
	coefs_timetr[1] = coef_bpos;

	SpaceTimeXTraceMassIntegrator<D,PAST> bfitimetr(coefs_timetr);
	bftau -> AddIntegrator (&bfitimetr);

	//DifferentialOperator * traceop = new SpaceTimeTimeTraceIntegrator<D,FUTURE>(new ConstantCoefficientFunction(1.0));

	DifferentialOperator * uptrace = new T_DifferentialOperator<DiffOpTimeTrace<D,FUTURE> >();
	GridFunctionCoefficientFunction coef_u_neg (*(gfu_vis->GetComponent(0)), uptrace); //, traceop);
	GridFunctionCoefficientFunction coef_u_pos (*(gfu_vis->GetComponent(1)), uptrace); //, traceop);

	LinearForm * lfrhs;
	Flags massflags2;
	// massflags2.SetFlag ("fespace", fesstr.c_str());
	// lfrhs = pde.AddLinearForm ("lfrhs", massflags2);

	FESpace * fes = const_cast<FESpace*>(&(gfu->GetFESpace()));
	CompoundFESpace & compfes = *dynamic_cast<CompoundFESpace * >(fes);
	XFESpace<D,D+1> & xfes = *dynamic_cast<XFESpace<D,D+1> * >(compfes[1]);
	
	CompoundFESpace & compfes2 = *dynamic_cast<CompoundFESpace * >(fesvis);
	LevelsetContainerFESpace & lcfes = *dynamic_cast<LevelsetContainerFESpace * >(compfes2[2]);

	lfrhs = CreateLinearForm(fes,"lfrhs",massflags2);

	
	Array<CoefficientFunction *> coef_upw(2);
	coef_upw[0] = coef_binineg;
	coef_upw[1] = coef_binipos;
	SpaceTimeXTraceSourceIntegrator<D,PAST> lfi_tr(coef_upw);
    lfrhs -> AddIntegrator (&lfi_tr);

	// just for testing
	Array<CoefficientFunction *> coef_rhs(4);
	coef_rhs[0] = coef_brhsneg;
	coef_rhs[1] = coef_brhspos;
	coef_rhs[2] = coef_told;
	coef_rhs[3] = coef_tnew;
	SpaceTimeXSourceIntegrator<D> lfi_rhs(coef_rhs);

	lfrhs -> AddIntegrator (&lfi_rhs);

	Array<CoefficientFunction* > bndcoefs(2);
	bndcoefs[0] = coef_bndneg;
	bndcoefs[1] = coef_bndpos;

	// time stepping
	double t;
	for (t = 0; t < tend; t += dt)
	{
	  HeapReset hr(lh);
	  TimeInterval ti(t,t+dt);
	  xfes.SetTimeInterval(ti);
	  lcfes.SetTime(ti.first,ti.second);
	  fes->Update(lh);
	  gfu->Update();

	  SetValuesX<D,double>( bndcoefs, ti, *gfu, true, lh);

	  BaseVector & vecu = gfu->GetVector();

	  bfilap.SetTimeInterval(ti);
	  bfixnit.SetTimeInterval(ti);
	  bfixconv.SetTimeInterval(ti);

	  bftau -> ReAssemble(lh,true);

	  BaseMatrix & mata = bftau->GetMatrix();
	  dynamic_cast<BaseSparseMatrix&> (mata) . SetInverseType (inversetype);
	  BaseMatrix & invmat = * dynamic_cast<BaseSparseMatrix&> (mata) . InverseMatrix(gfu->GetFESpace().GetFreeDofs());

	  lfi_rhs.SetTimeInterval(ti);
	  lfrhs -> Assemble(lh);

	  const BaseVector * vecf = &(lfrhs->GetVector());
	  
	  // vecu = 0.0;
	  BaseVector & d = *vecu.CreateVector();
	  BaseVector & w = *vecu.CreateVector();
	  
	  d = *vecf - mata * vecu;
	  w = invmat * d;
	  vecu += w;

	  // update status text
	  cout << "\rt = " << t;
	  cout << flush;
	  // update visualization
	  delete &invmat;
	  delete &d;
	  delete &w;

	  xfes.XToNegPos(*gfu,*gfu_vis);

	  lfi_tr.ChangeNegPosCoefficient(&coef_u_neg, &coef_u_pos, bneg, bpos);
	  Ng_Redraw ();
	  
	  if (userstepping)
	  	getchar();
	}
	cout << "\r               \rt = " << tend;
	cout << endl;
	// delete &d;
	// delete &w;
  }


  /*
	GetClassName
  */
  virtual string GetClassName () const
  {
	return "SolveInstatX Solver";
  }


  /*
	PrintReport
  */
  virtual void PrintReport (ostream & ost)
  {
	ost << GetClassName() << endl
		<< "NOT UP TO DATE\n " << endl;
		// << "Linear-form     = " << lff->GetName() << endl
		// << "Gridfunction    = " << gfu->GetName() << endl
		// << "dt              = " << dt << endl
		// << "tend            = " << tend << endl;
  }


  /*
	PrintDoc
  */
  static void PrintDoc (ostream & ost)
  {
	ost << 
	  "\n\n NOT UP TO DATE!!!\nNumproc SolveInstat:\n" \
	  "------------------\n" \
	  "Solves a solveinstat partial differential equation by an implicite Euler method\n\n" \
	  "Required flags:\n" 
	  "-linearform=<lfname>\n" \
	  "    linear-form providing the right hand side\n" \
	  "-gridfunction=<gfname>\n" \
	  "    grid-function to store the solution vector\n" 
	  "-dt=<value>\n"
	  "    time step\n"
	  "-tend=<value>\n"
	  "    total time\n"
		<< endl;
  }
};


static RegisterNumProc<NumProcSolveInstatX<2> > npst_solveinstat2dx("stx_solveinstat");
static RegisterNumProc<NumProcSolveInstatX<3> > npst_solveinstat3dx("stx_solveinstat");
