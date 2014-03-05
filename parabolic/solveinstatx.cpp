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
#include "../xfem/ghostpenalty.hpp"
#include "../utils/error.hpp"

#include <diffop_impl.hpp>

#include <unistd.h>

// #include "../utils/stcoeff.hpp"

using namespace ngsolve;
#include "../utils/calccond.hpp"


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

  BilinearForm * bftau;
  LinearForm * lfrhs;

  // time step
  double dt;
  // total time
  double tstart;
  double tend;
  // direct solver type
  string inversetype;
  string fesstr;
	
  bool userstepping;
  bool calccond;
  bool ghostpenalty;
  double sleep_time;

  double bneg = 1.0;
  double bpos = 1.0;
  double aneg = 1.0;
  double apos = 1.0;
  double lambda = 10.0;
  double delta = 0.1;

  CoefficientFunction* coef_bconvneg = NULL;
  CoefficientFunction* coef_bconvpos = NULL;

  CoefficientFunction* coef_bneg = NULL;
  CoefficientFunction* coef_bpos = NULL;
  CoefficientFunction* coef_aneg = NULL;
  CoefficientFunction* coef_apos = NULL;
  CoefficientFunction* coef_abneg = NULL;
  CoefficientFunction* coef_abpos = NULL;

  CoefficientFunction* coef_lambda = NULL; // Nitsche
  CoefficientFunction* coef_delta = NULL; // ghost penalty

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

  ErrorTable errtab;
  SolutionCoefficients<D> solcoef;

  SpaceTimeXTimeDerivativeIntegrator<D> * bfidt;
  SpaceTimeXLaplaceIntegrator<D> * bfilap; 
  LowOrderGhostPenaltyIntegrator<D> * bfigho;
  SpaceTimeXNitscheIntegrator<D,NITSCHE_VARIANTS::HANSBO> * bfixnit;
  SpaceTimeXConvectionIntegrator<D> * bfixconv;
  SpaceTimeXTraceMassIntegrator<D,PAST> * bfitimetr;

  SpaceTimeXTraceSourceIntegrator<D,PAST> * lfi_tr;
  SpaceTimeXSourceIntegrator<D> * lfi_rhs;

  Preconditioner * localprec;


public:
  /*
	In the constructor, the solver class gets the flags from the pde - input file.
	the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcSolveInstatX (PDE & apde, const Flags & flags) : NumProc (apde), errtab(), solcoef(apde, flags)
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
	calccond = flags.GetDefineFlag ("calccond");
	ghostpenalty = flags.GetDefineFlag ("ghostpenalty");

	dt = flags.GetNumFlag ("dt", 0.001);
	tstart = flags.GetNumFlag ("tstart", 0.0);
	tend = flags.GetNumFlag ("tend", 1);

	aneg = flags.GetNumFlag ("aneg", 1.0);
	apos = flags.GetNumFlag ("apos", 1.0);
	bneg = flags.GetNumFlag ("bneg", 1.0);
	bpos = flags.GetNumFlag ("bpos", 1.0);
	lambda = flags.GetNumFlag ("lambda", 10.0);
	delta = flags.GetNumFlag ("delta", 0.1);
	
	sleep_time = flags.GetNumFlag ("pause_after_step", 0.0);

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
	coef_delta = new ConstantCoefficientFunction(delta);

	coef_zero = new ConstantCoefficientFunction(0.0);
	coef_one = new ConstantCoefficientFunction(1.0);

	coef_told = new ConstantCoefficientFunction(0.0);
	coef_tnew = new ConstantCoefficientFunction(0.0+dt);

	Flags blflags;
	blflags.SetFlag ("fespace", fesstr.c_str());
	bftau = pde.AddBilinearForm("bftau", blflags);
	// bftau = CreateBilinearForm (fes, "bftau", blflags);
	bftau -> SetUnusedDiag (0);

	FESpace * fes = const_cast<FESpace*>(&(gfu->GetFESpace()));
	Flags lflags;
	lfrhs = CreateLinearForm(fes,"lfrhs",lflags);

	Flags empty;
	empty.SetFlag ("type", "local");
	empty.SetFlag ("fespace", fesstr.c_str());
	empty.SetFlag ("bilinearform", "bftau");
	empty.SetFlag ("laterupdate");
	// localprec = new LocalPreconditioner (&pde, empty, "l");
	// localprec = new BDDC (&pde, empty, "l");
	localprec = pde.AddPreconditioner("l",empty);
  }


  virtual ~NumProcSolveInstatX() 
  { 

	delete bfidt;
	delete bfilap; 
	delete bfigho;
	delete bfixnit;
	delete bfixconv;
	delete bfitimetr;

	delete lfi_tr;
	delete lfi_rhs;

	delete coef_delta;
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


  void AddBilinearFormIntegrators()
  {
	Array<CoefficientFunction*> coefs_timeder(2);
	coefs_timeder[0] = coef_bneg;
	coefs_timeder[1] = coef_bpos;

	bfidt = new SpaceTimeXTimeDerivativeIntegrator<D> (coefs_timeder);
	bftau -> AddIntegrator (bfidt);

	Array<CoefficientFunction*> coefs_xlaplace(4);
	coefs_xlaplace[0] = coef_abneg;
	coefs_xlaplace[1] = coef_abpos;
	coefs_xlaplace[2] = coef_told;
	coefs_xlaplace[3] = coef_tnew;

	bfilap = new SpaceTimeXLaplaceIntegrator<D> (coefs_xlaplace);
	bftau -> AddIntegrator (bfilap);

	Array<CoefficientFunction*> coefs_ghostpen(5);
	coefs_ghostpen[0] = coef_abneg;
	coefs_ghostpen[1] = coef_abpos;
	coefs_ghostpen[2] = coef_told;
	coefs_ghostpen[3] = coef_tnew;
	coefs_ghostpen[4] = coef_delta;

	bfigho = new LowOrderGhostPenaltyIntegrator<D> (coefs_ghostpen);
	if (ghostpenalty)
	  bftau -> AddIntegrator (bfigho);

	Array<CoefficientFunction *> coefs_xnitsche(7);
	coefs_xnitsche[0] = coef_aneg;
	coefs_xnitsche[1] = coef_apos;
	coefs_xnitsche[2] = coef_bneg;
	coefs_xnitsche[3] = coef_bpos;
	coefs_xnitsche[4] = coef_lambda;
	coefs_xnitsche[5] = coef_told;
	coefs_xnitsche[6] = coef_tnew;

	bfixnit = new SpaceTimeXNitscheIntegrator<D,NITSCHE_VARIANTS::HANSBO> (coefs_xnitsche);
	bftau -> AddIntegrator (bfixnit);

	Array<CoefficientFunction*> coefs_xconvection(4);
	coefs_xconvection[0] = coef_bconvneg;
	coefs_xconvection[1] = coef_bconvpos;
	coefs_xconvection[2] = coef_told;
	coefs_xconvection[3] = coef_tnew;

	bfixconv = new SpaceTimeXConvectionIntegrator<D> (coefs_xconvection);
	bftau -> AddIntegrator (bfixconv);

	Array<CoefficientFunction*> coefs_timetr(2);
	coefs_timetr[0] = coef_bneg;
	coefs_timetr[1] = coef_bpos;

	bfitimetr = new SpaceTimeXTraceMassIntegrator<D,PAST> (coefs_timetr);
	bftau -> AddIntegrator (bfitimetr);
  }

  void AddLinearFormIntegrators()
  {
	Array<CoefficientFunction *> coef_upw(2);
	coef_upw[0] = coef_binineg;
	coef_upw[1] = coef_binipos;
	lfi_tr = new SpaceTimeXTraceSourceIntegrator<D,PAST> (coef_upw);
    lfrhs -> AddIntegrator (lfi_tr);

	// just for testing
	Array<CoefficientFunction *> coef_rhs(4);
	coef_rhs[0] = coef_brhsneg;
	coef_rhs[1] = coef_brhspos;
	coef_rhs[2] = coef_told;
	coef_rhs[3] = coef_tnew;
	lfi_rhs = new SpaceTimeXSourceIntegrator<D> (coef_rhs);
	lfrhs -> AddIntegrator (lfi_rhs);
  }
  
  /*
	solve at one level
  */
  virtual void Do(LocalHeap & lh)
  {
	bftau -> SetNonAssemble(false);

    static int refinements = 0;
    cout << " This is the Do-call on refinement level " << refinements << std::endl;
    refinements++;
    errtab.Reset();

	cout << "solve instatx pde" << endl;

	if (refinements==1)
	  AddBilinearFormIntegrators();

	DifferentialOperator * uptrace = new T_DifferentialOperator<DiffOpTimeTrace<D,FUTURE> >();
	GridFunctionCoefficientFunction coef_u_neg (*(gfu_vis->GetComponent(0)), uptrace); //, traceop);
	GridFunctionCoefficientFunction coef_u_pos (*(gfu_vis->GetComponent(1)), uptrace); //, traceop);

	FESpace * fes = const_cast<FESpace*>(&(gfu->GetFESpace()));
	CompoundFESpace & compfes = *dynamic_cast<CompoundFESpace * >(fes);
	XFESpace<D,D+1> & xfes = *dynamic_cast<XFESpace<D,D+1> * >(compfes[1]);
	
	CompoundFESpace & compfes2 = *dynamic_cast<CompoundFESpace * >(fesvis);
	LevelsetContainerFESpace & lcfes = *dynamic_cast<LevelsetContainerFESpace * >(compfes2[2]);

	if (refinements==1)
	  AddLinearFormIntegrators();
	else
	  lfi_tr->ChangeNegPosCoefficient(coef_binineg, coef_binipos, bneg, bpos);

	Array<CoefficientFunction* > bndcoefs(2);
	bndcoefs[0] = coef_bndneg;
	bndcoefs[1] = coef_bndpos;

	// time stepping
	double t;
	for (t = tstart; t < tend; t += dt)
	{
	  HeapReset hr(lh);
	  TimeInterval ti(t,t+dt);
	  xfes.SetTimeInterval(ti);
	  lcfes.SetTime(ti.first,ti.second);
	  fes->Update(lh);
	  gfu->Update();

	  SetValuesX<D,double>( bndcoefs, ti, *gfu, true, lh);

	  BaseVector & vecu = gfu->GetVector();

	  bfilap->SetTimeInterval(ti);
	  if (ghostpenalty)
		bfigho->SetTimeInterval(ti);
	  bfixnit->SetTimeInterval(ti);
	  bfixconv->SetTimeInterval(ti);

	  bftau -> ReAssemble(lh,true);

	  BaseMatrix & mata = bftau->GetMatrix();
	  dynamic_cast<BaseSparseMatrix&> (mata) . SetInverseType (inversetype);
	  BaseMatrix * directinvmat = NULL;
	  if (calccond)
		directinvmat = dynamic_cast<BaseSparseMatrix&> (mata) . InverseMatrix(gfu->GetFESpace().GetFreeDofs());

	  localprec->Update();
	  GMRESSolver<double> invmat (mata, *localprec);
	  // invmat.SetPrintRates(true);
	  invmat.SetMaxSteps(10000);

	  lfi_tr->SetTime(ti.first);  //absolute time
	  lfi_rhs->SetTimeInterval(ti);
	  lfrhs -> Assemble(lh);

	  const BaseVector * vecf = &(lfrhs->GetVector());
	  
	  BaseVector & d = *vecu.CreateVector();
	  BaseVector & w = *vecu.CreateVector();
	  
	  d = *vecf - mata * vecu;
	  w = invmat * d;
	  vecu += w;

	  // update status text
	  cout << "\r          \rt = " << std::setw(6) << ti.first << " to t = " << std::setw(6) << ti.second;
	  cout << " - number of its.: " << std::setw(4) << invmat.GetSteps();
	  cout << flush;
	  
	  if (calccond)
		CalcCond(mata,*directinvmat,gfu->GetFESpace().GetFreeDofs(), false);

	  // update visualization
	  // delete &d;
	  delete &w;
	  if (calccond)
		directinvmat;

	  // *testout << " t = " << t << " \n vecu = \n " << vecu << endl;
      if (abs(ti.second - tend) < 1e-6*dt)
        CalcXError<D>(gfu, solcoef, 4, bneg, bpos, ti.second, errtab, lh, true);
      else
        CalcXError<D>(gfu, solcoef, 4, bneg, bpos, ti.second, errtab, lh, false);

	  xfes.XToNegPos(*gfu,*gfu_vis);

	  lfi_tr->ChangeNegPosCoefficient(&coef_u_neg, &coef_u_pos, bneg, bpos);
	  Ng_Redraw ();
	  
	  if (userstepping)
	  	getchar();
	  if (sleep_time>0)
	  	usleep(sleep_time*1000);
	  cout << endl;

	}
	cout << "\r          \rt = " << tend;
	cout << endl;

	// CalcXError<D>(gfu, solcoef, 4, bneg, bpos, tend, errtab, lh, true);

	// delete &d;
	// delete &w;
	bftau -> SetNonAssemble(true);
	delete uptrace;
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
