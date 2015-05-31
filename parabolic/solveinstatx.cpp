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
using namespace ngcomp;
// #include "../utils/calccond.hpp"


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
  
  shared_ptr<FESpace> fes;
  shared_ptr<FESpace> fesvis;

  shared_ptr<GridFunction> gfu;

  shared_ptr<GridFunction> gfu_vis = NULL;

  shared_ptr<BilinearForm> bftau;
  shared_ptr<LinearForm> lfrhs;

  // time step
  double dt;
  // total time
  double tstart;
  double tend;
  // direct solver type
  string inversetype;
  string fesstr;
	
  bool userstepping;
  // bool calccond;
  bool ghostpenalty;
  bool directsolve;
  bool minimal_stabilization;

  double sleep_time;

  double bneg = 1.0;
  double bpos = 1.0;
  double aneg = 1.0;
  double apos = 1.0;
  double lambda = 10.0;
  double delta = 0.1;

  shared_ptr<CoefficientFunction> coef_bconvneg = NULL;
  shared_ptr<CoefficientFunction> coef_bconvpos = NULL;

  shared_ptr<CoefficientFunction> coef_bneg = NULL;
  shared_ptr<CoefficientFunction> coef_bpos = NULL;
  shared_ptr<CoefficientFunction> coef_aneg = NULL;
  shared_ptr<CoefficientFunction> coef_apos = NULL;
  shared_ptr<CoefficientFunction> coef_abneg = NULL;
  shared_ptr<CoefficientFunction> coef_abpos = NULL;

  shared_ptr<CoefficientFunction> coef_lambda = NULL; // Nitsche
  shared_ptr<CoefficientFunction> coef_delta = NULL; // ghost penalty

  shared_ptr<CoefficientFunction> coef_zero = NULL;
  shared_ptr<CoefficientFunction> coef_one = NULL;

  shared_ptr<CoefficientFunction> coef_told = NULL;
  shared_ptr<CoefficientFunction> coef_tnew = NULL;

  shared_ptr<CoefficientFunction> coef_upwneg = NULL;
  shared_ptr<CoefficientFunction> coef_upwpos = NULL;

  shared_ptr<CoefficientFunction> coef_brhsneg = NULL;
  shared_ptr<CoefficientFunction> coef_brhspos = NULL;

  shared_ptr<CoefficientFunction> coef_binineg = NULL;
  shared_ptr<CoefficientFunction> coef_binipos = NULL;

  shared_ptr<CoefficientFunction> coef_bndneg = NULL;
  shared_ptr<CoefficientFunction> coef_bndpos = NULL;

  ErrorTable errtab;
  SolutionCoefficients<D> solcoef;

  shared_ptr<SpaceTimeXTimeDerivativeIntegrator<D> > bfidt;
  shared_ptr<SpaceTimeXLaplaceIntegrator<D> > bfilap; 
  shared_ptr<GhostPenaltyIntegrator<D,1> > bfigho;
  shared_ptr<SpaceTimeXNitscheIntegrator<D,NITSCHE_VARIANTS::HANSBO> > bfixnit;
  shared_ptr<SpaceTimeXConvectionIntegrator<D> > bfixconv;
  shared_ptr<SpaceTimeXTraceMassIntegrator<D,PAST> > bfitimetr;

  shared_ptr<SpaceTimeXTraceSourceIntegrator<D,PAST> > lfi_tr;
  shared_ptr<SpaceTimeXSourceIntegrator<D> > lfi_rhs;

  shared_ptr<Preconditioner> localprec;


public:
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcSolveInstatX (shared_ptr<PDE> apde, const Flags & flags) : NumProc (apde), errtab(), solcoef(apde, flags)
  {
    // in the input-file, you specify the bilinear-forms for the stiffness and for the mass-term
    // like  "-bilinearforma=k". Default arguments are 'a' and 'm'

    // bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    // bfm = pde.GetBilinearForm (flags.GetStringFlag ("bilinearformm", "m"));
    // lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));

    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    gfu_vis = apde->GetGridFunction (flags.GetStringFlag ("gf_vis", "u_vis"),true);

    inversetype = flags.GetStringFlag ("solver", "pardiso");
    fesstr = flags.GetStringFlag ("fespace", "fes_st");
    fes = apde->GetFESpace (fesstr.c_str());

    string fesvisstr = flags.GetStringFlag ("fespacevis", "fes_negpos");
    fesvis = apde->GetFESpace (fesvisstr.c_str());

    userstepping = flags.GetDefineFlag ("userstepping");
    // calccond = flags.GetDefineFlag ("calccond");
    ghostpenalty = flags.GetDefineFlag ("ghostpenalty");
    directsolve = flags.GetDefineFlag ("direct");
    minimal_stabilization = flags.GetDefineFlag ("minimal_stabilization");

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

    coef_bconvneg = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_conv_neg", "bconvneg"));
    coef_bconvpos = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_conv_pos", "bconvpos"));

    coef_brhsneg = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_rhs_neg", "brhsneg"));
    coef_brhspos = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_rhs_pos", "brhspos"));

    coef_binineg = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_ini_neg", "binineg"));
    coef_binipos = apde->GetCoefficientFunction (flags.GetStringFlag ("beta_ini_pos", "binipos"));

    coef_bndneg = apde->GetCoefficientFunction (flags.GetStringFlag ("boundary_neg", "bndneg"));
    coef_bndpos = apde->GetCoefficientFunction (flags.GetStringFlag ("boundary_pos", "bndpos"));

    coef_aneg = make_shared<ConstantCoefficientFunction>(aneg/bneg);
    coef_apos = make_shared<ConstantCoefficientFunction>(apos/bpos);
    coef_bneg = make_shared<ConstantCoefficientFunction>(1.0/bneg);
    coef_bpos = make_shared<ConstantCoefficientFunction>(1.0/bpos);
    coef_abneg = make_shared<ConstantCoefficientFunction>(aneg/bneg);
    coef_abpos = make_shared<ConstantCoefficientFunction>(apos/bpos);

    coef_lambda = make_shared<ConstantCoefficientFunction>(lambda);
    coef_delta = make_shared<ConstantCoefficientFunction>(delta);

    coef_zero = make_shared<ConstantCoefficientFunction>(0.0);
    coef_one = make_shared<ConstantCoefficientFunction>(1.0);

    coef_told = make_shared<ConstantCoefficientFunction>(0.0);
    coef_tnew = make_shared<ConstantCoefficientFunction>(0.0+dt);

    Flags blflags;
    blflags.SetFlag ("fespace", fesstr.c_str());
    bftau = apde->AddBilinearForm("bftau", blflags);
    // bftau = CreateBilinearForm (fes, "bftau", blflags);
    bftau -> SetUnusedDiag (0);

    shared_ptr<FESpace> fes = gfu->GetFESpace();
    Flags lflags;
    lfrhs = CreateLinearForm(fes,"lfrhs",lflags);

    Flags empty;
    empty.SetFlag ("type", "spacetime");
    // empty.SetFlag ("type", "local");
    // empty.SetFlag ("type", "bddc");
    empty.SetFlag ("block");
    empty.SetFlag ("fespace", fesstr.c_str());
    empty.SetFlag ("bilinearform", "bftau");
    empty.SetFlag ("laterupdate");
    // localprec = new LocalPreconditioner (&pde, empty, "l");
    // localprec = new BDDC (&pde, empty, "l");
    localprec = apde->AddPreconditioner("l",empty);
  }

  void DeleteBilinearFormIntegrators()
  {
    // if (bfidt) { delete bfidt; bfidt=NULL; }
    // if (bfilap) { delete bfilap;  bfilap=NULL; }
    // if (bfigho) { delete bfigho; bfigho=NULL; }
    // if (bfixnit) { delete bfixnit; bfixnit=NULL; }
    // if (bfixconv) { delete bfixconv; bfixconv=NULL; }
    // if (bfitimetr) { delete bfitimetr; bfitimetr=NULL; }
  }

  void DeleteLinearFormIntegrators()
  {
    // if (lfi_tr) { delete lfi_tr; lfi_tr=NULL; }
    // if (lfi_rhs) { delete lfi_rhs; lfi_rhs=NULL; }
  }

  virtual ~NumProcSolveInstatX() 
  { 
    DeleteLinearFormIntegrators();
    DeleteBilinearFormIntegrators();

  //   delete coef_delta;
  //   delete coef_lambda;
  //   delete coef_zero;
  //   delete coef_one;
  //   delete coef_told;
  //   delete coef_tnew;

  //   delete coef_aneg;
  //   delete coef_bneg;
  //   delete coef_abneg;
  //   delete coef_apos;
  //   delete coef_bpos;
  //   delete coef_abpos;
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
    Array<shared_ptr<CoefficientFunction>> coefs_timeder(2);
    coefs_timeder[0] = coef_bneg;
    coefs_timeder[1] = coef_bpos;

    bfidt = make_shared<SpaceTimeXTimeDerivativeIntegrator<D>> (coefs_timeder);
    bftau -> AddIntegrator (bfidt);

    Array<shared_ptr<CoefficientFunction>> coefs_xlaplace(4);
    coefs_xlaplace[0] = coef_abneg;
    coefs_xlaplace[1] = coef_abpos;
    coefs_xlaplace[2] = coef_told;
    coefs_xlaplace[3] = coef_tnew;

    bfilap = make_shared<SpaceTimeXLaplaceIntegrator<D>> (coefs_xlaplace);
    bftau -> AddIntegrator (bfilap);

    Array<shared_ptr<CoefficientFunction>> coefs_ghostpen(5);
    coefs_ghostpen[0] = coef_abneg;
    coefs_ghostpen[1] = coef_abpos;
    coefs_ghostpen[2] = coef_told;
    coefs_ghostpen[3] = coef_tnew;
    coefs_ghostpen[4] = coef_delta;

    bfigho = make_shared<GhostPenaltyIntegrator<D,1>>(coefs_ghostpen);
    if (ghostpenalty)
      bftau -> AddIntegrator (bfigho);

    Array<shared_ptr<CoefficientFunction>> coefs_xnitsche(minimal_stabilization ? 6 : 7);
    coefs_xnitsche[0] = coef_aneg;
    coefs_xnitsche[1] = coef_apos;
    coefs_xnitsche[2] = coef_one;
    coefs_xnitsche[3] = coef_one;
    if (minimal_stabilization)
      {
        coefs_xnitsche[4] = coef_told;
        coefs_xnitsche[5] = coef_tnew;
      }
    else
      {
        coefs_xnitsche[4] = coef_lambda;
        coefs_xnitsche[5] = coef_told;
        coefs_xnitsche[6] = coef_tnew;
      }

    bfixnit = make_shared<SpaceTimeXNitscheIntegrator<D,NITSCHE_VARIANTS::HANSBO> > (coefs_xnitsche);
    bftau -> AddIntegrator (bfixnit);

    Array<shared_ptr<CoefficientFunction>> coefs_xconvection(4);
    coefs_xconvection[0] = coef_bconvneg;
    coefs_xconvection[1] = coef_bconvpos;
    coefs_xconvection[2] = coef_told;
    coefs_xconvection[3] = coef_tnew;

    bfixconv = make_shared<SpaceTimeXConvectionIntegrator<D> > (coefs_xconvection);
    bftau -> AddIntegrator (bfixconv);

    Array<shared_ptr<CoefficientFunction>> coefs_timetr(2);
    coefs_timetr[0] = coef_bneg;
    coefs_timetr[1] = coef_bpos;

    bfitimetr = make_shared<SpaceTimeXTraceMassIntegrator<D,PAST> > (coefs_timetr);
    bftau -> AddIntegrator (bfitimetr);
  }

  void AddLinearFormIntegrators()
  {
    Array<shared_ptr<CoefficientFunction>> coef_upw(2);
    coef_upw[0] = coef_binineg;
    coef_upw[1] = coef_binipos;
    lfi_tr = make_shared<SpaceTimeXTraceSourceIntegrator<D,PAST>> (coef_upw);
    lfrhs -> AddIntegrator (lfi_tr);

    // just for testing
    Array<shared_ptr<CoefficientFunction>> coef_rhs(4);
    coef_rhs[0] = coef_brhsneg;
    coef_rhs[1] = coef_brhspos;
    coef_rhs[2] = coef_told;
    coef_rhs[3] = coef_tnew;
    lfi_rhs = make_shared<SpaceTimeXSourceIntegrator<D>> (coef_rhs);
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

    // if (refinements>1)
    //   dt /= 2.0;


    shared_ptr<DifferentialOperator> uptrace = make_shared<T_DifferentialOperator<DiffOpTimeTrace<D,FUTURE> > >();
    shared_ptr<GridFunctionCoefficientFunction> coef_u_neg = make_shared<GridFunctionCoefficientFunction> (gfu_vis->GetComponent(0), uptrace); //, traceop);
    shared_ptr<GridFunctionCoefficientFunction> coef_u_pos = make_shared<GridFunctionCoefficientFunction> (gfu_vis->GetComponent(1), uptrace); //, traceop);

    shared_ptr<FESpace> fes = gfu->GetFESpace();
    shared_ptr<CompoundFESpace> compfes = dynamic_pointer_cast<CompoundFESpace>(fes);
    shared_ptr<T_XFESpace<D,D+1> > xfes = dynamic_pointer_cast<T_XFESpace<D,D+1> >((*compfes)[1]);
	
    shared_ptr<CompoundFESpace> compfes2 = dynamic_pointer_cast<CompoundFESpace>(fesvis);
    shared_ptr<LevelsetContainerFESpace> lcfes = dynamic_pointer_cast<LevelsetContainerFESpace>((*compfes2)[2]);

    static ConstantCoefficientFunction copy1(1.0);
    static ConstantCoefficientFunction copy2(1.0);

    if (refinements==1)
      AddLinearFormIntegrators();
    else
      {
        lfi_tr->ChangeNegPosCoefficient(coef_binineg, coef_binipos, 1.0, 1.0);
      }

    std::cout << " coef_binineg = " << coef_binineg << std::endl;
    std::cout << " coef_binipos = " << coef_binipos << std::endl;

    Array<shared_ptr<CoefficientFunction> > bndcoefs(2);
    bndcoefs[0] = coef_bndneg;
    bndcoefs[1] = coef_bndpos;

    std::ofstream outf("condition.out");
    // time stepping
    double t;
    for (t = tstart; t < tend-dt+1e-12; t += dt)
      {
        HeapReset hr(lh);
        TimeInterval ti(t,t+dt);
        xfes->SetTimeInterval(ti);
        lcfes->SetTime(ti.first,ti.second);
        fes->Update(lh);
        gfu->Update();

        SetValuesX<D,double>( bndcoefs, ti, gfu, true, lh);
        
        BaseVector & vecu = gfu->GetVector();

        bfilap->SetTimeInterval(ti);
        if (ghostpenalty)
          bfigho->SetTimeInterval(ti);
        bfixnit->SetTimeInterval(ti);
        bfixconv->SetTimeInterval(ti);

        bftau -> ReAssemble(lh,true);

        BaseMatrix & mata = bftau->GetMatrix();
        dynamic_cast<BaseSparseMatrix&> (mata) . SetInverseType (inversetype);
        shared_ptr<BaseMatrix> directinvmat = NULL;
        GMRESSolver<double> * itinvmat = NULL;

        if (directsolve) // || calccond)
          directinvmat = dynamic_cast<BaseSparseMatrix&> (mata) . InverseMatrix(gfu->GetFESpace()->GetFreeDofs());

        if (!directsolve)
        {
            static Timer timer ("localprec->Update()");
            RegionTimer reg (timer);

            localprec->Update();
            itinvmat = new GMRESSolver<double> (mata, *localprec);
            itinvmat->SetPrintRates(true);
            itinvmat->SetMaxSteps(10000);
            itinvmat->SetPrecision(1e-6);
        }

        BaseMatrix & invmat = directsolve ? *directinvmat : *itinvmat;
        
        lfi_tr->SetTime(ti.first);  //absolute time
        lfi_rhs->SetTimeInterval(ti);
        lfrhs -> Assemble(lh.Borrow());

        const BaseVector * vecf = &(lfrhs->GetVector());

        AutoVector d = vecu.CreateVector();
        AutoVector w = vecu.CreateVector();

        d = *vecf - mata * vecu;
        w = invmat * d;
        vecu += w;

        // update status text
        cout << "\r          \rt = " << std::setw(6) << ti.first << " to t = " << std::setw(6) << ti.second;
        if (!directsolve)
          cout << " - number of its.: " << std::setw(4) << itinvmat->GetSteps();
        cout << flush;

        // if (calccond)
        // {
        //   static Timer timer ("CalcCond");
        //   RegionTimer reg (timer);
        //   outf << t << "\t";
        //   CalcCond(mata,*directinvmat,gfu->GetFESpace()->GetFreeDofs(), false, false, &outf);
        //   CalcCond(mata,*directinvmat,gfu->GetFESpace()->GetFreeDofs(), false, true, &outf);
        //   if (!directsolve)
        //     outf << itinvmat->GetSteps() << "\t";
        //   outf << endl;
        // }
        // if (!directsolve)
        // {
        // 	cout << " ... " << endl;
        // 	localprec->Test();
        // 	cout << " ... " << endl;
        // }

        // update visualization
        // delete &d;
        // delete &w;
        // if (calccond || directsolve)
        //   delete directinvmat;

        if (!directsolve)
          delete itinvmat;

        // *testout << " t = " << t << " \n vecu = \n " << vecu << endl;
        {
          static Timer timer ("CalcXError");
          Flags emptyflags;
          RegionTimer reg (timer);
          if (abs(ti.second - tend) < 1e-6*dt)
              CalcXError<D>(gfu, NULL, solcoef, 4, aneg, apos, bneg, bpos, ti.second, errtab, lh, true, emptyflags);
          else
              CalcXError<D>(gfu, NULL, solcoef, 4, aneg, apos, bneg, bpos, ti.second, errtab, lh, false, emptyflags);
        }

        xfes->XToNegPos(gfu,gfu_vis);

        lfi_tr->ChangeNegPosCoefficient(coef_u_neg, coef_u_pos, 1.0/bneg, 1.0/bpos);
        Ng_Redraw ();
	  
        if (userstepping)
          getchar();
        if (sleep_time>0)
          usleep(sleep_time*1000);
        cout << endl;

      }
    cout << "\r          \rt = " << tend;
    cout << endl;

    // CalcXError<D>(gfu, NULL, solcoef, 4, bneg, bpos, tend, errtab, lh, true);

    // delete &d;
    // delete &w;
    bftau -> SetNonAssemble(true);
    // delete uptrace;
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
