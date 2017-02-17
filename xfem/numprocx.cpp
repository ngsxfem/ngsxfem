/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
// #include "xintegration.hpp"
// #include "../spacetime/spacetimefespace.hpp"
#include "xfemIntegrators.hpp"
// #include "stxfemIntegrators.hpp"
#include "setvaluesx.hpp"
#include "../utils/error.hpp"
#include "../utils/output.hpp"
#include "../utils/calccond.hpp"
#include "../xfem/xFESpace.hpp"

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
namespace ngcomp
{ 
  template <int D, class SCAL>
  void SetValuesX (const Array<shared_ptr<CoefficientFunction>> & acoefs,
                   const TimeInterval & ti,
                   shared_ptr<GridFunction> bu,
                   bool bound,
                   LocalHeap & clh)
  {
    static Timer sv("timer setvaluesX"); RegionTimer r(sv);

    // S_GridFunction<SCAL> & u = dynamic_cast<S_GridFunction<SCAL> &> (bu);
    auto u = dynamic_pointer_cast<S_GridFunction<SCAL>> (bu);

    shared_ptr<FESpace> fes = u->GetFESpace();
    shared_ptr<MeshAccess> ma = fes->GetMeshAccess();
    
    ma->PushStatus ("setvalues");

    VorB vorb = VorB(bound);

    if (vorb == VOL)
      throw Exception ("not yet implemented");

    int dim   = fes->GetDimension();

    // const BilinearFormIntegrator & bli = *fes.GetIntegrator(bound);

    // if (&bli == NULL)
    //   throw Exception ("no integrator available");

    // int dimflux = 1; //diffop ? diffop->Dim() : bli.DimFlux(); 
    
    Array<int> cnti(fes->GetNDof());
    cnti = 0;

    u->GetVector() = 0.0;

    BilinearFormIntegrator * bfi = NULL;
    LinearFormIntegrator * lfi = NULL;

    shared_ptr<CompoundFESpace> cfes = dynamic_pointer_cast<CompoundFESpace>(fes);
    // shared_ptr<SpaceTimeFESpace> stfes = dynamic_pointer_cast<SpaceTimeFESpace>((*cfes)[0]);
    void * stfes = NULL;
    
    Array<shared_ptr<CoefficientFunction>> coefs(acoefs);
    Array<shared_ptr<CoefficientFunction>> coefs_one(2);
    auto one = make_shared<ConstantCoefficientFunction>(1.0);
    coefs_one[0] = one;
    coefs_one[1] = one;

    auto told = make_shared<ConstantCoefficientFunction>(ti.first);
    auto tnew = make_shared<ConstantCoefficientFunction>(ti.second);

    if (stfes != NULL)
    {
      throw Exception("Space time deactivated");
      // coefs.Append(told);
      // coefs.Append(tnew);
      // lfi = new (clh) SpaceTimeXNeumannIntegrator<D>(coefs);
      // coefs_one.Append(told);
      // coefs_one.Append(tnew);
      // bfi = new (clh) SpaceTimeXRobinIntegrator<D>(coefs_one);
    }
    else
    {
      lfi = new (clh) XNeumannIntegrator<D>(coefs);
      bfi = new (clh) XRobinIntegrator<D>(coefs_one);
    }

    ProgressOutput progress (ma, "setvaluesx element", ma->GetNE());

    IterateElements 
      (*fes, vorb, clh, 
       [&] (ElementId ei, LocalHeap & lh)
       {
         progress.Update ();

         // cout << "  ? on the boundary -> getchar(); " << endl; getchar();

         // if (bound && !fes->IsDirichletBoundary(ma->GetSElIndex(ei.Nr())))
         //   return;

         if (bound && !((*cfes)[0]->IsDirichletBoundary(ma->GetElIndex(ei))))
           return;

         const FiniteElement & bfel = fes->GetFE (ei, lh);

         const ElementTransformation & eltrans = ma->GetTrafo (ei, lh); 

         Array<int> dnums (bfel.GetNDof(), lh);
         fes->GetDofNrs (ei, dnums);

         FlatVector<SCAL> elvec(dnums.Size() * dim, lh);
         FlatVector<SCAL> elveci(dnums.Size() * dim, lh);
         FlatMatrix<double> elmat(dnums.Size(), lh);

         if (dim > 1)
           throw Exception ("dim > 1 not yet implemented");

         bfi->CalcElementMatrix (bfel, eltrans, elmat, lh);
         lfi->CalcElementVector (bfel, eltrans, elvec, lh);

         // for (int i = 0; i < elmat.Width(); ++i)
         //   elmat(i,i) += 1e-8;

         fes->TransformMat (ei, elmat, TRANSFORM_MAT_LEFT_RIGHT);
         fes->TransformVec (ei, elvec, TRANSFORM_RHS);

         bool regularize = false;
         for (int i = 0; i < dnums.Size(); ++i)
           if (abs(elmat(i,i)) < 1e-13)
             regularize = true;

         if (regularize)
         {
           for (int i = 0; i < dnums.Size(); ++i)
             for (int j = dnums.Size()/2; j < dnums.Size(); ++j)
             {
               elmat(i,j) = 0.0;
               elmat(j,i) = 0.0;
             }
           for (int j = dnums.Size()/2; j < dnums.Size(); ++j)
             elmat(j,j) = 1.0;
         }
              
         if (dnums.Size() < 50)
         {
           FlatCholeskyFactors<double> invelmat(elmat, lh);
           invelmat.Mult (elvec, elveci);
         }
         else
         {
           LapackInverse (elmat);
           elveci = elmat * elvec;
         }

         if (regularize)
         {
           const CompoundFiniteElement & cfel = 
             dynamic_cast<const CompoundFiniteElement&> (bfel);

           const XFiniteElement * xfe = NULL;

           if (xfe==NULL)
             xfe = dynamic_cast<const XFiniteElement* >(&cfel[1]);
           
           const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
           // const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();

           // DOMAIN_TYPE dt = xgeom.kappa[NEG] < 1e-12 ? POS : NEG;
           std::cout << " xgeom.kappa[NEG] = " << xgeom.kappa[NEG] << std::endl;
           std::cout << " xgeom.kappa[POS] = " << xgeom.kappa[POS] << std::endl;
           // for (int i = 0; i < dnums.Size()/2; ++i)
           // {
           //   if(xsign[i] == dt)
           //   {
           //     cnti[dnums[i]]--;
           //     elveci(dnums.Size()/2+i) = elveci(i);
           //     elveci(i) = 0.0;
           //   }
           //   else
           //     elveci(dnums.Size()/2+i) = 0.0;
           // }
         }

         bool printmore = regularize;
         // for (int i = 0; i < elveci.Size(); ++i)
         //   if (abs(elveci(i)) > 1.21)
         //     printmore = true;

         if (printmore)
         {
           std::cout << " elveci = " << elveci << std::endl;
           std::cout << " elmat = " << elmat << std::endl;
           std::cout << " elvec = " << elvec << std::endl;
         }

         // fes->TransformVec (i, bound, elveci, TRANSFORM_SOL);

         u->GetElementVector (dnums, elvec);
         elveci += elvec;
         u->SetElementVector (dnums, elveci);
	  
         for (int j = 0; j < dnums.Size(); j++)
           cnti[dnums[j]]++;
       });

    progress.Done();
    
#ifdef PARALLEL
    AllReduceDofData (cnti, MPI_SUM, fes->GetParallelDofs());
    u->GetVector().SetParallelStatus(DISTRIBUTED);
    u->GetVector().Cumulate(); 	 
#endif

    FlatVector<SCAL> fluxi(dim, clh);
    Array<int> dnums(1);
    for (int i = 0; i < cnti.Size(); i++)
      if (cnti[i])
      {
        dnums[0] = i;
        u->GetElementVector (dnums, fluxi);
        fluxi /= double (cnti[i]);
        u->SetElementVector (dnums, fluxi);
      }
    
    ma->PopStatus ();
  }
  
  void SetValuesX (const Array<shared_ptr<CoefficientFunction>> & coefs,
                   const TimeInterval & ti,
                   shared_ptr<GridFunction> u,
                   bool bound,
                   LocalHeap & clh)
  {
    if (u->GetFESpace()->IsComplex())
    {
      throw Exception("no complex yet");
      // if (u->GetMeshAccess().GetDimension() == 2)
      //   SetValuesX<2,Complex> (coefs, u, bound, clh);
      // else
      //   SetValuesX<3,Complex> (coefs, u, bound, clh);
    }
    else
      if (u->GetMeshAccess()->GetDimension() == 2)
        SetValuesX<2,double> (coefs, ti, u, bound, clh);
      else
        SetValuesX<3,double> (coefs, ti, u, bound, clh);
  }

///
  class NumProcSetValuesX : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
    bool boundary;
    bool coarsegridonly;
    int component;
    bool print;
    double told;
    double tnew;
  public:
    ///
    NumProcSetValuesX (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      coef_neg = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient_neg", ""));
      coef_pos = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient_pos", ""));
      told = flags.GetNumFlag ("told",0.0);
      tnew = flags.GetNumFlag ("new",1.0);
      boundary = flags.GetDefineFlag ("boundary");
      coarsegridonly = flags.GetDefineFlag ("coarsegridonly");
      component = int (flags.GetNumFlag ("component", 0))-1;
      print = flags.GetDefineFlag ("print");

      if (flags.NumFlagDefined ("component"))
      {
        cerr << "!!!!     numproc setvaluesX   ... -component   is depreciated and will be removed soon" << endl
             << "!!!!     please use  -gridfuncion=" << gfu->GetName() << "." << component << " instead" << endl;
      }
    }

    ///
    virtual ~NumProcSetValuesX() { ; }


    static void PrintDoc (ostream & ost)
    {
      ost << 
        "\n\nNumproc setvaluesX:\n"		\
        "-----------------\n"				\
        "Set a gridfunction to given values\n\n" \
        "Required flags:\n"
        "" \
        "-gridfunction=<gfname>\n"						\
        "    grid-function to be set\n"	\
        "-coefficient_neg=<coefname>\n"						\
        "    coefficient providing valuesX\n\n" \
        "-coefficient_pos=<coefname>\n"						\
        "    coefficient providing valuesX\n\n" \
        "-told=val\n"						\
        "    told for evaluating bnd function\n\n" \
        "-new=val\n"						\
        "    tnew for evaluating bnd function\n\n" \
        "\nOptional flags:\n"						\
        "-boundary\n only boundary valuesX are set\n" \
        "-component=<comp>\n set only this component (of CompoundFESpace)\n" \
        "-coarsegridonly\n" \
        "    set valuesX only on coarsest grid \n" 
          << endl;
    }
    
    ///
    virtual void Do(LocalHeap & lh)
    {
      if (coarsegridonly && ma->GetNLevels() > 1) return;
      shared_ptr<GridFunction> hgfu = gfu;
      if (component != -1)
        hgfu = gfu->GetComponent(component);

      Array<shared_ptr<CoefficientFunction> > coefs(2);
      coefs[0] = coef_neg;
      coefs[1] = coef_pos;
      TimeInterval ti(told,tnew);
      SetValuesX (coefs, ti, hgfu, boundary, lh);
      if (print) 
        *testout << "setvaluesX result:" << endl << hgfu->GetVector() << endl;
    }

    ///
    virtual string GetClassName () const
    {
      return "SetValuesX";
    }

    virtual void PrintReport (ostream & ost)
    {
      ost << GetClassName() << endl
          << "Gridfunction-Out = " << gfu->GetName() << endl;
    }
  };


 

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcXDifference : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    shared_ptr<GridFunction> gfu2;
    SolutionCoefficients<D> solcoef;
    double threshold;
    int intorder;
    double a_pos;
    double a_neg;
    double b_pos;
    double b_neg;
    double time;
    ErrorTable errtab;
    const Flags myflags;
    bool nooutput = false;
  public:


    NumProcXDifference (shared_ptr<PDE> apde, const Flags & flags)
        : NumProc (apde), solcoef(apde,flags), errtab(), myflags(flags)
    { 
      gfu  = apde->GetGridFunction (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));
      gfu2 = apde->GetGridFunction (flags.GetStringFlag ("solution2", flags.GetStringFlag("reference","")),true);
      threshold = flags.GetNumFlag ( "threshold", -0.1);
      intorder = (int) flags.GetNumFlag ( "intorder", 2);
      a_pos = flags.GetNumFlag ( "diffusion_p", 1.0);
      a_neg = flags.GetNumFlag ( "diffusion_n", 1.0);
      b_pos = flags.GetNumFlag ( "henryweight_p", 1.0);
      b_neg = flags.GetNumFlag ( "henryweight_n", 1.0);
      time = flags.GetNumFlag ( "time", 1.0);
      nooutput = flags.GetDefineFlag ("nooutput");
    }

    virtual ~NumProcXDifference()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcXDifference";
    }


    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;
      refinements++;
      CalcXError<D>(gfu, gfu2, solcoef, intorder, a_neg, a_pos, b_neg, b_pos, time, errtab, lh, !nooutput, myflags);
    }    
    

  };

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcSpecialOutput : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    SolutionCoefficients<D> * solcoef;
    int subdivision;
    bool onlygrid;
    Flags myflags;
  public:


    NumProcSpecialOutput (shared_ptr<PDE> apde, const Flags & flags)
        : NumProc (apde), solcoef(NULL), myflags(flags)
    { 
      gfu  = apde->GetGridFunction (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")),true);
      subdivision = (int) flags.GetNumFlag ( "subdivision", 2);
      onlygrid = flags.GetDefineFlag ("onlymesh");
      if (!onlygrid)
      {
          solcoef = new SolutionCoefficients<D>(apde,flags);
      }
    }

    virtual ~NumProcSpecialOutput()
    {
      if (solcoef)
        delete solcoef;
    }

    virtual string GetClassName () const
    {
      return "NumProcSpecialOutput";
    }


    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;
      refinements++;
      if (onlygrid)
        OutputMeshOnly(ma, lh);
      else
      {
        if (gfu==NULL)
          throw Exception("oh no oh nooo - give me a gridfunction next time! pleeeeaaase!");
        DoSpecialOutput<D>(gfu, *solcoef,subdivision, myflags, lh);
      }
    }    
    

  };


/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  class NumProcCalcCondition : public NumProc
  {
  protected:
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BilinearForm> bfid;
    string inversetype;
    bool symmetric;
    bool printmatrix;
    shared_ptr<GridFunction> gfu;
    shared_ptr<Preconditioner> pre;
  public:


    NumProcCalcCondition (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform","bfa"));
      bfid = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform_id","id"));
      inversetype = flags.GetStringFlag ("inverse", "pardiso");
      symmetric = flags.GetDefineFlag ("symmetric");
      printmatrix = flags.GetDefineFlag ("printmatrix");
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction","gfu"),true);
      pre = apde->GetPreconditioner (flags.GetStringFlag ("preconditioner","c"),true);
    }

    virtual ~NumProcCalcCondition()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcCalcCond";
    }


    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      cout << " This is the Do-call for CalcCondition on refinement level " << refinements << std::endl;
      refinements++;

      BaseMatrix & mata = bfa->GetMatrix();
      BaseMatrix & matid = bfid->GetMatrix();

	  dynamic_cast<BaseSparseMatrix&> (mata) . SetInverseType (inversetype);

      if (printmatrix)
      {
          ofstream outf("matrix.out");
          mata.Print(outf);

          ofstream outf2("precond.out");
          pre->GetMatrix().Print(outf2);
      }
      // shared_ptr<BaseMatrix> invmat = dynamic_cast<BaseSparseMatrix&> (mata) . InverseMatrix(bfa->GetFESpace()->GetFreeDofs());

      std::ofstream outf("condition_npcc.out");
      std::ofstream outjf("condition_jac_npcc.out");
      outf << refinements-1 << "\t";
      outjf << refinements-1 << "\t";
      CalcCond(mata, matid, bfa->GetFESpace()->GetFreeDofs(), true, false, &outf , symmetric);
      CalcCond(mata, matid, bfa->GetFESpace()->GetFreeDofs(), true, true, &outjf, symmetric, gfu);
      outjf << endl;

      // delete &invmat;
               
    }    
    

  };



  /* ---------------------------------------- 
                  numproc markinterface
   ---------------------------------------- */
  class NumProcMarkElementsOnInterface : public NumProc
  {
  protected:
    shared_ptr<FESpace> fes;
  public:
    
    NumProcMarkElementsOnInterface (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      fes = apde->GetFESpace(flags.GetStringFlag("fespace","xfes"));
      if (dynamic_pointer_cast<XFESpace>(fes) == nullptr)
      {
        auto compfes = dynamic_pointer_cast<CompoundFESpace>(fes);
        for (int i = 0; i < compfes->GetNSpaces(); ++i)
        {
          fes = dynamic_pointer_cast<XFESpace>((*compfes)[i]);
          if (fes != nullptr) break;
        }
      }
    }
  
    virtual string GetClassName () const
    {
      return "NPMarkElementsonInterface";
    }


    virtual void Do (LocalHeap & lh)
    {
      for (int i = 0; i < fes->GetMeshAccess()->GetNE(); ++i)
      {
        if (dynamic_pointer_cast<XFESpace>(fes)->IsElementCut(i))
          Ng_SetRefinementFlag (i+1, 1);
        else
          Ng_SetRefinementFlag (i+1, 0);
      }

      if (fes->GetMeshAccess()->GetDimension() == 3)
      {
        int nse = ma->GetNSE();
        for (int i = 0; i < nse; i++)
          Ng_SetSurfaceRefinementFlag (i+1, 0);
      }

    }
  };

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcTraceDiff : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    shared_ptr<CoefficientFunction> coef;
    int intorder;
    const Flags myflags;
    bool nooutput = false;
  public:


    NumProcTraceDiff (shared_ptr<PDE> apde, const Flags & flags)
        : NumProc (apde), myflags(flags)
    { 
      gfu  = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", flags.GetStringFlag("gridfunction","")));
      coef = apde->GetCoefficientFunction (flags.GetStringFlag ("coef", flags.GetStringFlag("coef","")),true);
      intorder = (int) flags.GetNumFlag ( "intorder", 2);
      nooutput = flags.GetDefineFlag ("nooutput");
    }

    virtual ~NumProcTraceDiff()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "NumProcTraceDiff";
    }


    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;
      refinements++;
      Array<double> errors;
      CalcTraceDiff<D>(gfu, coef, intorder, errors, lh);
      cout << " l2diff = " << errors[0] << endl;
      cout << " maxdiff = " << errors[1] << endl;

    }    
    

  };


}

static RegisterNumProc<NumProcSetValuesX> npinittestxfem2d("setvaluesx");
static RegisterNumProc<NumProcXDifference<2> > npxdiff("xdifference");
static RegisterNumProc<NumProcXDifference<3> > npxdiff3d("xdifference3d");
static RegisterNumProc<NumProcTraceDiff<2> > nptdiff("tracediff");
static RegisterNumProc<NumProcTraceDiff<3> > nptdiff3d("tracediff3d");
static RegisterNumProc<NumProcSpecialOutput<2> > npxoutp("xoutput");
static RegisterNumProc<NumProcCalcCondition> npinitcalccond("calccond");
static RegisterNumProc<NumProcMarkElementsOnInterface> npinitmark("markinterface");
