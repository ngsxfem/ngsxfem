/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
// #include "xintegration.hpp"
#include "../spacetime/spacetimefespace.hpp"
#include "xfemIntegrators.hpp"
#include "stxfemIntegrators.hpp"
#include "setvaluesx.hpp"
#include "../utils/error.hpp"
#include "../utils/output.hpp"
#include "../utils/calccond.hpp"

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
namespace ngcomp
{ 
  template <int D, class SCAL>
  void SetValuesX (const Array<CoefficientFunction *> & acoefs,
                   const TimeInterval & ti,
                   GridFunction & bu,
                   bool bound,
                   LocalHeap & clh)
  {
    static Timer sv("timer setvaluesX"); RegionTimer r(sv);

    S_GridFunction<SCAL> & u = dynamic_cast<S_GridFunction<SCAL> &> (bu);

    const FESpace & fes = u.GetFESpace();
    const MeshAccess & ma = fes.GetMeshAccess();
    
    ma.PushStatus ("setvalues");

    VorB vorb = VorB(bound);

    if (vorb == VOL)
      throw Exception ("not yet implemented");

    int dim   = fes.GetDimension();

    // const BilinearFormIntegrator & bli = *fes.GetIntegrator(bound);

    // if (&bli == NULL)
    //   throw Exception ("no integrator available");

    // int dimflux = 1; //diffop ? diffop->Dim() : bli.DimFlux(); 
    
    Array<int> cnti(fes.GetNDof());
    cnti = 0;

    u.GetVector() = 0.0;

    BilinearFormIntegrator * bfi = NULL;
    LinearFormIntegrator * lfi = NULL;

    const CompoundFESpace & cfes = dynamic_cast<const CompoundFESpace & >(fes);
    const SpaceTimeFESpace * stfes = dynamic_cast<const SpaceTimeFESpace * >(cfes[0]);
    
    Array<CoefficientFunction *> coefs(acoefs);
    Array<CoefficientFunction *> coefs_one(2);
    ConstantCoefficientFunction one(1.0);
    coefs_one[0] = &one;
    coefs_one[1] = &one;

    ConstantCoefficientFunction told(ti.first);
    ConstantCoefficientFunction tnew(ti.second);

    if (stfes != NULL)
    {
      coefs.Append(&told);
      coefs.Append(&tnew);
      lfi = new (clh) SpaceTimeXNeumannIntegrator<D>(coefs);
      coefs_one.Append(&told);
      coefs_one.Append(&tnew);
      bfi = new (clh) SpaceTimeXRobinIntegrator<D>(coefs_one);
    }
    else
    {
      lfi = new (clh) XNeumannIntegrator<D>(coefs);
      bfi = new (clh) XRobinIntegrator<D>(coefs_one);
    }

    ProgressOutput progress (ma, "setvaluesx element", ma.GetNE());

    IterateElements 
      (fes, vorb, clh, 
       [&] (ElementId ei, LocalHeap & lh)
       {
         progress.Update ();

         // cout << "  ? on the boundary -> getchar(); " << endl; getchar();

         // if (bound && !fes.IsDirichletBoundary(ma.GetSElIndex(ei.Nr())))
         //   return;

         if (bound && !(cfes[0]->IsDirichletBoundary(ma.GetSElIndex(ei.Nr()))))
           return;

         const FiniteElement & bfel = fes.GetFE (ei, lh);

         const ElementTransformation & eltrans = ma.GetTrafo (ei, lh); 

         Array<int> dnums (bfel.GetNDof(), lh);
         fes.GetDofNrs (ei, dnums);

         FlatVector<SCAL> elvec(dnums.Size() * dim, lh);
         FlatVector<SCAL> elveci(dnums.Size() * dim, lh);
         FlatMatrix<double> elmat(dnums.Size(), lh);

         if (dim > 1)
           throw Exception ("dim > 1 not yet implemented");

         bfi->CalcElementMatrix (bfel, eltrans, elmat, lh);
         lfi->CalcElementVector (bfel, eltrans, elvec, lh);

         // for (int i = 0; i < elmat.Width(); ++i)
         //   elmat(i,i) += 1e-8;

         fes.TransformMat (ei.Nr(), bound, elmat, TRANSFORM_MAT_LEFT_RIGHT);
         fes.TransformVec (ei.Nr(), bound, elvec, TRANSFORM_RHS);

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

         // fes.TransformVec (i, bound, elveci, TRANSFORM_SOL);

         u.GetElementVector (dnums, elvec);
         elveci += elvec;
         u.SetElementVector (dnums, elveci);
	  
         for (int j = 0; j < dnums.Size(); j++)
           cnti[dnums[j]]++;
       });

    progress.Done();
    
#ifdef PARALLEL
    AllReduceDofData (cnti, MPI_SUM, fes.GetParallelDofs());
    u.GetVector().SetParallelStatus(DISTRIBUTED);
    u.GetVector().Cumulate(); 	 
#endif

    FlatVector<SCAL> fluxi(dim, clh);
    Array<int> dnums(1);
    for (int i = 0; i < cnti.Size(); i++)
      if (cnti[i])
      {
        dnums[0] = i;
        u.GetElementVector (dnums, fluxi);
        fluxi /= double (cnti[i]);
        u.SetElementVector (dnums, fluxi);
      }
    
    ma.PopStatus ();
  }
  
  void SetValuesX (const Array<CoefficientFunction *> & coefs,
                   const TimeInterval & ti,
                   GridFunction & u,
                   bool bound,
                   LocalHeap & clh)
  {
    if (u.GetFESpace().IsComplex())
    {
      throw Exception("no complex yet");
      // if (u.GetMeshAccess().GetDimension() == 2)
      //   SetValuesX<2,Complex> (coefs, u, bound, clh);
      // else
      //   SetValuesX<3,Complex> (coefs, u, bound, clh);
    }
    else
      if (u.GetMeshAccess().GetDimension() == 2)
        SetValuesX<2,double> (coefs, ti, u, bound, clh);
      else
        SetValuesX<3,double> (coefs, ti, u, bound, clh);
  }

///
  class NumProcSetValuesX : public NumProc
  {
  protected:
    GridFunction * gfu;
    CoefficientFunction * coef_neg;
    CoefficientFunction * coef_pos;
    bool boundary;
    bool coarsegridonly;
    int component;
    bool print;
    double told;
    double tnew;
  public:
    ///
    NumProcSetValuesX (PDE & apde, const Flags & flags)
      : NumProc (apde)
    {
      gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      coef_neg = pde.GetCoefficientFunction (flags.GetStringFlag ("coefficient_neg", ""));
      coef_pos = pde.GetCoefficientFunction (flags.GetStringFlag ("coefficient_pos", ""));
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
      if (coarsegridonly && ma.GetNLevels() > 1) return;
      GridFunction * hgfu = gfu;
      if (component != -1)
        hgfu = gfu->GetComponent(component);

      Array<CoefficientFunction *> coefs(2);
      coefs[0] = coef_neg;
      coefs[1] = coef_pos;
      TimeInterval ti(told,tnew);
      SetValuesX (coefs, ti, *hgfu, boundary, lh);
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
    GridFunction * gfu;
    SolutionCoefficients<D> solcoef;
    double threshold;
    int intorder;
    double b_pos;
    double b_neg;
    double time;
    ErrorTable errtab;
  public:


    NumProcXDifference (PDE & apde, const Flags & flags)
      : NumProc (apde), solcoef(apde,flags), errtab()
    { 
      gfu  = pde.GetGridFunction (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));
      threshold = flags.GetNumFlag ( "threshold", -0.1);
      intorder = (int) flags.GetNumFlag ( "intorder", 2);
      b_pos = flags.GetNumFlag ( "henryweight_p", 1.0);
      b_neg = flags.GetNumFlag ( "henryweight_n", 1.0);
      time = flags.GetNumFlag ( "time", 1.0);
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
      CalcXError<D>(gfu, solcoef, intorder, b_neg, b_pos, time, errtab, lh, true);
    }    
    

  };

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcSpecialOutput : public NumProc
  {
  protected:
    GridFunction * gfu;
    SolutionCoefficients<D> solcoef;
    int subdivision;
  public:


    NumProcSpecialOutput (PDE & apde, const Flags & flags)
      : NumProc (apde), solcoef(apde,flags)
    { 
      gfu  = pde.GetGridFunction (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));
      subdivision = (int) flags.GetNumFlag ( "subdivision", 2);
    }

    virtual ~NumProcSpecialOutput()
    {
      ;
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
      DoSpecialOutput<D>(gfu, solcoef,subdivision, lh);
    }    
    

  };


/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  class NumProcCalcCondition : public NumProc
  {
  protected:
    BilinearForm * bfa;
    string inversetype;
    bool symmetric;
    bool printmatrix;
    GridFunction * gfu;
    Preconditioner * pre;
  public:


    NumProcCalcCondition (PDE & apde, const Flags & flags)
      : NumProc (apde)
    { 
      bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform","bfa"));
      inversetype = flags.GetStringFlag ("inverse", "pardiso");
      symmetric = flags.GetDefineFlag ("symmetric");
      printmatrix = flags.GetDefineFlag ("printmatrix");
      gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction","gfu"),true);
      pre = pde.GetPreconditioner (flags.GetStringFlag ("preconditioner","c"),true);
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

	  dynamic_cast<BaseSparseMatrix&> (mata) . SetInverseType (inversetype);

      if (printmatrix)
      {
          ofstream outf("matrix.out");
          mata.Print(outf);

          ofstream outf2("precond.out");
          pre->GetMatrix().Print(outf2);
      }
	  BaseMatrix & invmat = *dynamic_cast<BaseSparseMatrix&> (mata) . InverseMatrix(bfa->GetFESpace().GetFreeDofs());

      std::ofstream outf("condition_npcc.out");
      std::ofstream outjf("condition_jac_npcc.out");
      outf << refinements-1 << "\t";
      outjf << refinements-1 << "\t";
      CalcCond(mata, invmat, bfa->GetFESpace().GetFreeDofs(), true, false, &outf , symmetric);
      CalcCond(mata, invmat, bfa->GetFESpace().GetFreeDofs(), true, true, &outjf, symmetric, gfu);
      outjf << endl;

      delete &invmat;
               
    }    
    

  };

}

static RegisterNumProc<NumProcSetValuesX> npinittestxfem2d("setvaluesx");
static RegisterNumProc<NumProcXDifference<2> > npxdiff("xdifference");
static RegisterNumProc<NumProcSpecialOutput<2> > npxoutp("xoutput");
static RegisterNumProc<NumProcCalcCondition> npinitcalccond("calccond");
