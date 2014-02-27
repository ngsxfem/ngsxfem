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
// #include "../utils/fieldeval.hpp"
#include "xfemIntegrators.hpp"
#include "stxfemIntegrators.hpp"
#include "setvaluesx.hpp"

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

    int dimflux = 1; //diffop ? diffop->Dim() : bli.DimFlux(); 
    
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

    ProgressOutput progress (ma, "setvalues element", ma.GetNE());

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
         FlatArray<int> dnums = fes.GetDofNrs (ei, lh);

         FlatVector<SCAL> elvec(dnums.Size() * dim, lh);
         FlatVector<SCAL> elveci(dnums.Size() * dim, lh);
         FlatMatrix<double> elmat(dnums.Size(), lh);

         if (dim > 1)
           throw Exception ("dim > 1 not yet implemented");

         bfi->CalcElementMatrix (bfel, eltrans, elmat, lh);
         lfi->CalcElementVector (bfel, eltrans, elvec, lh);
         
         fes.TransformMat (ei.Nr(), bound, elmat, TRANSFORM_MAT_LEFT_RIGHT);
         fes.TransformVec (ei.Nr(), bound, elvec, TRANSFORM_RHS);
              
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

}

static RegisterNumProc<NumProcSetValuesX> npinittestxfem2d("setvaluesx");
