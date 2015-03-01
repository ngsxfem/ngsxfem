#include "xstokesspace.hpp"
#include "xfemVisInts.hpp"

namespace ngcomp
{


/*
  Evaluate (u_x, u_y) 
*/
  class DiffOpIdU : public DiffOp<DiffOpIdU>
  {
    // 2 components:
    // u1 u2

  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2 };
    enum { DIFFORDER = 0 };
  
    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      const CompoundFiniteElement & cfela = 
        dynamic_cast<const CompoundFiniteElement&> (bfel);
      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (cfela[0]);
      const ScalarFiniteElement<2> & fel_u = 
        dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
    
      int nd_u = fel_u.GetNDof();

      FlatVector<> vecu(nd_u, lh);
      fel_u.CalcShape (mip.IP(), vecu);

      mat = 0;
      mat.Row(0).Range(cfela.GetRange(0)) = vecu;
      mat.Row(1).Range(cfela.GetRange(1)) = vecu;
    }
  };


  class StokesUIntegrator 
    : public T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
  {
  public:
    ///
    StokesUIntegrator (const Array<shared_ptr<CoefficientFunction>> & /* coeffs */)
      :  T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
         (DiagDMat<2> (make_shared<ConstantCoefficientFunction>(1)))
    { ; }

    ///
    virtual string Name () const { return "Stokes IdU"; }
  };




  XStokesFESpace::XStokesFESpace   (shared_ptr<MeshAccess> ama, 		   
                                    const Flags & flags)
    : CompoundFESpace(ama, flags)
  {
    name="XStokesFESpace";

    Flags velflags;
    Flags pflags;

    // standard FE-type:

    velflags.SetFlag("type", "xstdfespace");
    pflags.SetFlag("type", "xstdfespace");

    velflags.SetFlag("type_std", "h1ho");
    pflags.SetFlag("type_std", "h1ho");

    // "empty" - flags

    bool empty_vel = flags.GetDefineFlag ("empty_vel"); 
    bool empty_p = flags.GetDefineFlag ("empty_p"); 
    bool empty = flags.GetDefineFlag ("empty"); // overwrites the others
    if (empty)
    {
      empty_vel = true; 
      empty_p = true;
    }

    cout << " XStokesFESpace with";
    if (empty_vel)
      cout << "out";
    cout << " Heaviside enrichment for the velocity" << endl;

    if (empty_vel)
      velflags.SetFlag ("empty");

    cout << " XStokesFESpace with";
    if (empty_p)
      cout << "out";
    cout << " Heaviside enrichment for the pressure" << endl;

    if (empty_p)
      pflags.SetFlag ("empty");

    // "dgjumps" - flags

    bool dgjumps = flags.GetDefineFlag ("dgjumps"); 

    cout << " XStokesFESpace with";
    if (!dgjumps)
      cout << "out";
    cout << " dgjumps" << endl;

    if (dgjumps)
    {
      velflags.SetFlag ("dgjumps");
      pflags.SetFlag ("dgjumps");
    }

    // refinement - flags

    int ref = flags.GetNumFlag ("ref",0); 
    int ref_space = max(ref,(int)flags.GetNumFlag ("ref_space",0)); 
    int ref_time = max(ref,(int)flags.GetNumFlag ("ref_time",0)); 

    velflags.SetFlag ("ref_space", ref_space);
    velflags.SetFlag ("ref_time", ref_time);
    pflags.SetFlag ("ref_space", ref_space);
    pflags.SetFlag ("ref_time", ref_time);

    std::cout << " XStokesFESpace with ref_space = " << ref_space << std::endl;

    // order - flags

    int order = int (flags.GetNumFlag ("order", 1));
    if (order < 1) order = 1;

    cout << " TaylorHood with order "<< order << " (pressure order = " << order 
         << ", velocity order = " << order + 1<< ") " << endl;

    pflags.SetFlag ("order", order);
    velflags.SetFlag ("order", order+1);

    // boundary conditions

    if (flags.NumListFlagDefined ("dirichlet_vel"))
      velflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet_vel"));

    if (flags.NumListFlagDefined ("dirichlet_p"))
      pflags.SetFlag ("dirichlet", flags.GetNumListFlag ("dirichlet_p"));

      
    if (flags.NumListFlagDefined("definedon")){ 
      Array<double> defon;
      defon = flags.GetNumListFlag("definedon");
      velflags.SetFlag("definedon", defon);
      pflags.SetFlag("definedon", defon);
    }
      
    /// SPACES:

    const int D = ma->GetDimension();

    spaces.SetSize(0);
    for (int d = 0; d < D; ++d)
      AddSpace (XStdFESpace::Create(ma, velflags));

    AddSpace (XStdFESpace::Create(ma, pflags));

    shared_ptr<CoefficientFunction> one = make_shared<ConstantCoefficientFunction>(1);
    Array<shared_ptr<CoefficientFunction>> arr(1); arr[0] = one;
    if (ma->GetDimension() == 2)
    {
      integrator = make_shared< StokesUIntegrator> (arr);
      //   integrator = make_shared<MassVecHDGIntegrator<2>> (one);
      //   // boundary_integrator = new RobinIntegrator<2> (&one);
    }
    else
    {
      shared_ptr<BilinearFormIntegrator> integrator_inner 
        = make_shared<XVisIntegrator<3>>(one) ;
      integrator = make_shared< CompoundBilinearFormIntegrator> (integrator_inner, 2);      //   integrator = make_shared<MassVecHDGIntegrator<3>> (one);
      //   evaluator = make_shared<T_DifferentialOperator<DiffOpVecIdHDG<3>>>();
      //   boundary_integrator = make_shared<RobinVecHDGIntegrator<3>>(one);
    }
  }

  static RegisterFESpace<XStokesFESpace> initxstokes ("xstokes");


  NumProcInformXStokesFESpace::NumProcInformXStokesFESpace (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  { 
    shared_ptr<FESpace> xstokes = apde->GetFESpace(flags.GetStringFlag("xstokesfespace","v"));
    shared_ptr<CoefficientFunction> coef_lset = apde->GetCoefficientFunction(flags.GetStringFlag("coef_levelset","coef_lset"));
    dynamic_pointer_cast<XStokesFESpace>(xstokes) -> SetLevelSet (coef_lset);
  }
  
  NumProcInformXStokesFESpace::~NumProcInformXStokesFESpace(){ ; }
  string NumProcInformXStokesFESpace::GetClassName () const {return "InformXStokesFESpace";  }
  void NumProcInformXStokesFESpace::Do (LocalHeap & lh)  { ; }
  
  static RegisterNumProc<NumProcInformXStokesFESpace> npinfoxfe("informxstokes");


}
