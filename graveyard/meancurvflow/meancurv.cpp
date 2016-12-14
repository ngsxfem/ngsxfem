/*

Solve the nonlinear problem mean curvature flow problem

1/g(|grad u|) * du/ dt = div ( 1/g(|grad u|) grad u )

with 

g(x) = sqrt( eps^2 + x^2)

*/


#include <solve.hpp>


using namespace ngfem;
using namespace ngsolve;


const double epsreg2 = 1e-8;

template <int D>
class MeanCurvatureMassIntegrator : public BilinearFormIntegrator
{
public:
  
  MeanCurvatureMassIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
  { ; }
  
  virtual string Name () const { return "MeanCurvMass"; }

  virtual VorB VB () const { return VOL; }
  virtual bool IsSymmetric () const { return true; }
  
  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    FlatVector<> elveclin(fel.GetNDof(), lh);
    elveclin = 0;
    CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
  }


  // compute the Hesse Matrix at point elveclin
  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> elveclin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    const ScalarFiniteElement<D> & fel = static_cast<const ScalarFiniteElement<D>&> (bfel);
    int ndof = fel.GetNDof();

    elmat = 0;
    

    FlatMatrixFixWidth<D> dshape(ndof, lh);
    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        MappedIntegrationPoint<D,D> mip(ir[i], eltrans); 

        fel.CalcMappedDShape (mip, dshape);
        fel.CalcShape (ir[i], shape);

        Vec<D> graduilin = Trans(dshape) * elveclin;

        double gradnorm = sqrt( L2Norm2(graduilin) + epsreg2); 
        // double gradnorm = 1.0;
        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();

        elmat += (fac / gradnorm) * shape * Trans(shape);
      }
  }
};

template <int D>
class MeanCurvatureStiffnessIntegrator : public BilinearFormIntegrator
{
public:
  
  MeanCurvatureStiffnessIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
  { ; }
  
  virtual string Name () const { return "MeanCurvStiff"; }

  virtual VorB VB () const { return VOL; }
  virtual bool IsSymmetric () const { return true; }
  
  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    FlatVector<> elveclin(fel.GetNDof(), lh);
    elveclin = 0;
    CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
  }


  // compute the Hesse Matrix at point elveclin
  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> elveclin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    const ScalarFiniteElement<D> & fel = static_cast<const ScalarFiniteElement<D>&> (bfel);
    int ndof = fel.GetNDof();

    elmat = 0;
    

    FlatMatrixFixWidth<D> dshape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        MappedIntegrationPoint<D,D> mip(ir[i], eltrans); 

        fel.CalcMappedDShape (mip, dshape);

        Vec<D> graduilin = Trans(dshape) * elveclin;

        // double gradnorm = 1.0;
        double gradnorm = sqrt( L2Norm2(graduilin) + epsreg2); 

        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();

        elmat += (fac / gradnorm) * dshape * Trans(dshape);
      }
  }
};





class NumProcMeanCurvature : public NumProc
{
protected:
  shared_ptr<BilinearForm> bfa;
  shared_ptr<BilinearForm> bfm;
  shared_ptr<GridFunction> gfu;

  double tend;
  double dt;

  bool stopeachtimestep;
public:
  NumProcMeanCurvature (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  { 
    bfm = apde->GetBilinearForm (flags.GetStringFlag ("bilinearformm", "m"));
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

    dt = flags.GetNumFlag ("dt", 0.001);
    tend = flags.GetNumFlag ("tend", 1);
    stopeachtimestep = flags.GetDefineFlag("stop_after_step");

  }

  virtual ~NumProcMeanCurvature()
  { ; }


  virtual void Do(LocalHeap & lh)
  {
    cout << "lset mean curvature solver called" << endl;

    BaseVector & vecu = gfu->GetVector();
    // const BaseVector & vecf = lff->GetVector();

    shared_ptr<BaseVector> d = vecu.CreateVector();
    shared_ptr<BaseVector> w = vecu.CreateVector();

    BaseMatrix& matm = bfm->GetMatrix();
    BaseMatrix& mata = bfa->GetMatrix();

    cout << " before time loop " << endl;
    Ng_Redraw ();
    getchar();
    
    for (double t = 0; t <= tend; t += dt)
    {
      cout << "\r t = " << t;
      *d = 0.0;
      bfm -> AssembleLinearization (vecu, lh);
      bfa -> AssembleLinearization (vecu, lh);

      shared_ptr<BaseMatrix> summat = matm.CreateMatrix();
      summat->AsVector() = (1.0/dt) * matm.AsVector() + mata.AsVector();

      shared_ptr<BaseMatrix> invmstar = summat->InverseMatrix();

      *d = mata * vecu;
      *w = *invmstar * *d;

      vecu -= *w;

      // update visualization
      Ng_Redraw ();
      if (stopeachtimestep)
        getchar();
    }
    cout << endl;
  }
};



static RegisterNumProc<NumProcMeanCurvature> npinit("meancurv");
static RegisterBilinearFormIntegrator<MeanCurvatureMassIntegrator<2> > initmcm2d ("meancurv_mass", 2, 0);
static RegisterBilinearFormIntegrator<MeanCurvatureStiffnessIntegrator<2> > initmca2d ("meancurv_stiff", 2, 0);
static RegisterBilinearFormIntegrator<MeanCurvatureMassIntegrator<3> > initmcm3d ("meancurv_mass", 3, 0);
static RegisterBilinearFormIntegrator<MeanCurvatureStiffnessIntegrator<3> > initmca3d ("meancurv_stiff", 3, 0);
