/*
*/


#include <solve.hpp>


using namespace ngfem;
using namespace ngsolve;


const double epsreg2 = 1e-8;
const double epsnumdiff =1e-8;

template <int D>
class LevelsetCorrectionIntegrator : public LinearFormIntegrator
{
public:
  shared_ptr<CoefficientFunction> coef_phi;
  Array<shared_ptr<CoefficientFunction>> coef_vel;
  
  LevelsetCorrectionIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
  {
    coef_phi = coeffs[0];
    coef_vel.SetSize(D);
    for (int d = 0; d < D; ++d)
      coef_vel[d] = coeffs[d+1];
  }
  
  virtual string Name () const { return "LevelsetCorrection"; }

  virtual bool BoundaryForm () const { return false; }
  
  // compute the Hesse Vector at point elveclin
  virtual void
  CalcElementVector (const FiniteElement & cfel,
                     const ElementTransformation & eltrans,
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
      const FiniteElement & bfel = 
        dynamic_cast<const CompoundFiniteElement&> (cfel) [0];


    auto & fel = static_cast<const ScalarFiniteElement<D>&> (bfel);
    int ndof = fel.GetNDof();

    elvec = 0;
    

    FlatMatrixFixWidth<D> dshape(ndof, lh);
    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        MappedIntegrationPoint<D,D> mip(ir[i], eltrans); 

        fel.CalcShape (ir[i], shape);

        Mat<D> trafo = Trans (mip.GetJacobianInverse());

        const double phi = coef_phi->Evaluate(mip);

        Vec<D> gradphi;
        Vec<D> gradphi_ref;
        for (int d = 0; d < D; ++d)
        {
          IntegrationPoint ipr (ir[i]); ipr(d) += epsnumdiff;
          IntegrationPoint ipl (ir[i]); ipl(d) -= epsnumdiff;
          MappedIntegrationPoint<D,D> mipr( ipr, eltrans);
          MappedIntegrationPoint<D,D> mipl( ipl, eltrans);
          gradphi_ref[d] = 0.5/epsnumdiff * ( coef_phi->Evaluate(mipr) - coef_phi->Evaluate(mipl));
        }

        gradphi = trafo * gradphi_ref;

        Mat<D> gradwmat;
        for (int k = 0; k < D; ++k)
          {
            Vec<D> gradw_ref;
            for (int d = 0; d < D; ++d)
              {
                IntegrationPoint ipr (ir[i]); ipr(d) += epsnumdiff;
                IntegrationPoint ipl (ir[i]); ipl(d) -= epsnumdiff;
                MappedIntegrationPoint<D,D> mipr( ipr, eltrans);
                MappedIntegrationPoint<D,D> mipl( ipl, eltrans);
                gradw_ref[d] = 0.5/epsnumdiff * ( coef_vel[k]->Evaluate(mipr) - coef_vel[k]->Evaluate(mipl));
              }
            Vec<D> gradw = trafo * gradw_ref;

            for (int d = 0; d < D; ++d)
              {
                gradwmat(k,d) = gradw(d);
              }
          }
        const double grad_wmat_grad = InnerProduct(gradphi, gradwmat * gradphi);
        const double grad_wmat_grad_normalized = grad_wmat_grad / (L2Norm2(gradphi) + epsreg2);
        
        // double gradnorm = 1.0;
        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();

        elvec += fac * grad_wmat_grad_normalized * phi * shape;
      }
  }
};

static RegisterLinearFormIntegrator<LevelsetCorrectionIntegrator<2> > initmcm2d ("lsetcorr", 2, 3);
static RegisterLinearFormIntegrator<LevelsetCorrectionIntegrator<3> > initmcm3d ("lsetcorr", 3, 4);
