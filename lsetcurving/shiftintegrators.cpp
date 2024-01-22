#include "shiftintegrators.hpp"

namespace ngfem
{

  template <int D>
  ShiftIntegrator<D> :: ShiftIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : coef_lset_p1(coeffs[0]),coef_lset_ho(coeffs[1]) 
  {
    
    if (coeffs.Size() > 2)
    {
      // std::cout << " ShiftIntegrator called with more than 2 arguments " << std::endl;
      max_deform = coeffs[2]->EvaluateConst();
    }
    if (coeffs.Size() > 3)
    {
      // std::cout << " ShiftIntegrator called with more than 3 arguments " << std::endl;
      lower_lset_bound = coeffs[3]->EvaluateConst();
    }
    if (coeffs.Size() > 4)
    {
      // std::cout << " ShiftIntegrator called with more than 4 arguments " << std::endl;
      upper_lset_bound = coeffs[4]->EvaluateConst();
    }
    if (coeffs.Size() > 5)
    {
      // std::cout << " ShiftIntegrator called with more than 5 arguments " << std::endl;
      qn = coeffs[5];
    }
    if (coeffs.Size() > 6)
    {
      // std::cout << " ShiftIntegrator called with more than 5 arguments " << std::endl;
      coef_blending = coeffs[6];
    }
  }
  
  template <int D>
  void ShiftIntegrator<D> :: CalcElementVector (const FiniteElement & fel,
                                                const ElementTransformation & eltrans,
                                                FlatVector<double> elvec,
                                                LocalHeap & lh,
                                                shared_ptr<LsetEvaluator<D>> lseteval) const
  {
    static Timer timer ("ShiftIntegrator<D>::CalcElementVector"); 
    RegionTimer reg (timer);
    
    elvec = 0.0;
    const ScalarFiniteElement<D> & scafe = dynamic_cast<const ScalarFiniteElement<D> &>(fel);

    if (!lseteval)
      lseteval = make_shared<LsetEvaluator<D>>(coef_lset_ho, eltrans);
    
    FlatMatrixFixWidth<D> elvecmat(scafe.GetNDof(),&elvec(0));
    elvecmat = 0.0;
    
    FlatVector<> shape (scafe.GetNDof(),lh);
      
    Vec<D> grad;
    if (!qn) //grad is constant on element...
    {
      IntegrationPoint ip(0.0,0.0,0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      CalcGradientOfCoeff(coef_lset_p1, mip, grad, lh);
    }
    
    IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*scafe.Order());
    for (int l = 0 ; l < ir.GetNIP(); l++)
    {
      MappedIntegrationPoint<D,D> mip(ir[l], eltrans);
      scafe.CalcShape(ir[l],shape);

      Mat<D> trafo_of_normals = mip.GetJacobianInverse() * Trans(mip.GetJacobianInverse());
        
      if (qn)
        qn->Evaluate(mip,grad);

      Vec<D> normal = mip.GetJacobianInverse() * grad;
      // double len = L2Norm(normal);
      // normal /= len;
        
      Vec<D> orig_point;
      for (int d = 0; d < D; ++d)
        orig_point(d) = ir[l](d);

      const double lsetp1val = coef_lset_p1->Evaluate(mip);
                                                                                 
      // const double h = D == 2 ? sqrt(2) * sqrt(mip.GetMeasure()) : sqrt(3) * cbrt(mip.GetMeasure());
      // double alpha = abs(lsetp1val / h) * abs(lsetp1val / h);

      double alpha = 0.0; // blending factor, 0.0 means: find phi_lin, 1.0 means: find phi (i.e. goal value = start value)
      if (coef_blending)
        alpha = coef_blending->Evaluate(mip);

      if (alpha > 1+1e-6)
        throw Exception("alpha should not be larger than 1");
      
      double goal_val = (1.0-alpha) * lsetp1val + alpha * lseteval->Evaluate(mip.IP(),lh);

      // double goal_val = coef_lset_p1->Evaluate(mip);
      
      Vec<D> final_point;
      SearchCorrespondingPoint<D>(*lseteval,
                                  orig_point, goal_val, 
                                  trafo_of_normals, normal, false,
                                  final_point, lh);
      Vec<D> ref_dist = (final_point - orig_point);
      const double ref_dist_size = L2Norm(ref_dist);
      if ((max_deform >= 0.0) && (ref_dist_size > max_deform))
      {
        ref_dist *= max_deform / ref_dist_size; 
      }

        
      Vec<D> deform = mip.GetJacobian() * ref_dist;
      // if (qn)
      //   deform = InnerProduct(deform,qnormal) * qnormal;

      elvecmat += mip.GetWeight() * shape * Trans(deform);
    }      
  }


  template class ShiftIntegrator<1>;
  template class ShiftIntegrator<2>;
  template class ShiftIntegrator<3>;

  
}

