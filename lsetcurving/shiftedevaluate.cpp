#define FILE_SHIFTEDEVALUATE_CPP
#include "shiftedevaluate.hpp"
#include <diffop_impl.hpp>
#include "../utils/ngsxstd.hpp"

namespace ngfem
{

  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {
    const MappedIntegrationPoint<SpaceD,SpaceD> & mip =
      static_cast<const MappedIntegrationPoint<SpaceD,SpaceD>&> (bmip);

    //const ScalarFiniteElement<SpaceD> & scafe =
            //dynamic_cast<const ScalarFiniteElement<SpaceD> & > (bfel);
    // const int ndof = bfel.GetNDof();

    //FlatVector<> shape (ndof,lh);

    IntegrationPoint ip(mip.IP());
    auto elid = mip.GetTransformation().GetElementId();
    Array<int> dnums;

    Vec<SpaceD> z = mip.GetPoint();

    if (forth)
    {
      forth->GetFESpace()->GetDofNrs(elid,dnums);
      FlatVector<> values_forth(dnums.Size()*SpaceD,lh);
      FlatMatrixFixWidth<SpaceD> vector_forth(dnums.Size(),&(values_forth(0)));
      forth->GetVector().GetIndirect(dnums,values_forth);

      FiniteElement& fe_forth = forth->GetFESpace()->GetFE(elid,lh);
      const ScalarFiniteElement<SpaceD> & scafe_forth =
        dynamic_cast<const ScalarFiniteElement<SpaceD> & > (fe_forth);
      FlatVector<> shape_forth(dnums.Size(),lh);
      scafe_forth.CalcShape(ip,shape_forth);
      Vec<SpaceD> dvec_forth = Trans(vector_forth)*shape_forth;
      z += dvec_forth;
    }
      

    // Solve the problem Theta(Phi(x)) = z

    if (back)
    {
      
      back->GetFESpace()->GetDofNrs(elid,dnums);
      FlatVector<> values_back(dnums.Size()*SpaceD,lh);
      FlatMatrixFixWidth<SpaceD> vector_back(dnums.Size(),&(values_back(0)));
      back->GetVector().GetIndirect(dnums,values_back);
    
      FiniteElement& fe_back = back->GetFESpace()->GetFE(elid,lh);
      const ScalarFiniteElement<SpaceD> & scafe_back =
        dynamic_cast<const ScalarFiniteElement<SpaceD> & > (fe_back);
      FlatVector<> shape_back(dnums.Size(),lh);
      scafe_back.CalcShape(ip,shape_back);
      Vec<SpaceD> dvec_back = Trans(vector_back)*shape_back;
    

      int its = 0;
      const double h = pow(abs(mip.GetJacobiDet()), 1./SpaceD);
      Vec<SpaceD> diff;
      IntegrationPoint ipx(ip);

      IntegrationPoint ipx0(0,0,0);
      MappedIntegrationPoint<SpaceD,SpaceD> mip_x0(ipx0,mip.GetTransformation());
      Vec<SpaceD> zdiff = z-mip_x0.GetPoint();
    
      // static atomic<int> cnt_its(0);
      // static atomic<int> cnt_calls(0);
      
      Vec<SpaceD> ipx_best_so_far;
      double diff_best_so_far;
      bool first = true; int idx_best;

      // Fixed point iteration
      while (its < globxvar.FIXED_POINT_ITER_TRESHOLD)
      {
        scafe_back.CalcShape(ipx,shape_back);
        dvec_back = Trans(vector_back)*shape_back;

        FlatVector<double> fv(SpaceD,&(ipx.Point())(0));

        diff = zdiff - dvec_back - mip.GetJacobian() * fv;
        //cout << "diff = " << diff << endl;
        //cout << "its = " << its << endl;
        if(first) {
            diff_best_so_far = L2Norm(diff);
            ipx_best_so_far = ipx.Point();
            idx_best = its;
            first = false;
        }
        else {
            if (L2Norm(diff) < diff_best_so_far){
                diff_best_so_far = L2Norm(diff);
                ipx_best_so_far = ipx.Point();
                idx_best = its;
            }
        }
        if ( L2Norm(diff) < globxvar.EPS_SHIFTED_EVAL*h ) break;
        ipx.Point() = mip.GetJacobianInverse() * (zdiff - dvec_back);

        its++;
        // cnt_its++;
      
      }
      if (its == globxvar.FIXED_POINT_ITER_TRESHOLD){
          if(diff_best_so_far < 1e0) {
              cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "In Shifted_eval: Not converged, but the "+to_string(idx_best)+"th iteration seems a reasonable candidate" << endl;
              ipx.Point() = ipx_best_so_far;
          }
          else {
              cout << "Last diff: " << diff << endl;
              cout << "Best diff: " << diff_best_so_far << endl;
              throw Exception(" shifted eval took FIXED_POINT_ITER_TRESHOLD = "+to_string(globxvar.FIXED_POINT_ITER_TRESHOLD)+" iterations and didn't (yet?) converge! In addition, the best interation step is no good fallback candidate.");
          }
      }
    
      // cnt_calls++;
      // cout << "cnt/calls = " << cnt_its/cnt_calls << endl;

      /* 
         FlatMatrixFixWidth<2> dshape_back(dnums.Size(),lh);
         FlatVector<> shape_back_new(dnums.Size(),lh);

         // not so robust Newton-type version (not working so well):
    
         while (its==0 || (L2Norm(diff) > 1e-8*h && its < 20000))
         {
         MappedIntegrationPoint<2,2> mip_x0(ipx,mip.GetTransformation());
         scafe_back.CalcShape(ipx,shape_back);
         dvec_back = Trans(vector_back)*shape_back;
         diff = z - dvec_back - mip_x0.GetPoint();
         Mat<2> dphi = mip_x0.GetJacobian();
         scafe_back.CalcDShape(ipx,dshape_back);
         Mat<2> dback_ref = Trans(dshape_back)*vector_back;
         Mat<2> dback = Trans(mip_x0.GetJacobianInverse())*dback_ref + dphi;
         Mat<2> invd;
         CalcInverse(dback,invd);

         // Backtracking-Armijo line search

         // scalars
         double alpha = 1;
         double beta = 0.1;
         int itsb = 0;
         // descent direction
         Vec<2> ddirection = invd * diff;
         // Jacobi matrix
         Mat<2> jacobi_old = dback;
         // function value at current Newton iterate
         Vec<2> value_old = -diff;
         // Stuff to calculate function value at
         // prospective new iterate
         IntegrationPoint ipx_new(x);
         // the difference in the Armijo condition of which all components
         // need to be <= 0
         Vec<2> diffb;

         while(itsb == 0 || ( ( (diffb(0) > 1e-8*h) || (diffb(1) > 1e-8*h) ) && itsb < 100) )
         {
         if(itsb != 0)
         alpha *= 0.98;
         ipx_new.Point() = x + alpha*ddirection;
         MappedIntegrationPoint<2,2> mip_xnew(ipx_new,mip.GetTransformation());
         scafe_back.CalcShape(ipx_new,shape_back_new);
         Vec<2> dvec_back_new = Trans(vector_back)*shape_back_new;
         Vec<2> value_new = -z + mip_xnew.GetPoint() + dvec_back_new;
         diffb = value_new - value_old - alpha*beta*jacobi_old*ddirection;
         itsb ++;
         // cout << "itsb = " << itsb << endl;
         // cout << "diffb(0) = " << diffb(0) << endl;
         // cout << "diffb(1) = " << diffb(1) << endl;
         }

         cout << "alpha = " << alpha << endl;


         its++;
         //alpha = 0.25;
         x += alpha*invd * diff;

         // cout << "x = " << x << endl;
         // cout << "z = " << z << endl;
         cout << "diff = " << diff << endl;
         cout << "its = " << its << endl;
         // cout << "dvec_back = " << dvec_back << endl;
         // cout << "dvec_forth = " << dvec_forth << endl;
         // cout << "dphi = " << dphi << endl;
         // cout << "dback_ref = " << dback_ref << endl;
         // cout << "dback = " << dback << endl;
         //getchar();

         cnt++;
         cout << "cnt = " << cnt << endl;
         ipx.Point() = x;

         }
      */
      MappedIntegrationPoint<SpaceD, SpaceD> mipx(ipx, mip.GetTransformation());
      evaluator->CalcMatrix(bfel, mipx, mat, lh);
      //scafe.CalcShape(ipx,shape);
      //mat = 0.0;
      //mat.Row(0) = shape;
      // for (int j = 0; j < D; j++)
      //   for (int k = 0; k < shape.Size(); k++)
      //     mat(j,k*D+j) = shape(k);      
    }
    else
    {
      int its = 0;
      const double h = pow(abs(mip.GetJacobiDet()), 1./SpaceD);
      Vec<SpaceD> diff;
      IntegrationPoint ipx(ip);
      IntegrationPoint ipx0(0,0,0);
      MappedIntegrationPoint<SpaceD,SpaceD> mip_x0(ipx0,mip.GetTransformation());
      Vec<SpaceD> zdiff = z-mip_x0.GetPoint();
    
      // Fixed point iteration
      while (its < globxvar.FIXED_POINT_ITER_TRESHOLD)
      {
        FlatVector<double> fv(SpaceD,&(ipx.Point())(0));
        diff = zdiff - mip.GetJacobian() * fv;
        if ( L2Norm(diff) < globxvar.EPS_SHIFTED_EVAL*h ) break;
        ipx.Point() = mip.GetJacobianInverse() * zdiff;
        its++;
      }
      if (its == globxvar.FIXED_POINT_ITER_TRESHOLD)
        throw Exception(" shifted eval took FIXED_POINT_ITER_TRESHOLD iterations and didn't (yet?) converge! ");

      MappedIntegrationPoint<SpaceD, SpaceD> mipx(ipx, mip.GetTransformation());
      evaluator->CalcMatrix(bfel, mipx, mat, lh);

      //scafe.CalcShape(ipx,shape);
      //mat = 0.0;
      // for (int j = 0; j < D; j++)
      //   for (int k = 0; k < shape.Size(); k++)
      //     mat(j,k*D+j) = shape(k);
      //mat.Row(0) = shape;
      //mat.Row(1) = shape;

      //cout << "mat:\n" << mat << endl;
    }
    
  }


  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof()*BlockDim(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              BareSliceVector<double> x,
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), fel.GetNDof()*BlockDim(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x.Range(0,fel.GetNDof()) = Trans(mat) * flux;
  }

  template class DiffOpShiftedEval<1>;
  template class DiffOpShiftedEval<2>;
  template class DiffOpShiftedEval<3>;

}


