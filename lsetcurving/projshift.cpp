#include "projshift.hpp"
#include "calcpointshift.hpp"
#include "shiftintegrators.hpp"

namespace ngcomp
{

  void ProjectShift (shared_ptr<GridFunction> lset_ho, shared_ptr<GridFunction> lset_p1,
                     shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn,
                     shared_ptr<BitArray> ba,
                     shared_ptr<CoefficientFunction> blending,
                     double lower_lset_bound, double upper_lset_bound, double threshold,
                     LocalHeap & clh)
  {
    static Timer timer("LsetCurv::ProjectShift"); RegionTimer reg (timer);
    
    auto ma = lset_p1->GetMeshAccess();
    int ne=ma->GetNE();
    int D =ma->GetDimension();

    shared_ptr<BilinearFormIntegrator> mass;
    if (D==1)
      mass = make_shared<MassIntegrator<1>>(make_shared<ConstantCoefficientFunction>(1.0));
    else if (D==2)
      mass = make_shared<MassIntegrator<2>>(make_shared<ConstantCoefficientFunction>(1.0));
    else
      mass = make_shared<MassIntegrator<3>>(make_shared<ConstantCoefficientFunction>(1.0));

    Array<shared_ptr<CoefficientFunction>> shift_array;

    shift_array.Append(lset_p1);
    shift_array.Append(lset_ho);
    shift_array.Append(make_shared<ConstantCoefficientFunction>(threshold));
    shift_array.Append(make_shared<ConstantCoefficientFunction>(lower_lset_bound));
    shift_array.Append(make_shared<ConstantCoefficientFunction>(upper_lset_bound));
    shift_array.Append(qn);
    shift_array.Append(blending);
    
    shared_ptr<ShiftIntegrator<1>> shift1D;
    shared_ptr<ShiftIntegrator<2>> shift2D;
    shared_ptr<ShiftIntegrator<3>> shift3D;
    if (D==1)
      shift1D = make_shared<ShiftIntegrator<1>>(shift_array);
    else if (D==2)
      shift2D = make_shared<ShiftIntegrator<2>>(shift_array);
    else
      shift3D = make_shared<ShiftIntegrator<3>>(shift_array);

    shared_ptr<BaseVector> factor = deform->GetVector().CreateVector();
    *factor = 0.0;
    deform->GetVector() = 0.0;

    ProgressOutput progress (ma, "project shift on element", ma->GetNE());

    IterateElements
      (*(deform->GetFESpace()), VOL, clh,  [&] (FESpace::Element el, LocalHeap & lh)
       {
         int elnr = el.Nr();
         HeapReset hr(lh);
         progress.Update();
      
         Array<int> p1_dofs;
         lset_p1->GetFESpace()->GetDofNrs(el,p1_dofs);
         FlatVector<> vals(p1_dofs.Size(),lh);
         lset_p1->GetVector().GetIndirect(p1_dofs,vals);
         const ElementTransformation & eltrans = el.GetTrafo();
      
         if ( (ba) && ( !(ba->Test(elnr)) ))
           return;
         if ( (!ba) && !ElementInRelevantBand(vals, lower_lset_bound, upper_lset_bound) )
           return;

         int ndofs = el.GetDofs().Size();
         const FiniteElement & fel_deform = el.GetFE();
         // FlatVector<> vals(def_dofs.Size(),lh);
         // lset_p1->GetVector().GetIndirect(def_dofs,vals);
         FlatMatrix<> massmat (ndofs,lh);
         FlatVector<> elvec (D*ndofs,lh);
         FlatVector<> elres (D*ndofs,lh);
         mass->CalcElementMatrix(fel_deform, eltrans, massmat, lh);
         CalcInverse(massmat);

      
         Array<int> lset_ho_dofs;
         lset_ho->GetFESpace()->GetDofNrs(el,lset_ho_dofs);
         FlatVector<> lset_ho_vals(lset_ho_dofs.Size(),lh);
         lset_ho->GetVector().GetIndirect(lset_ho_dofs,lset_ho_vals);
         const FiniteElement & fel_lset_ho = lset_ho->GetFESpace()->GetFE(el,lh);
      
         if (D==1)
         {
           const ScalarFiniteElement<1> & scafe_lset_ho = dynamic_cast< const ScalarFiniteElement<1> &>(fel_lset_ho);
           shared_ptr<LsetEvaluator<1>> lseteval = make_shared<LsetEvaluator<1>>(scafe_lset_ho,lset_ho_vals);
           shift1D->CalcElementVector(fel_deform, eltrans, elvec, lh, lseteval);
        
           FlatMatrixFixWidth<1> elvec_vec(ndofs,&elvec(0));
           FlatMatrixFixWidth<1> shift_vec(ndofs,&elres(0));
           for (int d = 0; d < D; ++d)
             shift_vec.Col(d) = massmat * elvec_vec.Col(d);
           // vertex values to zero
           for (int l = 0; l < D+1; ++l)
             shift_vec.Row(l) = 0.0;
         }
         else if (D==2)
         {
           const ScalarFiniteElement<2> & scafe_lset_ho = dynamic_cast< const ScalarFiniteElement<2> &>(fel_lset_ho);
           shared_ptr<LsetEvaluator<2>> lseteval = make_shared<LsetEvaluator<2>>(scafe_lset_ho,lset_ho_vals);
           shift2D->CalcElementVector(fel_deform, eltrans, elvec, lh, lseteval);
        
           FlatMatrixFixWidth<2> elvec_vec(ndofs,&elvec(0));
           FlatMatrixFixWidth<2> shift_vec(ndofs,&elres(0));
           for (int d = 0; d < D; ++d)
             shift_vec.Col(d) = massmat * elvec_vec.Col(d);
           // vertex values to zero
           for (int l = 0; l < D+1; ++l)
             shift_vec.Row(l) = 0.0;
         }
         else
         {
           const ScalarFiniteElement<3> & scafe_lset_ho = dynamic_cast< const ScalarFiniteElement<3> &>(fel_lset_ho);
           shared_ptr<LsetEvaluator<3>> lseteval = make_shared<LsetEvaluator<3>>(scafe_lset_ho,lset_ho_vals);
        
           shift3D->CalcElementVector(fel_deform, eltrans, elvec, lh, lseteval);
        
           FlatMatrixFixWidth<3> elvec_vec(ndofs,&elvec(0));
           FlatMatrixFixWidth<3> shift_vec(ndofs,&elres(0));
           for (int d = 0; d < D; ++d)
             shift_vec.Col(d) = massmat * elvec_vec.Col(d);
           // vertex values to zero
           for (int l = 0; l < D+1; ++l)
             shift_vec.Row(l) = 0.0;
         }

         FlatVector<> def_vals(D*ndofs,lh);
         deform->GetVector().GetIndirect(el.GetDofs(),def_vals);
         def_vals += elres;
         deform->GetVector().SetIndirect(el.GetDofs(),def_vals);

         // increase dof counter
         FlatVector<> factors(D*ndofs,lh);
         factor->GetIndirect(el.GetDofs(),factors);
         for (int k = 0; k < factors.Size(); k+=D)
           factors(k) += 1.0 ;
         factor->SetIndirect(el.GetDofs(),factors);

         
       });
    
    progress.Done();

    factor->SetParallelStatus(DISTRIBUTED);
    factor->Cumulate(); 	 
    deform->GetVector().SetParallelStatus(DISTRIBUTED);
    deform->GetVector().Cumulate(); 	 
    
    if (task_manager)
    {
      SharedLoop sl (Range (factor->Size()));
      task_manager->CreateJob
        ( [&] (const TaskInfo & ti) {
          LocalHeap lh = clh.Split();
          
          // averaging of the (summed) deformation
          Array<int> dnums(1);

          for (int i : sl)
          {
            HeapReset hr(lh);
            FlatVector<> val_fac(D,lh);
            FlatVector<> values(D,lh);
            dnums[0] = i;
            deform->GetVector().GetIndirect(dnums,values);
            factor->GetIndirect(dnums,val_fac);
            if (val_fac(0) > 0)
              values *= 1.0/val_fac(0);
            deform->GetVector().SetIndirect(dnums,values);
          }
        });
    }
    else
    {
      Array<int> dnums(1);
      for (int i : Range(factor->Size()) )
      {
        HeapReset hr(clh);
        FlatVector<> val_fac(D,clh);
        FlatVector<> values(D,clh);
        dnums[0] = i;
        deform->GetVector().GetIndirect(dnums,values);
        factor->GetIndirect(dnums,val_fac);
        if (val_fac(0) > 0)
              values *= 1.0/val_fac(0);
        deform->GetVector().SetIndirect(dnums,values);
      }
    }
  }
  
}

