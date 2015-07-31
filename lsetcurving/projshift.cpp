#include "projshift.hpp"
#include "calcpointshift.hpp"
#include "shiftintegrators.hpp"

namespace ngcomp
{

  void ProjectShift (shared_ptr<GridFunction> lset_ho, shared_ptr<GridFunction> lset_p1,
                     shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn,
                     double lower_lset_bound, double upper_lset_bound, double threshold,
                     LocalHeap & lh)
  {
    static Timer time_fct ("ProjectShift");
    RegionTimer reg (time_fct);
    
    auto ma = lset_p1->GetMeshAccess();
    int ne=ma->GetNE();
    int D =ma->GetDimension();

    shared_ptr<BilinearFormIntegrator> mass;
    if (D==2)
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
    
    shared_ptr<ShiftIntegrator<2>> shift2D;
    shared_ptr<ShiftIntegrator<3>> shift3D;
    if (D==2)
      shift2D = make_shared<ShiftIntegrator<2>>(shift_array);
    else
      shift3D = make_shared<ShiftIntegrator<3>>(shift_array);

    shared_ptr<BaseVector> factor = deform->GetVector().CreateVector();
    *factor = 0.0;
    deform->GetVector() = 0.0;

    ProgressOutput progress (ma, "project shift on element", ma->GetNE());
    
    for (int elnr = 0; elnr < ne; ++elnr)
    {
      HeapReset hr(lh);
      progress.Update();
      
      Ngs_Element ngel = ma->GetElement(elnr);
      // ELEMENT_TYPE eltype = ngel.GetType();
      Array<int> p1_dofs;
      lset_p1->GetFESpace()->GetDofNrs(elnr,p1_dofs);
      FlatVector<> vals(p1_dofs.Size(),lh);
      lset_p1->GetVector().GetIndirect(p1_dofs,vals);
      ElementTransformation & eltrans = ma->GetTrafo (ElementId(VOL,elnr), lh);
      
      if (!ElementInRelevantBand(vals, lower_lset_bound, upper_lset_bound))
        continue;


      Array<int> def_dofs;
      deform->GetFESpace()->GetDofNrs(elnr,def_dofs);
      const FiniteElement & fel_deform = deform->GetFESpace()->GetFE(elnr,lh);
      // FlatVector<> vals(def_dofs.Size(),lh);
      // lset_p1->GetVector().GetIndirect(def_dofs,vals);
      FlatMatrix<> massmat (def_dofs.Size(),lh);
      FlatVector<> elvec (D*def_dofs.Size(),lh);
      FlatVector<> elres (D*def_dofs.Size(),lh);
      mass->CalcElementMatrix(fel_deform, eltrans, massmat, lh);
      CalcInverse(massmat);

      
      Array<int> lset_ho_dofs;
      lset_ho->GetFESpace()->GetDofNrs(elnr,lset_ho_dofs);
      FlatVector<> lset_ho_vals(lset_ho_dofs.Size(),lh);
      lset_ho->GetVector().GetIndirect(lset_ho_dofs,lset_ho_vals);
      const FiniteElement & fel_lset_ho = lset_ho->GetFESpace()->GetFE(elnr,lh);
      
      if (D==2)
      {
        const ScalarFiniteElement<2> & scafe_lset_ho = dynamic_cast< const ScalarFiniteElement<2> &>(fel_lset_ho);
        shared_ptr<LsetEvaluator<2>> lseteval = make_shared<LsetEvaluator<2>>(scafe_lset_ho,lset_ho_vals);
        shift2D->CalcElementVector(fel_deform, eltrans, elvec, lh, lseteval);
        
        FlatMatrixFixWidth<2> elvec_vec(def_dofs.Size(),&elvec(0));
        FlatMatrixFixWidth<2> shift_vec(def_dofs.Size(),&elres(0));
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
        
        FlatMatrixFixWidth<3> elvec_vec(def_dofs.Size(),&elvec(0));
        FlatMatrixFixWidth<3> shift_vec(def_dofs.Size(),&elres(0));
        for (int d = 0; d < D; ++d)
          shift_vec.Col(d) = massmat * elvec_vec.Col(d);
        // vertex values to zero
        for (int l = 0; l < D+1; ++l)
          shift_vec.Row(l) = 0.0;
      }

      FlatVector<> def_vals(D*def_dofs.Size(),lh);
      deform->GetVector().GetIndirect(def_dofs,def_vals);
      def_vals += elres;
      deform->GetVector().SetIndirect(def_dofs,def_vals);

      // increase dof counter
      FlatVector<> factors(D*def_dofs.Size(),lh);
      factor->GetIndirect(def_dofs,factors);
      for (int k = 0; k < factors.Size(); k+=D)
        factors(k) += 1.0 ;
      factor->SetIndirect(def_dofs,factors);
      
    }
    progress.Done();
    
    // averaging of the (summed) deformation
    Array<int> dnums(1);
    for (int i = 0; i < factor->Size(); ++i)
    {
      FlatVector<> val_fac(D,lh);
      FlatVector<> values(D,lh);
      dnums[0] = i;
      deform->GetVector().GetIndirect(dnums,values);
      factor->GetIndirect(dnums,val_fac);
      if (val_fac(0) > 0)
        values *= 1.0/val_fac(0);
      deform->GetVector().SetIndirect(dnums,values);
    }

  }
    
  
  NumProcProjectShift::NumProcProjectShift (shared_ptr<PDE> apde, const Flags & flags)
  {
    lower_lset_bound = flags.GetNumFlag("lset_lower_bound",0.0);
    upper_lset_bound = flags.GetNumFlag("lset_upper_bound",0.0);
    threshold = flags.GetNumFlag("threshold",1.0);
    gf_lset_p1 = apde->GetGridFunction(flags.GetStringFlag("levelset_p1","gf_lset_p1"));
    lset = apde->GetGridFunction(flags.GetStringFlag("levelset","lset"));
    deform = apde->GetGridFunction(flags.GetStringFlag("deform","deform"));
    qn = apde->GetCoefficientFunction(flags.GetStringFlag("quasinormal","qn"));
  }

  void NumProcProjectShift::Do (LocalHeap & lh)
  {
    ProjectShift(lset, gf_lset_p1, deform, qn, lower_lset_bound, upper_lset_bound, threshold, lh);
  }

  static RegisterNumProc<NumProcProjectShift> npprojshift("projectshift");
}

