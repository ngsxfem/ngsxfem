#include "error.hpp"

namespace ngfem
{
  template<int D>
  void CalcDxShapeOfCoeff(const CoefficientFunction * coef, const MappedIntegrationPoint<D,D>& mip, 
                          double time,
                          Vec<D>& der, LocalHeap& lh)
  {
    HeapReset hr(lh);
    // bmatu = 0;
    // evaluate dshape by numerical diff
    //fel_u, eltrans, sip, returnval, lh

    const IntegrationPoint& ip = mip.IP();//volume_ir[i];
    const ElementTransformation & eltrans = mip.GetTransformation();
    
    Vec<D> der_ref;
  
    double eps = 1e-7;
    for (int j = 0; j < D; j++)   // d / dxj
	{
	  IntegrationPoint ipl(ip);
	  ipl(j) -= eps;
	  MappedIntegrationPoint<D,D> sipl(ipl, eltrans);
      DimMappedIntegrationPoint<D+1> mipl(ipl, eltrans);
      mipl.Point().Range(0,D) = sipl.GetPoint();
      mipl.Point()[D] = time;

      IntegrationPoint ipr(ip);
      ipr(j) += eps;
      MappedIntegrationPoint<D,D> sipr(ipr, eltrans);
      DimMappedIntegrationPoint<D+1> mipr(ipr, eltrans);
      mipr.Point().Range(0,D) = sipr.GetPoint();
      mipr.Point()[D] = time;

      const double valright = coef->Evaluate(mipr);
      const double valleft = coef->Evaluate(mipl);
      
	  der_ref[j] = (1.0/(2*eps)) * (valright-valleft);
	}
                              
    der = Trans(mip.GetJacobianInverse()) * der_ref;
  }

  // template<int D>
  // SolutionCoefficients<D>::SolutionCoefficients(const CoefficientFunction * a_coef_n,
  //                                               const CoefficientFunction * a_coef_p,
  //                                               const CoefficientFunction * a_coef_d_n,
  //                                               const CoefficientFunction * a_coef_d_p,
  //                                               const CoefficientFunction * a_lset)
  //   : coef_n(a_coef_n),
  //     coef_p(a_coef_p),
  //     coef_d_n(a_coef_d_n),
  //     coef_d_p(a_coef_d_p),
  //     lset(a_lset)
  // {; }

  template<int D>
  SolutionCoefficients<D>::SolutionCoefficients(PDE & pde, const Flags & flags)
  {
    made_coef_n = MakeHigherDimensionCoefficientFunction<D>(
      pde.GetCoefficientFunction(flags.GetStringFlag("function_n","")),
      coef_n);
    made_coef_p = MakeHigherDimensionCoefficientFunction<D>(
      pde.GetCoefficientFunction(flags.GetStringFlag("function_p","")),
      coef_p);
    made_coef_d_n = MakeHigherDimensionCoefficientFunction<D>(
      pde.GetCoefficientFunction(flags.GetStringFlag("derivative_n",""),true),
      coef_d_n);
    made_coef_d_p = MakeHigherDimensionCoefficientFunction<D>(
      pde.GetCoefficientFunction(flags.GetStringFlag("derivative_p",""),true),
      coef_d_p);
    made_lset = MakeHigherDimensionCoefficientFunction<D>(
      pde.GetCoefficientFunction(flags.GetStringFlag("levelset","")),
      lset);
  }
  
  template<int D>
  SolutionCoefficients<D>::~SolutionCoefficients()
  {
    if (made_coef_n) delete &(GetSolutionNeg());
    if (made_coef_p) delete &(GetSolutionPos());
    if (made_coef_d_n) delete &(GetSolutionDNeg());
    if (made_coef_d_p) delete &(GetSolutionDPos());
    if (made_lset) delete &(GetLevelSet());
  }

  template class SolutionCoefficients<2>;
  template class SolutionCoefficients<3>;

  ErrorTable::ErrorTable()
  {
    l2err_p.SetSize(0);
    l2err_n.SetSize(0);
    l2err.SetSize(0);
    h1err_p.SetSize(0);
    h1err_n.SetSize(0);
    h1err.SetSize(0);
    iferr.SetSize(0);
  }

  template<int D>
  void CalcXError (GridFunction * gfu, SolutionCoefficients<D> & solcoef, double intorder, 
                   double b_neg, double b_pos, double time, ErrorTable & errtab, LocalHeap & lh)
  {
    Array<int> dnums;
    int activeels = 0;
    double l2diff_n = 0;
    double l2diff_p = 0;
    double l2diff = 0;
    double h1diff_n = 0;
    double h1diff_p = 0;
    double h1diff = 0;
    double ifjumpl2 = 0;
    const MeshAccess & ma (gfu->GetFESpace().GetMeshAccess());
    for (int i = 0; i < ma.GetNE(); ++i)
    {
      HeapReset hr(lh);
      gfu -> GetFESpace().GetDofNrs (i, dnums);
      const int size = dnums.Size();

      FlatVector<double> elvec (size, lh);
      gfu -> GetVector().GetIndirect (dnums, elvec);   

      ElementTransformation & eltrans = ma.GetTrafo(i,false,lh);

      const FiniteElement & base_fel = gfu -> GetFESpace().GetFE(i,lh);
      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (base_fel);

      const XFiniteElement * xfe = NULL;
      const XDummyFE * dummfe = NULL;
      const ScalarFiniteElement<D> * scafe = NULL;
      const ScalarSpaceTimeFiniteElement<D> * scastfe = NULL;

      for (int j = 0; j < cfel.GetNComponents(); ++j)
      {
        if (xfe==NULL)
          xfe = dynamic_cast<const XFiniteElement* >(&cfel[j]);
        if (dummfe==NULL)
          dummfe = dynamic_cast<const XDummyFE* >(&cfel[j]);
        if (scafe==NULL)
          scafe = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel[j]);
        if (scastfe==NULL)
          scastfe = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel[j]);
      }

      bool spacetime = scastfe != NULL;

      int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
      int ndof = spacetime ? scastfe->GetNDof() : scafe->GetNDof();
      int ndof_total = ndof+ndof_x;
      FlatVector<> shape_total(ndof_total,lh);
      FlatVector<> shape(ndof,&shape_total(0));
      FlatVector<> shapex(ndof,&shape_total(ndof));

      FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
      FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
      FlatMatrixFixWidth<D> dshapex(ndof,&dshape_total(ndof,0));

      if (xfe)
      {

        DOMAIN_TYPE dt = POS;
        for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
        {

          const FlatXLocalGeometryInformation & xgeom( spacetime ? 
                                                       xfe->GetFlatLocalGeometryUpTrace() : 
                                                       xfe->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
          for (int i = 0; i < fquad.Size(); ++i)
          {
            IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
            MappedIntegrationPoint<D,D> mip(ip, eltrans);
            DimMappedIntegrationPoint<D+1> mipp(ip, eltrans);
            mipp.Point().Range(0,D) = mip.GetPoint();
            mipp.Point()[D] = time;

            double solval = dt == POS ? solcoef.GetSolutionPos().Evaluate(mipp) : solcoef.GetSolutionNeg().Evaluate(mipp);

            Vec<D> soldval;
            if (solcoef.HasSolutionDNeg() && solcoef.HasSolutionDPos())
            {
              if (dt == POS)
                solcoef.GetSolutionDPos().Evaluate(mipp,soldval);
              else
                solcoef.GetSolutionDNeg().Evaluate(mipp,soldval);
            }
            else
            {
              if (dt == POS)
                CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionPos()),mip,time,soldval,lh);
              else
                CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionNeg()),mip,time,soldval,lh);
            }

            if (!spacetime)
            {
              shape = scafe->GetShape(mip.IP(), lh);
              scafe->CalcMappedDShape(mip, dshape);
            }
            else
            {
              scastfe->CalcShapeSpaceTime(mip.IP(), 1.0, shape, lh);
              scastfe->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape, lh);
            }

            shapex = shape;
            dshapex = dshape;

            for (int l = 0; l < ndof_x; ++l)
            {
              if (xfe->GetSignsOfDof()[l] != dt)
              {
                shapex(l) = 0.0;
                dshapex.Row(l) = 0.0;
              }
            }

            double discval = InnerProduct(shape_total,elvec);
            Vec<D> discdval = Trans(dshape_total) * elvec;
            Vec<D> diffdval = soldval - discdval;
            double diffdsqr = L2Norm2(diffdval);

            double fac = mip.GetWeight();
            if (dt == POS)
            {
              l2diff_p += b_pos*fac*sqr(discval-solval);
              h1diff_p += b_pos*fac*diffdsqr;
            }
            else
            {
              l2diff_n += b_neg*fac*sqr(discval-solval);
              h1diff_n += b_neg*fac*diffdsqr;
            }
          } // quad rule

        }
      } // loop over els.
      else
      {
        DOMAIN_TYPE dt = dummfe->GetDomainType();


        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), intorder);
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          DimMappedIntegrationPoint<D+1> mipp(pir[i], eltrans);
          mipp.Point().Range(0,D) = mip.GetPoint();
          mipp.Point()[D] = time;

          double solval = dt == POS ? solcoef.GetSolutionPos().Evaluate(mipp) : solcoef.GetSolutionNeg().Evaluate(mipp);

          Vec<D> soldval;
          if (solcoef.HasSolutionDNeg() && solcoef.HasSolutionDPos())
          {
            if (dt == POS)
              solcoef.GetSolutionDPos().Evaluate(mipp,soldval);
            else
              solcoef.GetSolutionDNeg().Evaluate(mipp,soldval);
          }
          else
          {
            if (dt == POS)
              CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionPos()),mip,time,soldval,lh);
            else
              CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionNeg()),mip,time,soldval,lh);
          }
          if (!spacetime)
          {
            shape = scafe->GetShape(mip.IP(), lh);
            scafe->CalcMappedDShape(mip, dshape);
          }
          else
          {
            scastfe->CalcShapeSpaceTime(mip.IP(), 1.0, shape, lh);
            scastfe->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape, lh);
          }

          double discval = InnerProduct(shape,elvec);
          Vec<D> discdval = Trans(dshape) * elvec;
          Vec<D> diffdval = soldval - discdval;
          double diffdsqr = L2Norm2(diffdval);

          double fac = mip.GetWeight();
          if (dt == POS)
          {
            l2diff_p += b_pos*fac*sqr(discval-solval);
            h1diff_p += b_pos*fac*diffdsqr;
          }
          else
          {
            // cout << "discval, solval: " << discval << ", " << solval << endl;
            l2diff_n += b_neg*fac*sqr(discval-solval);
            h1diff_n += b_neg*fac*diffdsqr;
          }
        }
      }
    }

    ifjumpl2 = sqrt(ifjumpl2); errtab.iferr.Append(ifjumpl2);
    l2diff = l2diff_p + l2diff_n;
    h1diff = h1diff_p + h1diff_n;
    l2diff_p = sqrt(l2diff_p); errtab.l2err_p.Append(l2diff_p);
    h1diff_p = sqrt(h1diff_p); errtab.h1err_p.Append(h1diff_p);
    l2diff_n = sqrt(l2diff_n); errtab.l2err_n.Append(l2diff_n);
    h1diff_n = sqrt(h1diff_n); errtab.h1err_n.Append(h1diff_n);
    l2diff = sqrt(l2diff); errtab.l2err.Append(l2diff);
    h1diff = sqrt(h1diff); errtab.h1err.Append(h1diff);
    cout << " activeels = " << activeels << endl;
    cout << setw(12) << "l2_n" << "\t|";
    cout << setw(12) << "l2_p" << "\t|";
    cout << setw(12) << "l2" << "\t|";
    cout << setw(12) << "h1_n" << "\t|";
    cout << setw(12) << "h1_p" << "\t|";
    cout << setw(12) << "h1" << "\t|";
    cout << setw(12) << "ifl2" << endl;

    for (int i = 0; i < errtab.l2err_n.Size(); ++i)
    {
      cout << setw(12) << errtab.l2err_n[i] << "\t|";
      cout << setw(12) << errtab.l2err_p[i] << "\t|";
      cout << setw(12) << errtab.l2err[i] << "\t|";
      cout << setw(12) << errtab.h1err_n[i] << "\t|";
      cout << setw(12) << errtab.h1err_p[i] << "\t|";
      cout << setw(12) << errtab.h1err[i] << "\t|";
      cout << setw(12) << errtab.iferr[i] << endl;
    }
  }

  template void CalcXError<2>(GridFunction * gfu, SolutionCoefficients<2> & solcoef, double intorder, 
                              double b_neg, double b_pos, double time, ErrorTable & errtab, LocalHeap & lh);
  template void CalcXError<3>(GridFunction * gfu, SolutionCoefficients<3> & solcoef, double intorder, 
                              double b_neg, double b_pos, double time, ErrorTable & errtab, LocalHeap & lh);


}
