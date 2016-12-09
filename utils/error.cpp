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
  SolutionCoefficients<D>::SolutionCoefficients(shared_ptr<PDE> pde, const Flags & flags)
  {
    conv_n = pde->GetCoefficientFunction(flags.GetStringFlag("convection_n",""), true);
    conv_p = pde->GetCoefficientFunction(flags.GetStringFlag("convection_p",""), true);
    coef_n = pde->GetCoefficientFunction(flags.GetStringFlag("solution_n",""), true);
    coef_p = pde->GetCoefficientFunction(flags.GetStringFlag("solution_p",""), true);
    coef_d_n = pde->GetCoefficientFunction(flags.GetStringFlag("derivative_n",""), true);
    coef_d_p = pde->GetCoefficientFunction(flags.GetStringFlag("derivative_p",""), true);
    lset = pde->GetCoefficientFunction(flags.GetStringFlag("levelset",""), true);
    coef_jumprhs = pde->GetCoefficientFunction(flags.GetStringFlag("jumprhs",""), true);
  }
  
  template<int D>
  SolutionCoefficients<D>::~SolutionCoefficients()
  {
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
  void CalcXError (shared_ptr<GridFunction> gfu, 
                   shared_ptr<GridFunction> gfu2, 
                   SolutionCoefficients<D> & solcoef, 
                   int intorder, 
                   double a_neg, double a_pos, 
                   double b_neg, double b_pos, 
                   double time, 
                   ErrorTable & errtab, 
                   LocalHeap & clh,
                   bool output,
                   const Flags & flags)
  {
    static Timer time_fct ("CalcXError");
    RegionTimer reg (time_fct);
    
    double threshold  = flags.GetNumFlag ( "threshold", -1.0);

    // cout << " CalcXError at time = " << time << endl;
    // int activeels = 0;
    double l2diff_n = 0;
    double l2diff_p = 0;
    double l2diff = 0;
    double h1diff_n = 0;
    double h1diff_p = 0;
    double h1diff = 0;
    double ifjumpl2 = 0;
    double ifdudnl2 = 0;
    double ifsigmanl2 = 0;

    double mass_n = 0;
    double mass_p = 0;
    double mass_sol_n = 0;
    double mass_sol_p = 0;
    double vol = 0.0;
    
    shared_ptr<MeshAccess> ma (gfu->GetFESpace()->GetMeshAccess());

    ProgressOutput progress (ma, "CalcXError element", ma->GetNE());

    IterateElements 
      (*gfu->GetFESpace(), VOL, clh, 
       [&] (ElementId ei, LocalHeap & lh)
    {
         int elnr = ei.Nr();
         progress.Update ();
      HeapReset hr(lh);

         Array<int> dnums;
         Array<int> dnums2;
         
      gfu -> GetFESpace()->GetDofNrs (elnr, dnums);
      const int size = dnums.Size();

      FlatVector<double> elvec (size, lh);
      gfu -> GetVector().GetIndirect (dnums, elvec);   

      int size2 = 0;

      if (gfu2)
      {
        gfu2 -> GetFESpace()->GetDofNrs (elnr, dnums2);
        size2 = dnums2.Size();
      }

      FlatVector<double> elvec2 (size2, lh);

      if (gfu2)
      {
        gfu2 -> GetVector().GetIndirect (dnums2, elvec2);   
      }


      ElementTransformation & eltrans = ma->GetTrafo(elnr,false,lh);

      const FiniteElement & base_fel = gfu -> GetFESpace()->GetFE(elnr,lh);

      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (base_fel);

      const XFiniteElement * xfe = NULL;
      const XDummyFE * dummfe = NULL;
      const ScalarFiniteElement<D> * scafe = NULL;
      // const ScalarSpaceTimeFiniteElement<D> * scastfe = NULL;
      void * scastfe = NULL;

      for (int j = 0; j < cfel.GetNComponents(); ++j)
      {
        if (xfe==NULL)
          xfe = dynamic_cast<const XFiniteElement* >(&cfel[j]);
        if (dummfe==NULL)
          dummfe = dynamic_cast<const XDummyFE* >(&cfel[j]);
        if (scafe==NULL)
          scafe = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel[j]);
        // if (scastfe==NULL)
        //   scastfe = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel[j]);
      }

      const XFiniteElement * xfe2 = NULL;
      const XDummyFE * dummfe2 = NULL;
      const ScalarFiniteElement<D> * scafe2 = NULL;
      // const ScalarSpaceTimeFiniteElement<D> * scastfe2 = NULL;
      void * scastfe2 = NULL;


      if (gfu2)
      {
        const FiniteElement & base_fel2 = gfu2 -> GetFESpace()->GetFE(elnr,lh);

        const CompoundFiniteElement & cfel2 = 
          dynamic_cast<const CompoundFiniteElement&> (base_fel2);
        for (int j = 0; j < cfel2.GetNComponents(); ++j)
        {
          if (xfe2==NULL)
            xfe2 = dynamic_cast<const XFiniteElement* >(&cfel2[j]);
          if (dummfe2==NULL)
            dummfe2 = dynamic_cast<const XDummyFE* >(&cfel2[j]);
          if (scafe2==NULL)
            scafe2 = dynamic_cast<const ScalarFiniteElement<D>* >(&cfel2[j]);
          // if (scastfe2==NULL)
          //   scastfe2 = dynamic_cast<const ScalarSpaceTimeFiniteElement<D>* >(&cfel2[j]);
        }
      }

      bool spacetime = scastfe != NULL;

      int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
      // int ndof = spacetime ? scastfe->GetNDof() : scafe->GetNDof();
      int ndof = scafe->GetNDof();
      int ndof_total = ndof+ndof_x;
      FlatVector<> shape_total(ndof_total,lh);
      FlatVector<> shape(ndof,&shape_total(0));
      FlatVector<> shapex(ndof_x,&shape_total(ndof));

      FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
      FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
      FlatMatrixFixWidth<D> dshapex(ndof_x,&dshape_total(ndof,0));

      FlatVector<> dshapen_total(ndof_total,lh);
      FlatVector<> dshapen(ndof,&dshapen_total(0));
      FlatVector<> dshapenx(ndof_x,&dshapen_total(ndof));

      int ndof_x2 = xfe2!=NULL ? xfe2->GetNDof() : 0;
      // int ndof2 = gfu2 == NULL ? 0 : (spacetime ? scastfe2->GetNDof() : scafe2->GetNDof());
      int ndof2 = gfu2 == NULL ? 0 : scafe2->GetNDof();
      int ndof_total2 = gfu2 == NULL ? 0 : (ndof2+ndof_x2);
      FlatVector<> shape_total2(ndof_total2,lh);
      FlatVector<> shape2(ndof2,&shape_total2(0));
      FlatVector<> shapex2(ndof_x2,&shape_total2(ndof2));

      FlatMatrixFixWidth<D> dshape_total2(ndof_total2,lh);
      FlatMatrixFixWidth<D> dshape2(ndof2,&dshape_total2(0,0));
      FlatMatrixFixWidth<D> dshapex2(ndof_x2,&dshape_total2(ndof2,0));

      if (xfe)
      {
          
        DOMAIN_TYPE dt = POS;
        for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
        {
          // const FlatXLocalGeometryInformation & xgeom( spacetime ? 
          //                                              xfe->GetFlatLocalGeometryUpTrace() : 
          //                                              xfe->GetFlatLocalGeometry());
          const FlatXLocalGeometryInformation & xgeom( xfe->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
          for (int i = 0; i < fquad.Size(); ++i)
          {
            IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
            MappedIntegrationPoint<D,D> mip(ip, eltrans);
            DimMappedIntegrationPoint<D+1> mipp(ip, eltrans);
            mipp.Point().Range(0,D) = mip.GetPoint();
            mipp.Point()[D] = time;

            if (!spacetime)
            {
              shape = scafe->GetShape(mip.IP(), lh);
              scafe->CalcMappedDShape(mip, dshape);
              if (gfu2)
              {
                shape2 = scafe2->GetShape(mip.IP(), lh);
                scafe2->CalcMappedDShape(mip, dshape2);
              }
            }
            else
            {
              throw Exception("no spacetime anymore");
              // scastfe->CalcShapeSpaceTime(mip.IP(), 1.0, shape, lh);
              // scastfe->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape, lh);
              // if (gfu2)
              // {
              //   scastfe2->CalcShapeSpaceTime(mip.IP(), 1.0, shape2, lh);
              //   scastfe2->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape2, lh);
              // }
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

            double solval = 0.0;
            Vec<D> soldval;
            if (gfu2)
            {
              shapex2 = shape2;
              dshapex2 = dshape2;


              for (int l = 0; l < ndof_x2; ++l)
              {
                if (xfe2->GetSignsOfDof()[l] != dt)
                {
                  shapex2(l) = 0.0;
                  dshapex2.Row(l) = 0.0;
                }
              }
              solval = InnerProduct(shape_total2,elvec2);
              soldval = Trans(dshape_total2) * elvec2;

            }
            else
            {
              solval = dt == POS ? solcoef.GetSolutionPos().Evaluate(mipp) : solcoef.GetSolutionNeg().Evaluate(mipp);
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


            }

            double discval = InnerProduct(shape_total,elvec);
            Vec<D> discdval = Trans(dshape_total) * elvec;
            Vec<D> diffdval = soldval - discdval;
            Vec<D> conv;

            double diffdsqr = 0; 

            if (solcoef.HasConvectionNeg() && solcoef.HasConvectionPos())
            {
                if (dt == POS)
                  solcoef.GetConvectionPos().Evaluate(mip,conv);
                else
                  solcoef.GetConvectionNeg().Evaluate(mip,conv);
              diffdsqr = sqr(InnerProduct(diffdval,conv));
            }
            else
              diffdsqr = L2Norm2(diffdval);


            double fac = mip.GetWeight();
            if (dt == POS)
            {
              if (threshold <= 0.0)
              {
#pragma omp atomic
                l2diff_p += b_pos*fac*sqr(discval-solval);
#pragma omp atomic
                h1diff_p += b_pos*fac*diffdsqr;
              }
#pragma omp atomic
              mass_p += discval*fac;
#pragma omp atomic
                 mass_sol_p += solval*fac;
            }
            else
            {
              if (threshold <= 0.0)
              {
#pragma omp atomic
                l2diff_n += b_neg*fac*sqr(discval-solval);
#pragma omp atomic
                h1diff_n += b_neg*fac*diffdsqr;
              }
#pragma omp atomic
              mass_n += discval*fac;
#pragma omp atomic
                 mass_sol_n += solval*fac;
            }
#pragma omp atomic
               vol += fac;
          } // quad rule
        } // dt

        if (spacetime)
        {

        }
        else
        {

          IntegrationPoint ipc(0.0,0.0,0.0);
          MappedIntegrationPoint<D,D> mipc(ipc, eltrans);
          const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());

          FlatVector<> jump(ndof_total,lh);

          const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

          for (int i = 0; i < fquad.Size(); ++i)
          {
            IntegrationPoint ip(&fquad.points(i,0),0.0);
            MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
            Mat<D,D> Finv = mip.GetJacobianInverse();
            const double absdet = mip.GetMeasure();

            Vec<D> nref = fquad.normals.Row(i);
            Vec<D> normal = absdet * Trans(Finv) * nref ;
            double len = L2Norm(normal);
            normal /= len;

            const double weight = fquad.weights(i) * len; 
        
            shape = scafe->GetShape(mip.IP(), lh);
            jump.Range(0,ndof) = (b_pos-b_neg) * shape;
            jump.Range(ndof,ndof_total) = shape;

            for (int l = 0; l < ndof_x; ++l)
            {
              if (xfe->GetSignsOfDof()[l] == NEG)
                jump(ndof+l) *= -b_neg;
              else
                jump(ndof+l) *= b_pos;
            }
            
            double jumpval = InnerProduct(jump,elvec);
            if (solcoef.HasJumpRhs())
              jumpval -= solcoef.GetJumpRhs().Evaluate(mip);

#pragma omp atomic
            ifjumpl2 += weight * sqr(jumpval);


            if (!spacetime)
            {
              scafe->CalcMappedDShape(mip, dshape);
              if (gfu2)
                scafe2->CalcMappedDShape(mip, dshape2);
            }
            else
            {
              throw Exception("no space time");
              // scastfe->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape, lh);
              // if (gfu2)
              //   scastfe2->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape2, lh);
            }

            dshapex = dshape;

            double kappa_neg = xgeom.kappa[NEG];
            double kappa_pos = xgeom.kappa[POS];

            dshapen = (a_pos*kappa_pos+a_neg*kappa_neg) * (dshape * normal);
            dshapenx = dshape * normal;

            for (int l = 0; l < ndof_x; ++l)
              if (xfe->GetSignsOfDof()[l] == NEG)
                dshapenx(l) *= kappa_neg * a_neg;
              else
                dshapenx(l) *= kappa_pos * a_pos;


            double discdnval = InnerProduct(dshapen_total,elvec);

            Vec<D> soldvalneg;
            Vec<D> soldvalpos;
            if (gfu2)
            {
              cout << " no H^-1/2 - norm for gfu2 yet" << endl;
              // dshapex2 = dshape2;

              // for (int l = 0; l < ndof_x2; ++l)
              // {
              //   if (xfe2->GetSignsOfDof()[l] != dt)
              //   {
              //     dshapex2.Row(l) = 0.0;
              //   }
              // }
              // solval = InnerProduct(shape_total2,elvec2);
              // soldval = Trans(dshape_total2) * elvec2;

            }
            else
            {
              if (solcoef.HasSolutionDNeg() && solcoef.HasSolutionDPos())
              {
                // if (dt == POS)
                solcoef.GetSolutionDPos().Evaluate(mip,soldvalpos);
                // else
                solcoef.GetSolutionDNeg().Evaluate(mip,soldvalneg);
              }
              else
              {
                // if (dt == POS)
                CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionPos()),mip,time,soldvalpos,lh);
                // else
                CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionNeg()),mip,time,soldvalneg,lh);
              }
            }
            double soldnvalpos = a_pos * InnerProduct(soldvalpos,normal);
            double soldnvalneg = a_neg * InnerProduct(soldvalneg,normal);
            
            if (abs(soldnvalpos-soldnvalneg)/max(abs(soldnvalpos),abs(soldnvalneg)) > 1e-4)
            {
              static bool warndudn = false;
              if (!warndudn)
              {
                std::cout << " soldnvalpos = " << soldnvalpos << std::endl;
                std::cout << " soldnvalneg = " << soldnvalneg << std::endl;
                cout << "solutions neg/pos do not fit";
              }
              warndudn = true;
            }            
            
            double soldnval = kappa_neg * soldnvalneg + kappa_pos * soldnvalpos;
            
#pragma omp atomic
            ifdudnl2 += weight * sqr(soldnval - discdnval);

            const double ava = a_pos*kappa_pos+a_neg*kappa_neg;
            const double lam = 2;
            const double p = 1;
            // std::cout << " discdnval = " << discdnval << std::endl;
            discdnval -=  lam*(p+1)/p/h * ava * jumpval;
            // std::cout << " orig_jumpval = " << orig_jumpval << std::endl;
            // std::cout << " lam*(p+1)/p/h * ava * orig_jumpval = " << lam*(p+1)/p/h * ava * orig_jumpval << std::endl;

            // std::cout << " discdnval = " << discdnval << std::endl;
            // std::cout << " soldnval = " << soldnval << std::endl;

#pragma omp atomic
            ifsigmanl2 += weight * sqr(soldnval - discdnval);
          }
        }
      } // is xfe 
      else
      {
        DOMAIN_TYPE dt = dummfe->GetDomainType();


        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), intorder);

        bool skip = false;
        if (threshold > 0)
        {
          for (int i = 0 ; i < pir.GetNIP(); i++)
          {

            MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
            DimMappedIntegrationPoint<D+1> mipp(pir[i], eltrans);
            mipp.Point().Range(0,D) = mip.GetPoint();
            mipp.Point()[D] = time;
            if (threshold>0 && !solcoef.HasLevelSet())
              throw Exception("threshold>0 but no levelset...");
          
            const double lsetval = solcoef.GetLevelSet().Evaluate(mip);
            if (abs(lsetval) < threshold) skip = true;
          }
        }

        if (!skip)
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          DimMappedIntegrationPoint<D+1> mipp(pir[i], eltrans);
          mipp.Point().Range(0,D) = mip.GetPoint();
          mipp.Point()[D] = time;

          
          double solval = 0.0;
          Vec<D> soldval;
          if (gfu2)
          {
            if (!spacetime)
            {
              shape2 = scafe2->GetShape(mip.IP(), lh);
              scafe2->CalcMappedDShape(mip, dshape2);
            }
            else
            {
              throw Exception("no spacetime");
              // scastfe2->CalcShapeSpaceTime(mip.IP(), 1.0, shape2, lh);
              // scastfe2->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape2, lh);
            }

            solval = InnerProduct(shape_total2,elvec2);
            soldval = Trans(dshape_total2) * elvec2;

          }
          else
          {
            solval = dt == POS ? solcoef.GetSolutionPos().Evaluate(mipp) : solcoef.GetSolutionNeg().Evaluate(mipp);
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
          }

          if (!spacetime)
          {
            shape = scafe->GetShape(mip.IP(), lh);
            scafe->CalcMappedDShape(mip, dshape);
          }
          else
          {
            throw Exception("no spacetime");
            // scastfe->CalcShapeSpaceTime(mip.IP(), 1.0, shape, lh);
            // scastfe->CalcMappedDxShapeSpaceTime(mip, 1.0, dshape, lh);
          }

          double discval = InnerProduct(shape,elvec);
          Vec<D> discdval = Trans(dshape) * elvec;
          Vec<D> diffdval = soldval - discdval;

          Vec<D> conv;

          double diffdsqr = 0; 

          if (solcoef.HasConvectionNeg() && solcoef.HasConvectionPos())
          {
            if (dt == POS)
              solcoef.GetConvectionPos().Evaluate(mip,conv);
            else
              solcoef.GetConvectionNeg().Evaluate(mip,conv);
            diffdsqr = sqr(InnerProduct(diffdval,conv));
          }
          else
            diffdsqr = L2Norm2(diffdval);

          double fac = mip.GetWeight();
          if (dt == POS)
          {
#pragma omp atomic
            l2diff_p += b_pos*fac*sqr(discval-solval);
#pragma omp atomic
            h1diff_p += b_pos*fac*diffdsqr;
#pragma omp atomic
            mass_p += discval*fac;
#pragma omp atomic
                 mass_sol_p += solval*fac;
          }
          else
          {
            // cout << "discval, solval: " << discval << ", " << solval << endl;
#pragma omp atomic
            l2diff_n += b_neg*fac*sqr(discval-solval);
#pragma omp atomic
            h1diff_n += b_neg*fac*diffdsqr;
#pragma omp atomic
            mass_n += discval*fac;
#pragma omp atomic
                 mass_sol_n += solval*fac;
          }
#pragma omp atomic
               vol += fac;
        }
      }
       });

    ifjumpl2 = sqrt(ifjumpl2); errtab.iferr.Append(ifjumpl2);
    ifdudnl2 = sqrt(ifdudnl2); errtab.ifdudnerr.Append(ifdudnl2);
    ifsigmanl2 = sqrt(ifsigmanl2); errtab.ifsigmanerr.Append(ifsigmanl2);
    l2diff = b_pos * l2diff_p + b_neg * l2diff_n;
    h1diff = b_pos * h1diff_p + b_neg * h1diff_n;
    l2diff_p = sqrt(l2diff_p); errtab.l2err_p.Append(l2diff_p);
    h1diff_p = sqrt(h1diff_p); errtab.h1err_p.Append(h1diff_p);
    l2diff_n = sqrt(l2diff_n); errtab.l2err_n.Append(l2diff_n);
    h1diff_n = sqrt(h1diff_n); errtab.h1err_n.Append(h1diff_n);
    l2diff = sqrt(l2diff); errtab.l2err.Append(l2diff);
    h1diff = sqrt(h1diff); errtab.h1err.Append(h1diff);
    // cout << " activeels = " << activeels << endl;
    if (!output) //hack for stokes
    {
      shared_ptr<BaseVector> corr = gfu->GetComponent(0)->GetVector().CreateVector();
      corr->FVDouble() = (mass_sol_n + mass_sol_p - mass_n - mass_p)/vol;
      cout << " vol = " << vol << endl;
      cout << " correction = " << (mass_sol_n + mass_sol_p - mass_n - mass_p)/vol << endl;
      gfu->GetComponent(0)->GetVector() += *corr;
    }
    if (output)
    {
      cout << endl;
      cout << " mass_n = " << mass_n << endl;
      cout << " mass_p = " << mass_p << endl;
      cout << " mass_sol_n = " << mass_sol_n << endl;
      cout << " mass_sol_p = " << mass_sol_p << endl;
      cout << " total mass = " << mass_p + mass_n << endl;
      cout << " total mass_sol = " << mass_sol_p + mass_sol_n << endl;

      cout << endl;
      cout << setw(12) << "l2_n" << "       |";
      cout << setw(12) << "l2_p" << "       |";
      cout << setw(12) << "l2" << "       |";
      cout << setw(12) << "h1_n" << "       |";
      cout << setw(12) << "h1_p" << "       |";
      cout << setw(12) << "h1" << "       |";
      cout << setw(12) << "ifl2" << "       |";
      cout << setw(12) << "ifdudnl2" << "        |";
      cout << setw(12) << "ifsigml2" << endl;
      for (int i = 0; i < 174; ++i)
        cout << "-";
      cout << endl;

      std::streamsize p = cout.precision();
      const int N = errtab.l2err_n.Size();
      for (int i = 0; i < N; ++i)
      {
        cout << setw(12) << errtab.l2err_n[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.l2err_n[i]/errtab.l2err_n[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.l2err_p[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.l2err_p[i]/errtab.l2err_p[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.l2err[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.l2err[i]/errtab.l2err[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.h1err_n[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.h1err_n[i]/errtab.h1err_n[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.h1err_p[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.h1err_p[i]/errtab.h1err_p[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.h1err[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.h1err[i]/errtab.h1err[i-1])/log(0.5) << "]|";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.iferr[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "      |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.iferr[i]/errtab.iferr[i-1])/log(0.5) << "]| ";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.ifdudnerr[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "       |";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.ifdudnerr[i]/errtab.ifdudnerr[i-1])/log(0.5) << "]| ";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << setw(12) << errtab.ifsigmanerr[i] << " ";
        cout.setf(ios::fixed);
        cout.precision(2);
        if (N>0 && i==0)
          cout << "       ";
        if (i>0)
          cout << "[" << setw(4) << log(errtab.ifsigmanerr[i]/errtab.ifsigmanerr[i-1])/log(0.5) << "] ";
        cout.precision(p);
        cout.unsetf(ios::fixed);

        cout << endl;
      }
    }
  }

  template void CalcXError<2>(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfu2, SolutionCoefficients<2> & solcoef, int intorder, double a_neg, double a_pos,
                              double b_neg, double b_pos, double time, ErrorTable & errtab, LocalHeap & lh, bool output, const Flags & flags);
  template void CalcXError<3>(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfu2, SolutionCoefficients<3> & solcoef, int intorder, double a_neg, double a_pos,
                              double b_neg, double b_pos, double time, ErrorTable & errtab, LocalHeap & lh, bool output, const Flags & flags);


  template<int D>
  void CalcTraceDiff (shared_ptr<GridFunction> gfu, 
                      shared_ptr<CoefficientFunction> coef, 
                      int intorder, Array<double> & errors,
                      LocalHeap & clh)
  {
    static Timer time_fct ("CalcTraceDiff");
    RegionTimer reg (time_fct);
    
    double l2diff = 0;
    double maxdiff = 0;
    // double h1diff = 0;
    // double surf = 0.0;
    
    shared_ptr<MeshAccess> ma (gfu->GetFESpace()->GetMeshAccess());

    ProgressOutput progress (ma, "CalcTraceDiff element", ma->GetNE());
    
    IterateElements 
      (*gfu->GetFESpace(), VOL, clh, 
       [&] (ElementId ei, LocalHeap & lh)
       {
         int elnr = ei.Nr();
         progress.Update ();
         HeapReset hr(lh);

         Array<int> dnums;
         gfu -> GetFESpace()->GetDofNrs (elnr, dnums);
         const int ndof = dnums.Size();
         FlatVector<double> elvec (ndof, lh);
         gfu -> GetVector().GetIndirect (dnums, elvec);   

         ElementTransformation & eltrans = ma->GetTrafo(ei,lh);

         const XFiniteElement * xfe
           = dynamic_cast<const XFiniteElement* >(&(gfu->GetFESpace()->GetFE(elnr,lh)) );

         if (xfe)
         {
           FlatVector<> shape(ndof,lh);

           const ScalarFiniteElement<D> & scafe =
             dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());
           IntegrationPoint ipc(0.0,0.0,0.0);
           MappedIntegrationPoint<D,D> mipc(ipc, eltrans);
           // const double h = D == 2 ? sqrt(mipc.GetMeasure()) : cbrt(mipc.GetMeasure());
//              FlatVector<> jump(ndof_total,lh);

           const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
           const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
           const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

             for (int i = 0; i < fquad.Size(); ++i)
             {
               IntegrationPoint ip(&fquad.points(i,0),0.0);
               MappedIntegrationPoint<D,D> mip(ip, eltrans);
      
               Mat<D,D> Finv = mip.GetJacobianInverse();
               const double absdet = mip.GetMeasure();

               Vec<D> nref = fquad.normals.Row(i);
               Vec<D> normal = absdet * Trans(Finv) * nref ;
               double len = L2Norm(normal);
               normal /= len;

               const double weight = fquad.weights(i) * len; 

// #pragma omp atomic
               // surf += weight;
               shape = scafe.GetShape(mip.IP(), lh);
               const double discrete_val = InnerProduct(shape,elvec);
               const double sol_val = coef->Evaluate(mip);

               const double error_contrib_l2 = sqr(discrete_val - sol_val) * weight;
#pragma omp atomic
               l2diff += error_contrib_l2;
               const double curr_maxdiff = abs(discrete_val - sol_val);
#pragma omp critical (max)
               {
                 if (curr_maxdiff > maxdiff)
                   maxdiff = curr_maxdiff;
               }
               
             }
         } // is xfe 
       }); //iterate elements end
    l2diff = sqrt(l2diff);
    // cout << " surf = " << surf << endl;
    errors.SetSize(2);
    errors[0] = l2diff;
    errors[1] = maxdiff;
  }

  template void CalcTraceDiff<2>(shared_ptr<GridFunction> gfu, shared_ptr<CoefficientFunction> coef, int intorder, Array<double>& errors, LocalHeap & lh);
  template void CalcTraceDiff<3>(shared_ptr<GridFunction> gfu, shared_ptr<CoefficientFunction> coef, int intorder, Array<double>& errors, LocalHeap & lh);


}
