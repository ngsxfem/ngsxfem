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

    // int dimflux = 1; //diffop ? diffop->Dim() : bli.DimFlux(); 
    
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

         // for (int i = 0; i < elmat.Width(); ++i)
         //   elmat(i,i) += 1e-8;

         fes.TransformMat (ei.Nr(), bound, elmat, TRANSFORM_MAT_LEFT_RIGHT);
         fes.TransformVec (ei.Nr(), bound, elvec, TRANSFORM_RHS);

         bool regularize = false;
         for (int i = 0; i < dnums.Size(); ++i)
           if (abs(elmat(i,i)) < 1e-13)
             regularize = true;

         if (regularize)
         {
           for (int i = 0; i < dnums.Size(); ++i)
             for (int j = dnums.Size()/2; j < dnums.Size(); ++j)
             {
               elmat(i,j) = 0.0;
               elmat(j,i) = 0.0;
             }
           for (int j = dnums.Size()/2; j < dnums.Size(); ++j)
             elmat(j,j) = 1.0;
         }
              
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

         if (regularize)
         {
           const CompoundFiniteElement & cfel = 
             dynamic_cast<const CompoundFiniteElement&> (bfel);

           const XFiniteElement * xfe = NULL;

           if (xfe==NULL)
             xfe = dynamic_cast<const XFiniteElement* >(&cfel[1]);
           
           const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
           // const FlatArray<DOMAIN_TYPE>& xsign = xfe->GetSignsOfDof();

           // DOMAIN_TYPE dt = xgeom.kappa[NEG] < 1e-12 ? POS : NEG;
           std::cout << " xgeom.kappa[NEG] = " << xgeom.kappa[NEG] << std::endl;
           std::cout << " xgeom.kappa[POS] = " << xgeom.kappa[POS] << std::endl;
           // for (int i = 0; i < dnums.Size()/2; ++i)
           // {
           //   if(xsign[i] == dt)
           //   {
           //     cnti[dnums[i]]--;
           //     elveci(dnums.Size()/2+i) = elveci(i);
           //     elveci(i) = 0.0;
           //   }
           //   else
           //     elveci(dnums.Size()/2+i) = 0.0;
           // }
         }

         bool printmore = regularize;
         // for (int i = 0; i < elveci.Size(); ++i)
         //   if (abs(elveci(i)) > 1.21)
         //     printmore = true;

         if (printmore)
         {
           std::cout << " elveci = " << elveci << std::endl;
           std::cout << " elmat = " << elmat << std::endl;
           std::cout << " elvec = " << elvec << std::endl;
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





/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class NumProcXDifference : public NumProc
  {
  protected:
    GridFunction * gfu;
    CoefficientFunction * coef_n;
    CoefficientFunction * coef_p;
    CoefficientFunction * coef_d_n;
    CoefficientFunction * coef_d_p;
    CoefficientFunction * lset;
    double threshold;
    int intorder;
    double b_pos;
    double b_neg;
    Array<double> l2err_p;
    Array<double> l2err_n;
    Array<double> l2err;
    Array<double> h1xerr_p;
    Array<double> h1xerr_n;
    Array<double> h1yerr_p;
    Array<double> h1yerr_n;
    Array<double> h1xerr;
    Array<double> h1yerr;
    Array<double> h1err;
    Array<double> iferr;
    double time;
  public:
    
    NumProcXDifference (PDE & apde, const Flags & flags)
      : NumProc (apde)
    { 
      gfu  = pde.GetGridFunction (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));
      coef_n= pde.GetCoefficientFunction(flags.GetStringFlag("function_n",""));
      coef_p= pde.GetCoefficientFunction(flags.GetStringFlag("function_p",""));
      coef_d_n= pde.GetCoefficientFunction(flags.GetStringFlag("derivative_n",""));
      coef_d_p= pde.GetCoefficientFunction(flags.GetStringFlag("derivative_p",""));
      lset= pde.GetCoefficientFunction(flags.GetStringFlag("levelset",""));
      threshold = flags.GetNumFlag ( "threshold", -0.1);
      intorder = (int) flags.GetNumFlag ( "intorder", 2);
      b_pos = flags.GetNumFlag ( "henryweight_p", 1.0);
      b_neg = flags.GetNumFlag ( "henryweight_n", 1.0);
      time = flags.GetNumFlag ( "time", 1.0);
      l2err_p.SetSize(0);
      l2err_n.SetSize(0);
      l2err.SetSize(0);
      h1xerr_p.SetSize(0);
      h1xerr_n.SetSize(0);
      h1yerr_p.SetSize(0);
      h1yerr_n.SetSize(0);
      h1xerr.SetSize(0);
      h1yerr.SetSize(0);
      h1err.SetSize(0);
      iferr.SetSize(0);
    }
  
    virtual string GetClassName () const
    {
      return "NumProcXDifference";
    }


    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      refinements++;
      cout << " This is the Do-call on refinement level " << refinements << std::endl;
      Array<int> dnums;
      int activeels = 0;
      double l2diff_n = 0;
      double l2diff_p = 0;
      double l2diff = 0;
      double h1xdiff_n = 0;
      double h1xdiff_p = 0;
      double h1ydiff_n = 0;
      double h1ydiff_p = 0;
      double h1xdiff = 0;
      double h1ydiff = 0;
      double h1diff = 0;
      double ifjumpl2 = 0;
      for (int i = 0; i < pde.GetMeshAccess().GetNE(); ++i)
      {
        HeapReset hr(lh);
        gfu -> GetFESpace().GetDofNrs (i, dnums);
        const int size = dnums.Size();

        FlatVector<double> elvec (size, lh);
        gfu -> GetVector().GetIndirect (dnums, elvec);   

        // ElementTransformation eltrans;
        // ma.GetElementTransformation (i, eltrans, lh);

        ElementTransformation & eltrans = ma.GetTrafo(i,false,lh);

        // IntegrationRule pir;
        // pir = SelectIntegrationRule (eltrans.GetElementType(), intorder);
                  
        // bool skip = false;
        // for (int j = 0 ; j < pir.GetNIP(); j++)
        // {
        //   MappedIntegrationPoint<D,D> mip(pir[j], eltrans);
        //   if (abs(lset->Evaluate(mip)) <= threshold)
        //   {
        //     skip = true;
        //     break;
        //   }
        // }

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

              double solval = dt == POS ? coef_p->Evaluate(mipp) : coef_n->Evaluate(mipp);

              Vec<D> soldval;
              if (dt == POS)
                coef_d_p->Evaluate(mipp,soldval);
              else
                coef_d_n->Evaluate(mipp,soldval);
          
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
                h1xdiff_p += b_pos*fac*diffdsqr;
              }
              else
              {
                l2diff_n += b_neg*fac*sqr(discval-solval);
                h1xdiff_n += b_neg*fac*diffdsqr;
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
            double solval = dt == POS ? coef_p->Evaluate(mipp) : coef_n->Evaluate(mipp);

            Vec<D> soldval;
            if (dt == POS)
              coef_d_p->Evaluate(mipp,soldval);
            else
              coef_d_n->Evaluate(mipp,soldval);

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
              h1xdiff_p += b_pos*fac*diffdsqr;
            }
            else
            {
              // cout << "discval, solval: " << discval << ", " << solval << endl;
              l2diff_n += b_neg*fac*sqr(discval-solval);
              h1xdiff_n += b_neg*fac*diffdsqr;
            }
          }
          

        }
            /*
            if (threshold<=0.0)
            {
              activeels++;
              bool sign=false;
              FlatMatrixFixWidth<D> dshape_h1(ndof_h1,&dshape(0,0)); //flat overlay
              FlatMatrixFixWidth<D> dshape_x(ndof_x,&dshape(ndof_h1,0)); //flat overlay
              int p = scafe->Order();

              for (int k=0; k<2; sign=!sign, k++)
              {
                CoefficientFunction * sol = sign ? coef_p : coef_n;
                CoefficientFunction * derx = sign ? coef_dx_p : coef_dx_n;
                CoefficientFunction * dery = sign ? coef_dy_p : coef_dy_n;

                IntegrationRule pir; //partial integration rule
                if(!xfe)
                { 
                  if (dummfe->GetSign() != sign)
                    continue;
                  pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
                }
                else{
                  xfe->GetMasterElement().FillVolumeIntegrationRule(intorder,sign,pir);
                }

                for (int i = 0 ; i < pir.GetNIP(); i++)
                {
                  MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
                  scafe->CalcMappedDShape (mip,dshape_h1);
                  shape.Range(0,ndof_h1) = scafe->GetShape(mip.IP(), lh);
                  if (xfe)
                  {
                    xfe->CalcMappedDShape (mip,dshape_x);
                    shape.Range(ndof_h1,ndof) = xfe->GetShape(mip.IP(), lh);
                    for (int l = 0; l < ndof_x; ++l)
                    {
                      if (xfe->GetSignsOfDof()[l] == sign){
                        dshape.Row(ndof_h1+l) = 0.0;
                        shape(ndof_h1+l) = 0.0;
                      }
                    }
                  }

                  const double solval = sol->Evaluate(mip);
                  const double discval = InnerProduct(elvec,shape);
                  // const double absdet = mip.GetJacobiDet();
                  const double derxval = derx->Evaluate(mip);
                  const double deryval = dery->Evaluate(mip);
                  const double weight = mip.GetWeight();
                  Vec<2> gradu = Trans(dshape) * elvec;
                  if (sign)
                    l2diff_p += (solval-discval)*(solval-discval)*weight*b_pos;
                  else
                    l2diff_n += (solval-discval)*(solval-discval)*weight*b_neg;

                  if (sign){
                    h1xdiff_p += (derxval-gradu(0))*(derxval-gradu(0))*weight*b_pos;
                    h1ydiff_p += (deryval-gradu(1))*(deryval-gradu(1))*weight*b_pos;
                  }
                  else
                  {
                    h1xdiff_n += (derxval-gradu(0))*(derxval-gradu(0))*weight*b_pos;
                    h1ydiff_n += (deryval-gradu(1))*(deryval-gradu(1))*weight*b_pos;
                  }
                }
                
              }
            }
            */
            
          }
          /*
          else if (!skip)
          {
            activeels++;
            if (!dummfe) throw Exception(" not an extended fespace?! ");

            bool sign = dummfe->GetSign();
            CoefficientFunction * sol = sign ? coef_p : coef_n;
            CoefficientFunction * derx = sign ? coef_dx_p : coef_dx_n;
            CoefficientFunction * dery = sign ? coef_dy_p : coef_dy_n;

            for (int j = 0 ; j < pir.GetNIP(); j++)
            {
              MappedIntegrationPoint<D,D> mip(pir[j], eltrans);
              const double solval = sol->Evaluate(mip);
              const double discval = InnerProduct(elvec,scafe->GetShape(mip.IP(), lh));
              // const double absdet = mip.GetJacobiDet();
              scafe->CalcMappedDShape (mip,dshape);
              const double derxval = derx->Evaluate(mip);
              const double deryval = dery->Evaluate(mip);
              const double weight = mip.GetWeight();
              Vec<2> gradu = Trans(dshape) * elvec;
              if (sign)
                l2diff_p += (solval-discval)*(solval-discval)*weight*b_pos;
              else
                l2diff_n += (solval-discval)*(solval-discval)*weight*b_neg;

                if (sign){
                h1xdiff_p += (derxval-gradu(0))*(derxval-gradu(0))*weight*b_pos;
                h1ydiff_p += (deryval-gradu(1))*(deryval-gradu(1))*weight*b_pos;
                }
                else
                {
                h1xdiff_n += (derxval-gradu(0))*(derxval-gradu(0))*weight*b_pos;
                h1ydiff_n += (deryval-gradu(1))*(deryval-gradu(1))*weight*b_pos;
                }
                }
                }*/ // if xfe

      ifjumpl2 = sqrt(ifjumpl2); iferr.Append(ifjumpl2);
      l2diff = l2diff_p + l2diff_n;
      h1xdiff = h1xdiff_p + h1xdiff_n;
      l2diff_p = sqrt(l2diff_p); l2err_p.Append(l2diff_p);
      h1xdiff_p = sqrt(h1xdiff_p); h1xerr_p.Append(h1xdiff_p);
      l2diff_n = sqrt(l2diff_n); l2err_n.Append(l2diff_n);
      h1xdiff_n = sqrt(h1xdiff_n); h1xerr_n.Append(h1xdiff_n);
      l2diff = sqrt(l2diff); l2err.Append(l2diff);
      h1xdiff = sqrt(h1xdiff); h1xerr.Append(h1xdiff);
      cout << " activeels = " << activeels << endl;
      cout << setw(12) << "l2_n" << "\t|";
      cout << setw(12) << "l2_p" << "\t|";
      cout << setw(12) << "l2" << "\t|";
      cout << setw(12) << "h1x_n" << "\t|";
      cout << setw(12) << "h1x_p" << "\t|";
      cout << setw(12) << "h1x" << "\t|";
      cout << setw(12) << "ifl2" << endl;

      for (int i = 0; i < refinements; ++i)
      {
        cout << setw(12) << l2err_n[i] << "\t|";
        cout << setw(12) << l2err_p[i] << "\t|";
        cout << setw(12) << l2err[i] << "\t|";
        cout << setw(12) << h1xerr_n[i] << "\t|";
        cout << setw(12) << h1xerr_p[i] << "\t|";
        cout << setw(12) << h1xerr[i] << "\t|";
        cout << setw(12) << iferr[i] << endl;
      }
      // cout << " l2diff_p = " << l2diff_p << endl;
      // cout << " l2diff_n = " << l2diff_n << endl;
      // cout << " l2diff = " << l2diff << endl;
      // cout << " h1diff_p = " << h1diff_p << endl;
      // cout << " h1diff_n = " << h1diff_n << endl;
      // cout << " h1diff = " << h1diff << endl;
      // cout << " ifjumpl2  = " << ifjumpl2  << endl;




    }
  };


/* ---------------------------------------- 
                  numproc
   ---------------------------------------- */
/*
    class NumProcMarkElementsOnInterface : public NumProc
    {
    protected:
        string gciname;
    public:
    
        NumProcMarkElementsOnInterface (PDE & apde, const Flags & flags)
            : NumProc (apde)
            { 
                gciname = flags.GetStringFlag("cutname","cut");
            }
  
        virtual string GetClassName () const
            {
                return "NPMarkElementsonInterface";
            }


        virtual void Do (LocalHeap & lh)
            {
                CutInfoContainer & cic = CutInfoContainer::getInstance();
                AdLinCutTriang & gci = *cic.Get(gciname);
                if (!gci.Finalized())
                  throw Exception("gci not ready yet!");
                for (int i = 0; i < pde.GetMeshAccess().GetNE(); ++i)
                {
                  if (gci.IsElementCut(i))
                    Ng_SetRefinementFlag (i+1, 1);
                  else
                    Ng_SetRefinementFlag (i+1, 0);
                }
            }
    };
*/




}

// static RegisterNumProc<NumProcMarkElementsOnInterface> npinitmark("markinterface");


static RegisterNumProc<NumProcSetValuesX> npinittestxfem2d("setvaluesx");
static RegisterNumProc<NumProcXDifference<2> > npxdiff("xdifference");
