#include "output.hpp"

namespace ngfem
{

  template<int D>
  void DoSpecialOutput (GridFunction * gfu, 
                        SolutionCoefficients<D> & solcoef, 
                        int subdivision, 
                        LocalHeap & lh)
  {
    ofstream outf("special.output");

    const double time = 0.0;
    // cout << " CalcXError at time = " << time << endl;
    Array<int> dnums;
    int activeels = 0;
    const MeshAccess & ma (gfu->GetFESpace().GetMeshAccess());
    for (int elnr = 0; elnr < ma.GetNE(); ++elnr)
    {
      HeapReset hr(lh);
      gfu -> GetFESpace().GetDofNrs (elnr, dnums);
      const int size = dnums.Size();

      FlatVector<double> elvec (size, lh);
      gfu -> GetVector().GetIndirect (dnums, elvec);   

      ElementTransformation & eltrans = ma.GetTrafo(elnr,false,lh);

      const FiniteElement & base_fel = gfu -> GetFESpace().GetFE(elnr,lh);
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
      FlatVector<> shapex(ndof_x,&shape_total(ndof));

      FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
      FlatMatrixFixWidth<D> dshape(ndof,&dshape_total(0,0));
      FlatMatrixFixWidth<D> dshapex(ndof_x,&dshape_total(ndof,0));
      if (xfe)
      {
          /*
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
        } // dt

        if (spacetime)
        {

        }
        else
        {

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
            
            const double jumpval = InnerProduct(jump,elvec);
            ifjumpl2 += weight * sqr(jumpval);
          }
        }
          */
      } // is xfe 
      else
      {
        DOMAIN_TYPE dt = dummfe->GetDomainType();

        Array<IntegrationPoint*> ips(4);
        IntegrationPoint ip1(0.0);
        ip1(0) = 1.0;
        IntegrationPoint ip2(0.0);
        ip2(1) = 1.0;
        IntegrationPoint ip3(0.0);
        ip3(2) = 1.0;
        
        ips[0] = &ip1;
        ips[1] = &ip2;
        ips[2] = &ip3;
        ips[3] = &ip1;
        
        // IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), intorder);
        for (int i = 0 ; i < ips.Size(); i++)
        {
          MappedIntegrationPoint<D,D> mip(*ips[i], eltrans);
          DimMappedIntegrationPoint<D+1> mipp(*ips[i], eltrans);
          mipp.Point().Range(0,D) = mip.GetPoint();
          mipp.Point()[D] = time;

          // double solval = dt == POS ? solcoef.GetSolutionPos().Evaluate(mipp) : solcoef.GetSolutionNeg().Evaluate(mipp);
          // Vec<D> soldval;
          // if (solcoef.HasSolutionDNeg() && solcoef.HasSolutionDPos())
          // {
          //   if (dt == POS)
          //     solcoef.GetSolutionDPos().Evaluate(mipp,soldval);
          //   else
          //     solcoef.GetSolutionDNeg().Evaluate(mipp,soldval);
          // }
          // else
          // {
          //   if (dt == POS)
          //     CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionPos()),mip,time,soldval,lh);
          //   else
          //     CalcDxShapeOfCoeff<D>(&(solcoef.GetSolutionNeg()),mip,time,soldval,lh);
          // }
          
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

          outf << mip.GetPoint() << "\t" << discval << endl;
        //   Vec<D> discdval = Trans(dshape) * elvec;
        //   Vec<D> diffdval = soldval - discdval;
        //   double diffdsqr = L2Norm2(diffdval);

        //   double fac = mip.GetWeight();
        //   if (dt == POS)
        //   {
        //     l2diff_p += b_pos*fac*sqr(discval-solval);
        //     h1diff_p += b_pos*fac*diffdsqr;
        //   }
        //   else
        //   {
        //     l2diff_n += b_neg*fac*sqr(discval-solval);
        //     h1diff_n += b_neg*fac*diffdsqr;
        //   }
        }
        outf << endl << endl;
      }
    }

  }

  template void DoSpecialOutput<2>(GridFunction * gfu, SolutionCoefficients<2> & solcoef, 
                                   int subdivision, LocalHeap & lh);
  template void DoSpecialOutput<3>(GridFunction * gfu, SolutionCoefficients<3> & solcoef, 
                                   int subdivision, LocalHeap & lh);


}
