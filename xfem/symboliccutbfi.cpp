/*********************************************************************/
/* File:   symboliccutbfi.cpp                                        */
/* Author: Christoph Lehrenfeld based on symbolicintegrator.cpp      */
/*         from Joachim Schoeberl (in NGSolve)                       */
/* Date:   September 2016                                            */
/*********************************************************************/
/*
   Symbolic cut integrators
*/

#include <fem.hpp>
#include "../xfem/symboliccutbfi.hpp"
#include "../cutint/xintegration.hpp"
namespace ngfem
{

  SymbolicCutBilinearFormIntegrator ::
  SymbolicCutBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                     shared_ptr<CoefficientFunction> acf,
                                     DOMAIN_TYPE adt,
                                     int aforce_intorder,
                                     int asubdivlvl,
                                     SWAP_DIMENSIONS_POLICY apol,
                                     VorB vb)
    : SymbolicBilinearFormIntegrator(acf,vb,VOL),
    cf_lset(acf_lset),
    dt(adt),
    force_intorder(aforce_intorder),
    subdivlvl(asubdivlvl),
    pol(apol)
  {
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(cf_lset,subdivlvl);
  }


  void 
  SymbolicCutBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    elmat = 0.0;
    T_CalcElementMatrixAdd<double,double> (fel, trafo, elmat, lh);
  }

  void
  SymbolicCutBilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & trafo,
                        FlatMatrix<double> elmat,
                        LocalHeap & lh) const
  {
    T_CalcElementMatrixAdd<double> (fel, trafo, elmat, lh);
  }

  void
  SymbolicCutBilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & trafo,
                        FlatMatrix<Complex> elmat,
                        LocalHeap & lh) const
  {
    if (fel.ComplexShapes() || trafo.IsComplex())
      T_CalcElementMatrixAdd<Complex,Complex> (fel, trafo, elmat, lh);
    else
      T_CalcElementMatrixAdd<Complex,double> (fel, trafo, elmat, lh);
  }

  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicCutBilinearFormIntegrator ::
  T_CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo,
                          FlatMatrix<SCAL> elmat,
                          LocalHeap & lh) const

  {
    static Timer t("symbolicCutBFI - CalcElementMatrix", 2);
    HeapReset hr(lh);
    // static Timer tstart("symbolicCutBFI - CalcElementMatrix startup", 2);
    // static Timer tstart1("symbolicCutBFI - CalcElementMatrix startup 1", 2);
    // static Timer tmain("symbolicCutBFI - CalcElementMatrix main", 2);

    /*
    static Timer td("symbolicCutBFI - CalcElementMatrix dmats", 2);
    static Timer tb("symbolicCutBFI - CalcElementMatrix diffops", 2);
    static Timer tdb("symbolicCutBFI - CalcElementMatrix D * B", 2);
    static Timer tlapack("symbolicCutBFI - CalcElementMatrix lapack", 2);
    */
    // tstart.Start();
    if (element_vb != VOL)
      {
        switch (trafo.SpaceDim())
          {
          // case 1:
          //   T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
          //   return;
          // case 2:
          //   T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
          //   return;
          // case 3:
          //   T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
          //   return;
          default:
            throw Exception ("EB not yet implemented");
            // throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }

    RegionTimer reg(t);

    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = fel_trial.Order()+fel_test.Order();

    auto et = trafo.GetElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;

    if (! (et == ET_SEGM || et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX) )
      throw Exception("SymbolicCutBFI can only treat simplices or hyperrectangulars right now");

    if (force_intorder >= 0)
      intorder = force_intorder;


    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    // elmat = 0;

    const IntegrationRule * ir1 = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, intorder, time_order, lh, subdivlvl, pol);

    if (ir1 == nullptr)
      return;
    ///
    const IntegrationRule * ir = nullptr;
    if (false && time_order > -1) //simple tensor product rule (no moving cuts with this..) ...
    {
       static bool warned = false;
       if (!warned)
       {
         cout << "WARNING: This is a pretty simple tensor product rule in space-time.\n";
         cout << "         A mapped integration rule of this will not see the time,\n";
         cout << "         but the underlying integration rule will." << endl;
         warned = true;
       }
       auto ir1D = SelectIntegrationRule (ET_SEGM, time_order);
       ir = new (lh) IntegrationRule(ir1->Size()*ir1D.Size(),lh);
       for (int i = 0; i < ir1D.Size(); i ++)
         for (int j = 0; j < ir->Size(); j ++)
           (*ir)[i*ir1->Size()+j] = IntegrationPoint((*ir1)[j](0),(*ir1)[j](1),ir1D[i](0),(*ir1)[j].Weight()*ir1D[i].Weight());
       //cout << *ir<< endl;
    }
    else
      ir = ir1;

    BaseMappedIntegrationRule & mir = trafo(*ir, lh);

    bool symmetric_so_far = false; // not reasonable for the Add-version
    /// WHAT FOLLOWS IN THIS FUNCTION IS COPY+PASTE FROM NGSOLVE !!!

    int k1 = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        for (auto proxy2 : test_proxies)
          {
            bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
            bool is_nonzero = false;

            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                if (nonzeros(l1+l, k1+k))
                  {
                    if (k != l) is_diagonal = false;
                    is_nonzero = true;
                  }


            if (is_nonzero)
              {
                HeapReset hr(lh);
                bool samediffop = *(proxy1->Evaluator()) == *(proxy2->Evaluator());
                // td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);


                if (!is_diagonal)
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        if (nonzeros(l1+l, k1+k))
                          {
                            if (k != l) is_diagonal = false;
                            is_nonzero = true;
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;

                            cf -> Evaluate (mir, val);
                            proxyvalues(STAR,k,l) = val.Col(0);
                          }
                        else
                          proxyvalues(STAR,k,l) = 0.0;
                      }
                else
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = k;

                      if (!elementwise_constant)
                        {
                          cf -> Evaluate (mir, val);
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val.Col(0);
                        }
                      else
                        {
                          cf -> Evaluate (mir[0], val.Row(0));
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val(0,0);
                        }
                    }

                // td.Stop();

                if (!mir.IsComplex())
                  {
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetWeight();
                  }
                else
                  { // pml
                    throw Exception("not treated yet (interface-weights!)");
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *=
                          static_cast<const ScalMappedIntegrationPoint<SCAL>&> (mir[i]).GetJacobiDet()*(*ir)[i].Weight();
                  }
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);


                enum { BS = 16 };
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min(int(BS), int(mir.Size())-i);

                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);

                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        /*
                        for (int i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        */
                        MultMatDiagMat(bbmat1, diagd, bdbmat1);
                        // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                      }
                    else
                      {
                        for (int j = 0; j < bs; j++)
                          {
                            int ii = i+j;
                            IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                            // bdbmat1.Cols(r2) = bbmat1.Cols(r1) * proxyvalues(ii,STAR,STAR);
                            MultMatMat (bbmat1.Cols(r1), proxyvalues(ii,STAR,STAR), bdbmat1.Cols(r2));
                          }
                        // tdb.AddFlops (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                      }
                    //  tdb.Stop();

                    // tlapack.Start();
                    // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                    // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                    symmetric_so_far &= samediffop && is_diagonal;
                    if (symmetric_so_far)
                      AddABtSym (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    else
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);

                    // tlapack.Stop();
                    // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                  }

                if (symmetric_so_far)
                  for (int i = 0; i < part_elmat.Height(); i++)
                    for (int j = i+1; j < part_elmat.Width(); j++)
                      part_elmat(i,j) = part_elmat(j,i);
              }

            l1 += proxy2->Dimension();
          }
        k1 += proxy1->Dimension();
      }
  }



  SymbolicCutFacetBilinearFormIntegrator ::
  SymbolicCutFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf_lset,
                                          shared_ptr<CoefficientFunction> acf,
                                          DOMAIN_TYPE adt,
                                          int aforce_intorder,
                                          int asubdivlvl)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false),
      cf_lset(acf_lset), dt(adt),
      force_intorder(aforce_intorder), subdivlvl(asubdivlvl)
  {
    simd_evaluate=false;
  }

  void  SymbolicCutFacetBilinearFormIntegrator::CalcFacetMatrix (
    const FiniteElement & fel1, int LocalFacetNr1,
    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
    const FiniteElement & fel2, int LocalFacetNr2,
    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
  {
    static Timer t_all("SymbolicCutFacetBilinearFormIntegrator::CalcFacetMatrix", 2);
    RegionTimer reg(t_all);
    elmat = 0.0;
    
    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());
    
    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    Facet2ElementTrafo transform1(eltype1, ElVertices1); 

    if (etfacet != ET_SEGM)
      if (etfacet != ET_TRIG || dt != IF)
        throw Exception("cut facet bilinear form can only do ET_SEGM-facets or IF right now");

    IntegrationRule * ir_facet = nullptr;

    if (etfacet != ET_SEGM && dt == IF)
    {
      static Timer t("symbolicCutBFI - CoDim2-hack", 2);
      RegionTimer reg(t);
      static bool first = true;
      if (first)
      {
        cout << "WARNING: unfitted codim-2 integrals are experimental!" << endl;
        cout << "         (and not performance-tuned)" << endl;
      }
      first = false;

      // Determine vertex values of the level set function:
      double lset[3];
      int v = 0;

      const POINT3D * verts_pts = ElementTopology::GetVertices(etfacet);

      Vec<2> verts[] {{verts_pts[0][0],verts_pts[0][1]},
        {verts_pts[1][0],verts_pts[1][1]},
        {verts_pts[2][0],verts_pts[2][1]}};
      bool haspos = false;
      bool hasneg = false;
      for (int i = 0; i < 3; i++)
      {
        IntegrationPoint ip = *(new (lh) IntegrationPoint(verts_pts[i][0],verts_pts[i][1]));
        
        const IntegrationPoint & ip_in_tet = transform1( LocalFacetNr1, ip);
        MappedIntegrationPoint<3,3> & mip = *(new (lh) MappedIntegrationPoint<3,3>(ip_in_tet,trafo1));
        lset[v] = cf_lset->Evaluate(mip);
        haspos = lset[v] > 0 ? true : haspos;
        hasneg = lset[v] < 0 ? true : hasneg;
        v++;
      }

      if (!hasneg || !haspos) return;

      // Determine the two cut positions:
      IntegrationRule cut(2,lh);
      IntegrationRule tetcut(2,lh);
      int ncut = 0;
      int edges[3][2] = {{0,1},{1,2},{2,0}};
      for (int edge = 0; edge < 3; edge++)
      {
        int node1 = edges[edge][0];
        int node2 = edges[edge][1];
        if (((lset[node1] > 0) && (lset[node2] < 0)) || ((lset[node1] < 0) && (lset[node2] > 0)))
        {
          Vec<2> cutpoint = 0.0;
          cutpoint += lset[node2] / (lset[node2]-lset[node1]) * verts[node1];
          cutpoint += lset[node1] / (lset[node1]-lset[node2]) * verts[node2];
          
          cut[ncut] = IntegrationPoint(cutpoint);
          tetcut[ncut] = transform1( LocalFacetNr1, cut[ncut]);
          ncut++;
        }
      }

      // Make an integration along the two cut points:

      // direction of the edge to integrate along
      auto tetdiffvec = tetcut[1].Point() - tetcut[0].Point();
      // rule for reference segm (not a  copy)
      IntegrationRule ir_tmp(ET_SEGM, 2*maxorder);
      const int npoints = ir_tmp.Size();
      // new rule
      ir_facet = new (lh) IntegrationRule(npoints,lh);

      for (int i = 0; i < npoints; i++)
      {
        IntegrationPoint & ip = (*ir_facet)[i];
        const double s = ir_tmp[i].Point()[0];
        ip.SetNr(ir_tmp[i].Nr());
        ip.Point() = (1-s) * cut[0].Point() + s * cut[1].Point();
      }

      MappedIntegrationRule<3,3> mir(*ir_facet,trafo1,lh);
      for (int i = 0; i < npoints; i++)
      {
        IntegrationPoint & ip = (*ir_facet)[i];
        auto F = mir[i].GetJacobian();
        auto mapped_diffvec = F * tetdiffvec;
        const double meas1D = L2Norm(mapped_diffvec);
        ip.SetWeight(ir_tmp[i].Weight() * meas1D);
      }
    }
    else //ET_SEGM
    {
      IntegrationPoint ipl(0,0,0,0);
      IntegrationPoint ipr(1,0,0,0);
      const IntegrationPoint & facet_ip_l = transform1( LocalFacetNr1, ipl);
      const IntegrationPoint & facet_ip_r = transform1( LocalFacetNr1, ipr);
      MappedIntegrationPoint<2,2> mipl(facet_ip_l,trafo1);
      MappedIntegrationPoint<2,2> mipr(facet_ip_r,trafo1);
      double lset_l = cf_lset->Evaluate(mipl);
      double lset_r = cf_lset->Evaluate(mipr);

      
      if ((lset_l > 0 && lset_r > 0) && dt != POS) return;
      if ((lset_l < 0 && lset_r < 0) && dt != NEG) return;

      if (dt == IF)
      {
        ir_facet = new (lh) IntegrationRule(1,lh);
        double xhat = - lset_l / (lset_r - lset_l );
        (*ir_facet)[0] = IntegrationPoint(xhat, 0, 0, 1.0);
        
      }
      else if ((lset_l > 0) != (lset_r > 0))
      {
        IntegrationRule ir_tmp (etfacet, 2*maxorder);
        ir_facet = new (lh) IntegrationRule(ir_tmp.Size(),lh);
        ///....CutIntegrationRule(cf_lset, trafo, dt, intorder, subdivlvl, lh);
        double x0 = 0.0;
        double x1 = 1.0;
        double xhat = - lset_l / (lset_r - lset_l );
        if ( ((lset_l > 0) && dt == POS) || ((lset_l < 0) && dt == NEG))
          x1 = xhat;
        else
          x0 = xhat;
        double len = x1-x0;
        for (int i = 0; i < ir_tmp.Size(); i++)
          (*ir_facet)[i] = IntegrationPoint(x0 + ir_tmp[i].Point()[0] * len, 0, 0, len*ir_tmp[i].Weight());
      }
      else
      {
        ir_facet = new (lh) IntegrationRule(etfacet, 2*maxorder);
      }
    }

    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, (*ir_facet), lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);

    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, (*ir_facet), lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          /*
          FlatVector<> measure(mir1.Size(), lh);
          switch (trafo1.SpaceDim())
            {
	    case 1:
              {
                Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                    Mat<1> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            case 2:
              {
                Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
                for (int i = 0; i < mir1.Size(); i++)
                  {
                    auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                    Mat<2> inv_jac = mip.GetJacobianInverse();
                    double det = mip.GetMeasure();
                    Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                    double len = L2Norm (normal);    // that's the surface measure 
                    normal /= len;                   // normal vector on physical element
                    const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                    measure(i) = len;
                  }
                break;
              }
            default:
              cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
            }
          */

          mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }
          if (dt == IF)
            for (int i = 0; i < mir1.Size(); i++)
              // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
              proxyvalues(i,STAR,STAR) *= (*ir_facet)[i].Weight();
          else
            for (int i = 0; i < mir1.Size(); i++)
              // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
              proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * (*ir_facet)[i].Weight();

          IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }
  SymbolicFacetBilinearFormIntegrator2 ::
  SymbolicFacetBilinearFormIntegrator2 (shared_ptr<CoefficientFunction> acf,
                                        int aforce_intorder)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false),
      force_intorder(aforce_intorder)
  {
    simd_evaluate=false;
  }

  void SymbolicFacetBilinearFormIntegrator2 ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    Facet2ElementTrafo transform2(eltype2, ElVertices2);
    
    IntegrationRule & ir_facet_vol1_tmp = transform1(LocalFacetNr1, ir_facet, lh);
    IntegrationRule & ir_facet_vol2_tmp = transform2(LocalFacetNr2, ir_facet, lh);

    IntegrationRule * ir_facet_vol1 = nullptr;
    IntegrationRule * ir_facet_vol2 = nullptr;

    
    if (time_order >= 0)
    {
      FlatVector<> st_point(3,lh);
      const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_order);

      auto ir_spacetime1 = new (lh) IntegrationRule (ir_facet_vol1_tmp.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_facet_vol1_tmp.Size(); j++)
        {
          const int ij = i*ir_facet_vol1_tmp.Size()+j;
          (*ir_spacetime1)[ij].SetWeight( ir_time[i].Weight() * ir_facet_vol1_tmp[j].Weight() );
          st_point = ir_facet_vol1_tmp[j].Point();
          st_point(2) = ir_time[i](0);
          (*ir_spacetime1)[ij].Point() = st_point;
        }
      }
      auto ir_spacetime2 = new (lh) IntegrationRule (ir_facet_vol2_tmp.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_facet_vol2_tmp.Size(); j++)
        {
          const int ij = i*ir_facet_vol2_tmp.Size()+j;
          (*ir_spacetime2)[ij].SetWeight( ir_time[i].Weight() * ir_facet_vol2_tmp[j].Weight() );
          st_point = ir_facet_vol2_tmp[j].Point();
          st_point(2) = ir_time[i](0);
          (*ir_spacetime2)[ij].Point() = st_point;
        }
      }
      ir_facet_vol1 = ir_spacetime1;
      ir_facet_vol2 = ir_spacetime2;
    }
    else
    {
      ir_facet_vol1 = &ir_facet_vol1_tmp;
      ir_facet_vol2 = &ir_facet_vol2_tmp;
    }

    // cout << *ir_facet_vol1 << endl;
    // cout << *ir_facet_vol2 << endl;
    // getchar();
    
    BaseMappedIntegrationRule & mir1 = trafo1(*ir_facet_vol1, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(*ir_facet_vol2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

          mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              
              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }


  SymbolicFacetPatchBilinearFormIntegrator ::
  SymbolicFacetPatchBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf,
                                            int aforce_intorder)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false),
      force_intorder(aforce_intorder)
  {
    simd_evaluate=false;
  }

  void SymbolicFacetPatchBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;
    if (trafo1.SpaceDim () > 2)
      throw Exception ("Patch integrator only implemented for 2D right now");
    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetPatchBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();

    IntegrationRule ir_vol1(eltype1, 2*maxorder);
    IntegrationRule ir_vol2(eltype2, 2*maxorder);

    IntegrationRule ir_patch1 (ir_vol1.Size()+ir_vol2.Size(),lh);
    IntegrationRule ir_patch2 (ir_vol1.Size()+ir_vol2.Size(),lh);

    Vec<2> vec;
    Vec<2> diff;
    Vec<2> update;

    
    for (int l = 0; l < ir_patch1.Size(); l++)
    { /// TODO : D == 2 or D == 3
      if (l<ir_vol1.Size())
      {
        ir_patch1[l] = ir_vol1[l];
        MappedIntegrationPoint<2,2> mip(ir_vol1[l], trafo1);
        // const double h = D==2 ? sqrt(mip.GetJacobiDet()) : cbrt(mip.GetJacobiDet());
        const double h = sqrt(mip.GetJacobiDet());
        IntegrationPoint ip_x0;
        vec = mip.GetPoint();
        int its = 0;
        double w = 0;
        while (its==0 || (L2Norm(diff) > 1e-8*h && its < 20))
        {
          MappedIntegrationPoint<2,2> mip_x0(ip_x0,trafo2);
          diff = vec - mip_x0.GetPoint();
          update = mip_x0.GetJacobianInverse() * diff;
          for (int d = 0; d < 2; ++d)
            ip_x0(d) += update(d);
          its++;
          w = mip_x0.GetMeasure();
        }
        ir_patch2[l] = ip_x0;
        ir_patch2[l].SetWeight(mip.GetWeight()/w);
      }
      else
      {
        ir_patch2[l] = ir_vol2[l-ir_vol1.Size()];
        MappedIntegrationPoint<2,2> mip(ir_vol2[l-ir_vol1.Size()], trafo2);
        // const double h = D==2 ? sqrt(mip.GetJacobiDet()) : cbrt(mip.GetJacobiDet());
        const double h = sqrt(mip.GetJacobiDet());
        IntegrationPoint ip_x0;
        vec = mip.GetPoint();
        int its = 0;
        double w = 0;
        while (its==0 || (L2Norm(diff) > 1e-8*h && its < 20))
        {
          MappedIntegrationPoint<2,2> mip_x0(ip_x0,trafo1);
          diff = vec - mip_x0.GetPoint();
          update = mip_x0.GetJacobianInverse() * diff;
          for (int d = 0; d < 2; ++d)
            ip_x0(d) += update(d);
          its++;
          w = mip_x0.GetMeasure();
        }
        ir_patch1[l] = ip_x0;
        ir_patch1[l].SetWeight(mip.GetWeight()/w);
      }
    }
    
    IntegrationRule * ir1 = nullptr;
    IntegrationRule * ir2 = nullptr;

    
    if (time_order >= 0)
    {
      FlatVector<> st_point(3,lh);
      const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_order);

      auto ir_spacetime1 = new (lh) IntegrationRule (ir_patch1.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_patch1.Size(); j++)
        {
          const int ij = i*ir_patch1.Size()+j;
          (*ir_spacetime1)[ij].SetWeight( ir_time[i].Weight() * ir_patch1[j].Weight() );
          st_point = ir_patch1[j].Point();
          st_point(2) = ir_time[i](0);
          (*ir_spacetime1)[ij].Point() = st_point;
        }
      }
      auto ir_spacetime2 = new (lh) IntegrationRule (ir_patch2.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_patch2.Size(); j++)
        {
          const int ij = i*ir_patch2.Size()+j;
          (*ir_spacetime2)[ij].SetWeight( ir_time[i].Weight() * ir_patch2[j].Weight() );
          st_point = ir_patch2[j].Point();
          st_point(2) = ir_time[i](0);
          (*ir_spacetime2)[ij].Point() = st_point;
        }
      }
      ir1 = ir_spacetime1;
      ir2 = ir_spacetime2;
    }
    else
    {
      ir1 = &ir_patch1;
      ir2 = &ir_patch2;
    }

    BaseMappedIntegrationRule & mir1 = trafo1(*ir1, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(*ir2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

          // mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          // mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            // proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();
            proxyvalues(i,STAR,STAR) *= mir1[i].GetWeight();

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              
              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }

}
