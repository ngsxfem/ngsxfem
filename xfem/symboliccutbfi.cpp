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
#include <variant>
#include "../xfem/symboliccutbfi.hpp"
#include "../cutint/xintegration.hpp"
#include "../cutint/straightcutrule.hpp"
#include "../utils/ngsxstd.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"
#include "../cutint/spacetimecutrule.hpp"

namespace ngfem
{


  SymbolicCutBilinearFormIntegrator ::
  SymbolicCutBilinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                     shared_ptr<CoefficientFunction> acf,
                                     VorB avb,
                                     VorB aelement_vb)
    : SymbolicBilinearFormIntegrator(acf,avb,aelement_vb)
  {
    lsetintdom = make_shared<LevelsetIntegrationDomain>(lsetintdom_in);

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
                        bool & symmetric_so_far,
                        LocalHeap & lh) const
  {
    symmetric_so_far = false;
    T_CalcElementMatrixAdd<double,double,double> (fel, trafo, elmat, lh);
  }

  void
  SymbolicCutBilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & trafo,
                        FlatMatrix<Complex> elmat,
                        bool & symmetric_so_far,
                        LocalHeap & lh) const
  {
    symmetric_so_far = false;
    if (fel.ComplexShapes() || trafo.IsComplex())
      T_CalcElementMatrixAdd<Complex,Complex> (fel, trafo, elmat, lh);
    else
      T_CalcElementMatrixAdd<Complex,double> (fel, trafo, elmat, lh);
  }


  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void SymbolicCutBilinearFormIntegrator ::
  T_CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<SCAL_RES> elmat,
                          LocalHeap & lh) const
    
  {
    static Timer timer("SymbolicCutBFI::CalcElementMatrixAdd");
    RegionTimer reg (timer);

    if (element_vb != VOL)
      {
        //throw Exception ("EB not yet implemented");
        T_CalcElementMatrixEBAdd<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo, elmat, lh);
        return;
      }

    bool is_mixedfe = typeid(fel) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = is_mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = is_mixedfe ? mixedfe->FETest() : fel;
    // size_t first_std_eval = 0;

    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = fel_trial.Order()+fel_test.Order();

    auto et = trafo.GetElementType();
    if (et == ET_SEGM || et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;

    if (! (et == ET_SEGM || et == ET_TRIG || et == ET_TET || et == ET_QUAD || et == ET_HEX) )
      throw Exception("SymbolicCutBFI can only treat simplices or hyperrectangulars right now");

    LevelsetIntegrationDomain lsetintdom_local(*lsetintdom);    
    if (lsetintdom_local.GetIntegrationOrder() < 0) // integration order shall not be enforced by lsetintdom
      lsetintdom_local.SetIntegrationOrder(intorder);
     bool symmetric_so_far = false; //we don't check for symmetry in the formulatin so far (TODO)!
    
		ProxyUserData ud;
    const_cast<ElementTransformation &>(trafo).userdata = &ud;

    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    if (ir == nullptr)
      return;

    if (simd_evaluate && globxvar.SIMD_EVAL)  
			try
      {
        SIMD_IntegrationRule & simd_ir = *(new (lh) SIMD_IntegrationRule(*ir, lh));
        FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(wei_arr, lh);
        auto &simd_mir = trafo(simd_ir, lh);

				PrecomputeCacheCF(cache_cfs, simd_mir, lh);

        // bool symmetric_so_far = true;
        int k1 = 0;
        int k1nr = 0;

        for (auto proxy1 : trial_proxies)
        {
          int l1 = 0;
          int l1nr = 0;
          for (auto proxy2 : test_proxies)
          {
            size_t dim_proxy1 = proxy1->Dimension();
            size_t dim_proxy2 = proxy2->Dimension();

            size_t tt_pair = l1nr * trial_proxies.Size() + k1nr;
            // first_std_eval = k1nr*test_proxies.Size()+l1nr;  // in case of SIMDException
            bool is_nonzero = nonzeros_proxies(tt_pair);
            bool is_diagonal = diagonal_proxies(tt_pair);

            if (is_nonzero)
            {
              HeapReset hr(lh);
              bool samediffop = same_diffops(tt_pair) && !is_mixedfe;
              // td.Start();

              FlatMatrix<SIMD<SCAL>> proxyvalues(dim_proxy1 * dim_proxy2, simd_ir.Size(), lh);
              FlatMatrix<SIMD<SCAL>> diagproxyvalues(dim_proxy1, simd_ir.Size(), lh);
              FlatMatrix<SIMD<SCAL>> val(1, simd_ir.Size(), lh);
              {
                // RegionTimer regdmat(timer_SymbBFIdmat);
                if (!is_diagonal)
                  for (size_t k = 0, kk = 0; k < dim_proxy1; k++)
                    for (size_t l = 0; l < dim_proxy2; l++, kk++)
                    {
                      if (nonzeros(l1 + l, k1 + k))
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;

                        cf->Evaluate(simd_mir, proxyvalues.Rows(kk, kk + 1));
                      }
                      else
                        ;
                      // proxyvalues.Row(kk) = 0.0;
                    }
                else
                  for (size_t k = 0; k < dim_proxy1; k++)
                  {
                    ud.trialfunction = proxy1;
                    ud.trial_comp = k;
                    ud.testfunction = proxy2;
                    ud.test_comp = k;

                    cf->Evaluate(simd_mir, diagproxyvalues.Rows(k, k + 1));
                  }
                // td.Stop();
              }
              // NgProfiler::StartThreadTimer (timer_SymbBFIscale, TaskManager::GetThreadId());
              FlatVector<SIMD<double>> weights(simd_ir.Size(), lh);
/*              if (!is_diagonal)
                for (size_t i = 0; i < ir.Size(); i++)
                  // proxyvalues.Col(i) *= mir[i].GetWeight();
                  weights(i) = ns_wei_arr[i]*mir[i].GetMeasure();
              else
                for (size_t i = 0; i < ir.Size(); i++)
                  diagproxyvalues.Col(i) *= ns_wei_arr[i]*mir[i].GetMeasure();*/

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
              SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);

              FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width() * dim_proxy1, simd_ir.Size(), lh);
              FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width() * dim_proxy2, simd_ir.Size(), lh);
              FlatMatrix<SIMD<SCAL_SHAPES>> bbmat2 = samediffop ? bbmat1 : FlatMatrix<SIMD<SCAL_SHAPES>>(elmat.Height() * dim_proxy2, simd_ir.Size(), lh);

              FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), dim_proxy2 * simd_ir.Size(),
                                              bdbmat1.Data());
              FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat2(elmat.Height(), dim_proxy2 * simd_ir.Size(),
                                                    bbmat2.Data());

              {
                proxy1->Evaluator()->CalcMatrix(fel_trial, simd_mir, bbmat1);
                if (!samediffop)
                  proxy2->Evaluator()->CalcMatrix(fel_test, simd_mir, bbmat2);
              }

              if (is_diagonal)
              {

                for (size_t j = 0; j < dim_proxy1; j++)
                {
                  auto hbbmat1 = bbmat1.RowSlice(j, dim_proxy1).Rows(r1);
                  auto hbdbmat1 = bdbmat1.RowSlice(j, dim_proxy1).Rows(r1);

                  for (size_t k = 0; k < bdbmat1.Width(); k++)
                    hbdbmat1.Col(k).Range(0, r1.Size()) = diagproxyvalues(j, k)*simd_wei_arr[k]*simd_mir[k].GetMeasure() * hbbmat1.Col(k);
                }
              }
              else
              {
                // static Timer t("DB", NoTracing);
                // RegionTracer reg(TaskManager::GetThreadId(), t);

                bdbmat1 = 0.0;
                //hbdbmat1.Rows(r1) = 0.0;

                for (size_t j = 0; j < dim_proxy2; j++)
                  for (size_t k = 0; k < dim_proxy1; k++)
                    if (nonzeros(l1 + j, k1 + k))
                    {
                      auto proxyvalues_jk = proxyvalues.Row(k * dim_proxy2 + j);
                      auto bbmat1_k = bbmat1.RowSlice(k, dim_proxy1).Rows(r1);
                      auto bdbmat1_j = bdbmat1.RowSlice(j, dim_proxy2).Rows(r1);

                      for (size_t i = 0; i < simd_ir.Size(); i++)
                        bdbmat1_j.Col(i).Range(0, r1.Size()) += proxyvalues_jk(i) *simd_wei_arr[i]*simd_mir[i].GetMeasure() * bbmat1_k.Col(i);
                    }
              }
                  symmetric_so_far &= samediffop && is_diagonal;
                      /*
                      if (symmetric_so_far)
                        AddABtSym (AFlatMatrix<double>(hbbmat2.Rows(r2)),
                                   AFlatMatrix<double> (hbdbmat1.Rows(r1)), part_elmat);
                      else
                        AddABt (AFlatMatrix<double> (hbbmat2.Rows(r2)),
                                AFlatMatrix<double> (hbdbmat1.Rows(r1)), part_elmat);
                      */

                      {
                        // static Timer t("AddABt", NoTracing);
                        // RegionTracer reg(TaskManager::GetThreadId(), t);
                        
                        if (symmetric_so_far)
                        {
                          /*
                            RegionTimer regdmult(timer_SymbBFImultsym);
                            NgProfiler::AddThreadFlops(timer_SymbBFImultsym, TaskManager::GetThreadId(),
                            SIMD<double>::Size()*2*r2.Size()*(r1.Size()+1)*hbbmat2.Width() / 2);
                          */
                          AddABtSym (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                        }
                      else
                        {
                          /*
                          RegionTimer regdmult(timer_SymbBFImult);
                          NgProfiler::AddThreadFlops(timer_SymbBFImult, TaskManager::GetThreadId(),
                                                     SIMD<double>::Size()*2*r2.Size()*r1.Size()*hbbmat2.Width());
                          */
                          AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                        }
                      }
                

              // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
              // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));
            }

            l1 += proxy2->Dimension();
            l1nr++;
          }
          k1 += proxy1->Dimension();
          k1nr++;
        }

        // ir.NothingToDelete();
        return;
      }
      catch (ExceptionNOSIMD e)
      {
        cout << IM(6) << e.What() << endl
             << "switching to scalar evaluation" << endl;
        simd_evaluate = false;
        throw ExceptionNOSIMD("in TCalcElementMatrixAdd");
        // T_CalcElementMatrixAdd<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo, elmat, lh);
        // return;
      }
			
    ///
    BaseMappedIntegrationRule & mir = trafo(*ir, lh);
    PrecomputeCacheCF(cache_cfs, mir, lh);
    
    // tstart.Stop();
    int k1 = 0;
    int k1nr = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        int l1nr = 0;
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

            if (is_nonzero) //   && k1nr*test_proxies.Size()+l1nr >= first_std_eval)
              {
                HeapReset hr(lh);
                bool samediffop = (*(proxy1->Evaluator()) == *(proxy2->Evaluator())) && !is_mixedfe;
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
                        proxyvalues(i,STAR,STAR) *= mir[i].GetMeasure()*wei_arr[i];
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetMeasure()*wei_arr[i];
                  }
                else
                  { // pml
                    throw Exception("not treated yet (interface-weights!)");
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetMeasure()*wei_arr[i];
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *=
                          static_cast<const ScalMappedIntegrationPoint<SCAL>&> (mir[i]).GetJacobiDet()*wei_arr[i];
                  }
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                
                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(size_t(BS), mir.Size()-i);
                    
                    FlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    FlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    FlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : FlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);

                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);

                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        FlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        for (size_t i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        // MultMatDiagMat(bbmat1, diagd, bdbmat1);
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
                    // tdb.Stop();
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
            l1nr++;
          }
        k1 += proxy1->Dimension();
        k1nr++;
      }
  }

  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void SymbolicCutBilinearFormIntegrator ::
    T_CalcElementMatrixEBAdd (const FiniteElement & fel,
                                   const ElementTransformation & trafo,
                                   FlatMatrix<SCAL_RES> elmat,
                                   LocalHeap & lh) const

                                       {
      static int timer = NgProfiler::CreateTimer ("symbolicBFI - CalcElementMatrix EB");
      if (lsetintdom->IsMultiLevelsetDomain())
        throw Exception("cut element boundary integrals not implemented for multi level sets");
      /*
      static Timer tir("symbolicBFI - CalcElementMatrix EB - intrules", 2);
      static Timer td("symbolicBFI - CalcElementMatrix EB - dmats", 2);
      static Timer tdb("symbolicBFI - CalcElementMatrix EB - b*d", 2);
      static Timer tb("symbolicBFI - CalcElementMatrix EB - bmats", 2);
      static Timer tmult("symbolicBFI - CalcElementMatrix EB - mult", 2);
      */
      NgProfiler::RegionTimer reg (timer);

      //elmat = 0;

      const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
      const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
      const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

      auto eltype = trafo.GetElementType();

      Facet2ElementTrafo transform(eltype, element_vb);
      int nfacet = transform.GetNFacets();
      
      const int order_sum = fel_trial.Order()+fel_test.Order();

      if (simd_evaluate && globxvar.SIMD_EVAL)
        // if (false)  // throwing the no-simd exception after some terms already added is still a problem
        {
          try
            {
              for (int k = 0; k < nfacet; k++)
                {
                  HeapReset hr(lh);
                  ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);

                  const IntegrationRule * ir_facet_tmp = nullptr;

                  if(etfacet == ET_SEGM){
                    IntegrationPoint ipl(0,0,0,0);
                    IntegrationPoint ipr(1,0,0,0);
                    const IntegrationPoint & facet_ip_l = transform( k, ipl);
                    const IntegrationPoint & facet_ip_r = transform( k, ipr);
                    MappedIntegrationPoint<2,2> mipl(facet_ip_l,trafo);
                    MappedIntegrationPoint<2,2> mipr(facet_ip_r,trafo);
                    double lset_l = lsetintdom->GetLevelsetGF()->Evaluate(mipl); //TODO: Not sure why that is seemingly better than cf_lset....
                    double lset_r = lsetintdom->GetLevelsetGF()->Evaluate(mipr);
      
                    if ((lset_l > 0 && lset_r > 0) && lsetintdom->GetDomainType() != POS) continue;
                    if ((lset_l < 0 && lset_r < 0) && lsetintdom->GetDomainType() != NEG) continue;
      
                    ir_facet_tmp = StraightCutIntegrationRuleUntransformed(Vec<2>{lset_r, lset_l}, ET_SEGM, lsetintdom->GetDomainType(), order_sum, FIND_OPTIMAL, lh);
                  } else if((etfacet == ET_TRIG) || (etfacet == ET_QUAD)){
                    int nverts = ElementTopology::GetNVertices(etfacet);
                    // Determine vertex values of the level set function:
                    vector<double> lset(nverts);
                    const POINT3D * verts_pts = ElementTopology::GetVertices(etfacet);
      
                    vector<Vec<2>> verts;
                    for(int i=0; i<nverts; i++) verts.push_back(Vec<2>{verts_pts[i][0], verts_pts[i][1]});
                    bool haspos = false;
                    bool hasneg = false;
                    for (int i = 0; i < nverts; i++) {
                      IntegrationPoint ip = *(new (lh) IntegrationPoint(verts_pts[i][0],verts_pts[i][1]));
      
                      const IntegrationPoint & ip_in_tet = transform( k, ip);
                      MappedIntegrationPoint<3,3> & mip = *(new (lh) MappedIntegrationPoint<3,3>(ip_in_tet,trafo));
      
                      //cout << "mip : " << mip.GetPoint() << endl;
                      lset[i] = lsetintdom->GetLevelsetGF()->Evaluate(mip);
                      //cout << "lset[i] : " << lset[i] << endl;
                      haspos = lset[i] > 0 ? true : haspos;
                      hasneg = lset[i] < 0 ? true : hasneg;
                    }

                    if(lsetintdom->GetDomainType() != POS && !hasneg) continue;
                    if(lsetintdom->GetDomainType() != NEG && !haspos) continue;
                    FlatVector<double> lset_fv(nverts, lh);
                    for(int i=0; i<nverts; i++){
                        lset_fv[i] = lset[i];
                        if(abs(lset_fv[i]) < 1e-16) throw Exception("lset val 0 in SymbolicCutFacetBilinearFormIntegrator");
                    }

                    LevelsetWrapper lsw(lset, etfacet);
                    ir_facet_tmp = StraightCutIntegrationRuleUntransformed(lset_fv, etfacet, lsetintdom->GetDomainType(), order_sum, FIND_OPTIMAL, lh);
                    //cout << "ir_facet_tmp: " << *ir_facet_tmp << endl;
                    Vec<3> tetdiffvec2(0.);

                    IntegrationRule & ir_scr_intet2 = transform( k, (*ir_facet_tmp), lh);
                    MappedIntegrationRule<3,3> mir3(ir_scr_intet2,trafo,lh);
                    int npoints = ir_facet_tmp->Size();
                    for (int i = 0; i < npoints; i++) {
                      IntegrationPoint & ip = (*ir_facet_tmp)[i];
                      Vec<3> normal = lsw.GetNormal(ip.Point());
                      Vec<2> tang = {normal[1],-normal[0]};

                      tetdiffvec2 = transform.GetJacobian( k, lh) * tang;
                      auto F = mir3[i].GetJacobian();
                      auto mapped_tang = F * tetdiffvec2;
                      const double ratio_meas1D = L2Norm(mapped_tang);
                      ip.SetWeight((*ir_facet_tmp)[i].Weight() * ratio_meas1D);
                    }
                  } else {
                    throw Exception("Only ET_SEGM, ET_TRIG and ET_QUAD are implemented yet.");
                  }

                  // SIMD_IntegrationRule ir_facet( (*ir_facet_tmp).Size(), lh);
                  SIMD_IntegrationRule ir_facet(*ir_facet_tmp, lh);
                  // NON-SIMD conform:
                  //for(int i=0; i<(*ir_facet_tmp).Size(); i++) ir_facet[i] = (*ir_facet_tmp)[i];
                  auto & ir_facet_vol = transform(k, ir_facet, lh);

                  SIMD_BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);

                  ProxyUserData ud;
                  const_cast<ElementTransformation&>(trafo).userdata = &ud;

                  mir.ComputeNormalsAndMeasure(eltype, k);

                  for (int k1 : Range(trial_proxies))
                    for (int l1 : Range(test_proxies))
                      {
                        if (!nonzeros_proxies(l1, k1)) continue;

                        auto proxy1 = trial_proxies[k1];
                        auto proxy2 = test_proxies[l1];
                        size_t dim_proxy1 = proxy1->Dimension();
                        size_t dim_proxy2 = proxy2->Dimension();
                        HeapReset hr(lh);
                        FlatMatrix<SIMD<SCAL>> proxyvalues(dim_proxy1*dim_proxy2, ir_facet.Size(), lh);

                        // td.Start();
                        for (int ck = 0; ck < dim_proxy1; ck++)
                          for (int l = 0; l < dim_proxy2; l++)
                            {
                              ud.trialfunction = proxy1;
                              ud.trial_comp = ck;
                              ud.testfunction = proxy2;
                              ud.test_comp = l;

                              auto kk = l + ck*dim_proxy2;
                              cf->Evaluate (mir, proxyvalues.Rows(kk, kk+1));
                              if (lsetintdom->GetDomainType() != IF) {
                                 for (size_t i = 0; i < mir.Size(); i++)
                                   proxyvalues(kk, i) *= mir[i].GetWeight();
                              }
                              else
                                  for (size_t i = 0; i < mir.Size(); i++)
                                    proxyvalues(kk, i) *= ir_facet[i].Weight();
                            }


                        IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                        IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                        SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);

                        FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width()*dim_proxy1, mir.Size(), lh);
                        FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width()*dim_proxy2, mir.Size(), lh);
                        bool samediffop = false; // not yet available
                        FlatMatrix<SIMD<SCAL_SHAPES>> bbmat2 = samediffop ?
                          bbmat1 : FlatMatrix<SIMD<SCAL_SHAPES>>(elmat.Height()*dim_proxy2, mir.Size(), lh);

                        FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), dim_proxy2*mir.Size(),
                                                        &bdbmat1(0,0));
                        FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat2(elmat.Height(), dim_proxy2*mir.Size(),
                                                              &bbmat2(0,0));

                        {
                          // ThreadRegionTimer regbmat(timer_SymbBFIbmat, TaskManager::GetThreadId());
                          proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);
                          if (!samediffop)
                            proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);
                        }

                        bdbmat1 = 0.0;
                        for (auto i : r1)
                          for (size_t j = 0; j < dim_proxy2; j++)
                            for (size_t ck = 0; ck < dim_proxy1; ck++)
                              {
                                auto res = bdbmat1.Row(i*dim_proxy2+j);
                                auto a = bbmat1.Row(i*dim_proxy1+ck);
                                auto b = proxyvalues.Row(ck*dim_proxy2+j);
                                res += pw_mult(a,b);
                              }

                        AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                      }
                }
              return;
            }

          catch (ExceptionNOSIMD e)
            {
              cout << IM(4) << e.What() << endl
                   << "switching to non-SIMD evaluation (in T_CalcElementMatrixEBAdd)" << endl;
              simd_evaluate = false;
              throw ExceptionNOSIMD("disabled simd-evaluate in AddElementMatrixEB");
            }
        }

      for (int k = 0; k < nfacet; k++)
        {
          // tir.Start();
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
          //if(etfacet != ET_SEGM) throw Exception("Only ET_SEGM support yet!");
          const IntegrationRule * ir_facet_tmp = nullptr;

          if(etfacet == ET_SEGM){
              IntegrationPoint ipl(0,0,0,0);
              IntegrationPoint ipr(1,0,0,0);
              const IntegrationPoint & facet_ip_l = transform( k, ipl);
              const IntegrationPoint & facet_ip_r = transform( k, ipr);
              MappedIntegrationPoint<2,2> mipl(facet_ip_l,trafo);
              MappedIntegrationPoint<2,2> mipr(facet_ip_r,trafo);
              double lset_l = lsetintdom->GetLevelsetGF()->Evaluate(mipl); //TODO: Not sure why that is seemingly better than cf_lset....
              double lset_r = lsetintdom->GetLevelsetGF()->Evaluate(mipr);

              if ((lset_l > 0 && lset_r > 0) && lsetintdom->GetDomainType() != POS) continue;
              if ((lset_l < 0 && lset_r < 0) && lsetintdom->GetDomainType() != NEG) continue;

              ir_facet_tmp = StraightCutIntegrationRuleUntransformed(Vec<2>{lset_r, lset_l}, ET_SEGM, lsetintdom->GetDomainType(), order_sum, FIND_OPTIMAL, lh);
          }
          else if((etfacet == ET_TRIG) || (etfacet == ET_QUAD)){
              int nverts = ElementTopology::GetNVertices(etfacet);
              // Determine vertex values of the level set function:
              vector<double> lset(nverts);
              const POINT3D * verts_pts = ElementTopology::GetVertices(etfacet);

              vector<Vec<2>> verts;
              for(int i=0; i<nverts; i++) verts.push_back(Vec<2>{verts_pts[i][0], verts_pts[i][1]});
              bool haspos = false;
              bool hasneg = false;
              for (int i = 0; i < nverts; i++)
              {
                IntegrationPoint ip = *(new (lh) IntegrationPoint(verts_pts[i][0],verts_pts[i][1]));

                const IntegrationPoint & ip_in_tet = transform( k, ip);
                MappedIntegrationPoint<3,3> & mip = *(new (lh) MappedIntegrationPoint<3,3>(ip_in_tet,trafo));

                //cout << "mip : " << mip.GetPoint() << endl;
                lset[i] = lsetintdom->GetLevelsetGF()->Evaluate(mip);
                //cout << "lset[i] : " << lset[i] << endl;
                haspos = lset[i] > 0 ? true : haspos;
                hasneg = lset[i] < 0 ? true : hasneg;
              }

              //if (!hasneg || !haspos) continue;
              if(lsetintdom->GetDomainType() != POS && !hasneg) continue;
              if(lsetintdom->GetDomainType() != NEG && !haspos) continue;
              FlatVector<double> lset_fv(nverts, lh);
              for(int i=0; i<nverts; i++){
                  lset_fv[i] = lset[i];
                  if(abs(lset_fv[i]) < 1e-16) throw Exception("lset val 0 in SymbolicCutFacetBilinearFormIntegrator");
              }

              LevelsetWrapper lsw(lset, etfacet);
              ir_facet_tmp = StraightCutIntegrationRuleUntransformed(lset_fv, etfacet, lsetintdom->GetDomainType(), order_sum, FIND_OPTIMAL, lh);
              //cout << "ir_facet_tmp: " << *ir_facet_tmp << endl;
              Vec<3> tetdiffvec2(0.);

              IntegrationRule & ir_scr_intet2 = transform( k, (*ir_facet_tmp), lh);
              MappedIntegrationRule<3,3> mir3(ir_scr_intet2,trafo,lh);
              int npoints = ir_facet_tmp->Size();
              for (int i = 0; i < npoints; i++)
              {
                  IntegrationPoint & ip = (*ir_facet_tmp)[i];
                  Vec<3> normal = lsw.GetNormal(ip.Point());
                  Vec<2> tang = {normal[1],-normal[0]};

                  tetdiffvec2 = transform.GetJacobian( k, lh) * tang;
                  auto F = mir3[i].GetJacobian();
                  auto mapped_tang = F * tetdiffvec2;
                  const double ratio_meas1D = L2Norm(mapped_tang);
                  ip.SetWeight((*ir_facet_tmp)[i].Weight() * ratio_meas1D);
              }
          }

          IntegrationRule ir_facet( (*ir_facet_tmp).Size(), lh);
          for(int i=0; i<(*ir_facet_tmp).Size(); i++) ir_facet[i] = (*ir_facet_tmp)[i];

          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);

          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);

          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;

          mir.ComputeNormalsAndMeasure(eltype, k);

          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                if (!nonzeros_proxies(l1, k1)) continue;

                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];

                HeapReset hr(lh);
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);

                // td.Start();
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;

                      cf->Evaluate (mir, val);
                      if(lsetintdom->GetDomainType() != IF){
                        for (int i = 0; i < mir.Size(); i++)
                          val(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
                      }
                      else
                          for (int i = 0; i < mir.Size(); i++)
                            val(i) *= ir_facet[i].Weight();
                      proxyvalues(STAR,k,l) = val.Col(0);
                    }
                // td.Stop();
                /*
                for (int i = 0; i < mir.Size(); i++)
                  {
                    tb.Start();
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                    proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
                    proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
                    tb.Stop();
                    tmult.Start();
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);

                    dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
                    elmat.Rows(r2).Cols(r1) += Trans (bmat2.Cols(r2)) * dbmat1.Cols(r1);
                    tmult.Stop();
                  }
                */

                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);

                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(BS, mir.Size()-i);

                    //TODO: Changed AFlatMatrix into FlatMatrix here... Fine?
                    FlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    FlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    FlatMatrix<SCAL_SHAPES> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);
                    proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();
                    bdbmat1 = 0.0;
                    // tdb.Start();

                    auto part_bbmat1 = bbmat1.Rows(r1);
                    auto part_bdbmat1 = bdbmat1.Rows(r1);
                    auto part_bbmat2 = bbmat2.Rows(r2);

                    for (int j = 0; j < bs; j++)
                      {
                        IntRange rj1 = proxy1->Dimension() * IntRange(j,j+1);
                        IntRange rj2 = proxy2->Dimension() * IntRange(j,j+1);
                        MultMatMat (part_bbmat1.Cols(rj1), proxyvalues(i+j,STAR,STAR), part_bdbmat1.Cols(rj2));
                      }

                    // tdb.Stop();

                    // tmult.Start();
                    AddABt (part_bbmat2, part_bdbmat1, part_elmat);
                    // part_elmat += part_bbmat2 * Trans(part_bdbmat1);
                    // tmult.Stop();
                    // tmult.AddFlops (r2.Size() * r1.Size() * bbmat2.Width());
                  }
              }
        }
    }


  SymbolicCutFacetBilinearFormIntegrator ::
  SymbolicCutFacetBilinearFormIntegrator (LevelsetIntegrationDomain & lsetintdom_in,
                                          shared_ptr<CoefficientFunction> acf)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false)
  {
    lsetintdom = make_shared<LevelsetIntegrationDomain>(lsetintdom_in);
    simd_evaluate=false;
    time_order = lsetintdom_in.GetTimeIntegrationOrder();
  }

  template <typename SCAL, typename SCAL_SHAPES>
  void  SymbolicCutFacetBilinearFormIntegrator::T_CalcFacetMatrix (
    const FiniteElement & fel1, int LocalFacetNr1,
    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
    const FiniteElement & fel2, int LocalFacetNr2,
    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
    FlatMatrix<SCAL> elmat,
    LocalHeap & lh) const
  {
    static int timer = NgProfiler::CreateTimer ("SymbolicCutFacetBilinearFormIntegrator::CalcFacetMatrix");
    NgProfiler::RegionTimer reg(timer);
    elmat = SCAL(0.0);

    if (lsetintdom->IsMultiLevelsetDomain())
      throw Exception("cut element boundary integrals not implemented for multi level sets");

    
    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());
    
    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    Facet2ElementTrafo transform1(eltype1, ElVertices1); 

    //if (etfacet != ET_SEGM){
        //if (time_order > -1) throw Exception("Time order > -1 not allowed in 3D.");
        //if (lsetintdom->GetDomainType() != IF) throw Exception("cut facet bilinear form can only do volume ints on ET_SEGM");
        //if (etfacet != ET_TRIG && etfacet != ET_QUAD) throw Exception("cut facet bilinear form can do IF ints only on ET_SEGM, ET_TRIG and ET_QUAD");
    //}
    //if(etfacet == ET_POINT) throw Exception("ET_POINT not implemented/ tested in SymbolicCutFacetBilinearFormIntegrator");

    // IntegrationRule * ir_facet = nullptr;
    const IntegrationRule * ir_scr = nullptr;

    Array<double> wei_arr;
    auto gflset = lsetintdom->GetLevelsetGF();
    if(gflset == nullptr) throw Exception("No gf in SymbolicCutFacetBilinearFormIntegrator::T_CalcFacetMatrix :(");

    Array<DofId> dnums(0,lh);
    gflset->GetFESpace()->GetDofNrs(trafo1.GetElementId(),dnums);
    FlatVector<> elvec(dnums.Size(),lh);
    gflset->GetVector().GetIndirect(dnums,elvec);
    
    shared_ptr<SpaceTimeFESpace> st_FE = nullptr;
    shared_ptr<NodalTimeFE> time_FE = nullptr;
    if(time_order > -1){
      st_FE = dynamic_pointer_cast<SpaceTimeFESpace >(gflset->GetFESpace());
      if(st_FE == nullptr) 
        throw Exception("Unable to cast SpaceTimeFESpace in SymbolicCutFacetBilinearFormIntegrator::T_CalcFacetMatrix");
      time_FE = dynamic_pointer_cast< NodalTimeFE>(st_FE->GetTimeFE());
      if(time_FE == nullptr) 
        throw Exception("Unable to cast time finite element in SymbolicCutFacetBilinearFormIntegrator::T_CalcFacetMatrix");
    }

    int NV = ElementTopology::GetNVertices(etfacet); // number of vertices of et_facet

    if (Dim(etfacet) == 2) 
    {
      //if (lsetintdom->GetDomainType() != IF) throw Exception("No Neg/pos facet ints in 3D yet");
      // so far only : Codim 2 special case (3D -> 1D)
      //if (time_order > -1) throw Exception("This is not possible for space_time");
      
      static Timer t("symbolicCutBFI - CoDim2", NoTracing);
      RegionTimer reg (t);

      IVec<4> int_tuple = 
        SwitchET<ET_HEX,ET_TET,ET_PRISM,ET_PYRAMID> (eltype1, [&LocalFacetNr1, &ElVertices1] (auto et)
         { return ET_trait<et>::GetFaceSort(LocalFacetNr1,ElVertices1); });

      // for(int i=0; i<NV; i++)
      //     if(abs(lset_fv[i]) < globxvar.EPS_INTERPOLATE_TO_P1 ) throw Exception("lset val 0 in SymbolicCutFacetBilinearFormIntegrator");

      FlatVector<> lset_fv(NV,lh);
      if(time_order < 0) {
        for (int i = 0; i < NV; i++)
          lset_fv(i) = elvec[int_tuple[i]];
        ir_scr = StraightCutIntegrationRuleUntransformed(lset_fv, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
      }
      else {
        int M = time_FE->GetNDof(); //time nodes
        int L = dnums.Size()/M; // dofs (or verts) per vol element 
        FlatVector<> cf_lset_at_element(NV*M, lh);
        FlatMatrix<> cf_lset_at_element_as_mat(M,NV,&cf_lset_at_element(0));

        for (int i = 0; i < NV; i++)
          cf_lset_at_element_as_mat.Col(i) = elvec.Slice(int_tuple[i],L);

        tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(cf_lset_at_element, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
      }
      if(ir_scr == nullptr) return;
      //cout << "ir_scr: " << *ir_scr << endl;
      if (lsetintdom->GetDomainType() == IF) // correct measure for 3D -> 1D case ( codim1 - adjustment is done automatically)
      {
        if (time_order > -1)
          throw Exception("correcting the IF-measure in space-time case is still missing");
        vector<double> lset(NV);
        for(int i=0; i<NV; i++) 
          lset[i] = lset_fv(i);
        LevelsetWrapper lsw(lset, etfacet);
        Vec<3> tetdiffvec2(0.);

        IntegrationRule & ir_scr_intet2 = transform1( LocalFacetNr1, (*ir_scr), lh);
        MappedIntegrationRule<3,3> mir3(ir_scr_intet2,trafo1,lh);
        int npoints = ir_scr->Size();

        for (int i = 0; i < npoints; i++)
        {
            IntegrationPoint & ip = (*ir_scr)[i];
            Vec<3> normal = lsw.GetNormal(ip.Point());
            Vec<2> tang = {normal[1],-normal[0]};

            tetdiffvec2 = transform1.GetJacobian( LocalFacetNr1, lh) * tang;
            auto F = mir3[i].GetJacobian();
            auto mapped_tang = F * tetdiffvec2;
            const double ratio_meas1D = L2Norm(mapped_tang);
            ip.SetWeight((*ir_scr)[i].Weight() * ratio_meas1D);
        }
      }
    }
    else if (Dim(etfacet) == 1) 
    {

        IVec<2> int_tuple = 
          SwitchET<ET_TRIG,ET_QUAD> (eltype1, [&LocalFacetNr1, &ElVertices1] (auto et) { return ET_trait<et>::GetEdgeSort(LocalFacetNr1,ElVertices1); });

        if(time_order < 0) {
            Vec<2> lset_vals_edge = {elvec[int_tuple[0]],elvec[int_tuple[1]]}; 
            ir_scr = StraightCutIntegrationRuleUntransformed(lset_vals_edge, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
        }
        else {
            int M = time_FE->GetNDof(); //time nodes
            int L = dnums.Size()/M; // dofs (or verts) per vol element 
            FlatVector<> cf_lset_at_element(NV*M, lh);
            FlatMatrix<> cf_lset_at_element_as_mat(M,NV,&cf_lset_at_element(0));

            for (int i = 0; i < NV; i++)
              cf_lset_at_element_as_mat.Col(i) = elvec.Slice(int_tuple[i],L);

            tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(cf_lset_at_element, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
        }
        if (ir_scr == nullptr) return;
    }
    else if (Dim(etfacet) == 0) {
      FlatVector <> lset_on_facet(elvec.Size()/2,lh);
      lset_on_facet = elvec.Slice(LocalFacetNr1,2);
      //LocalFacetNr1
        if(time_order < 0) {
            ir_scr = StraightCutIntegrationRuleUntransformed(lset_on_facet, etfacet, lsetintdom->GetDomainType(), 2*maxorder, FIND_OPTIMAL, lh);
        }
        else {
          tie( ir_scr, wei_arr) = SpaceTimeCutIntegrationRuleUntransformed(lset_on_facet, etfacet, time_FE.get(), lsetintdom->GetDomainType(), time_order, 2*maxorder, FIND_OPTIMAL,lh);
        }
        if (ir_scr == nullptr) return;

        //throw Exception("no 1D cut facet integration provided yet");
    }
    else
      throw Exception("no 1D cut facet integration provided yet");

    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, (*ir_scr), lh);

    Facet2ElementTrafo transform2(eltype2, ElVertices2);
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, (*ir_scr), lh);
    MarkAsSpaceTimeIntegrationRule(ir_facet_vol1);
    MarkAsSpaceTimeIntegrationRule(ir_facet_vol2);

    //The function operator () of Facet2ElementTrafo has in relation to ET_POINT the property of just transferring the first point correct.
    //The following lines of code repair this.
    if(transform1.FacetType(LocalFacetNr1) == ET_POINT){
        Vec<3> right_pnt;
        right_pnt = ir_facet_vol1[0].Point();
        for(int i=1; i<ir_facet_vol1.Size(); i++) ir_facet_vol1[i].Point() = right_pnt;
    }
    if(transform2.FacetType(LocalFacetNr2) == ET_POINT){
        Vec<3> right_pnt;
        right_pnt = ir_facet_vol2[0].Point();
        for(int i=1; i<ir_facet_vol2.Size(); i++) ir_facet_vol2[i].Point() = right_pnt;
    }

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    if (simd_evaluate && globxvar.SIMD_EVAL) {
      try {
        SIMD_IntegrationRule simd_ir_facet_vol1(ir_facet_vol1);
        SIMD_IntegrationRule simd_ir_facet_vol2(ir_facet_vol2);
        SIMD_IntegrationRule simd_ir_scr(*ir_scr);
        FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(wei_arr, lh);

        SIMD_BaseMappedIntegrationRule & simd_mir1 = trafo1(simd_ir_facet_vol1, lh);
        SIMD_BaseMappedIntegrationRule & simd_mir2 = trafo1(simd_ir_facet_vol2, lh);

        simd_mir1.SetOtherMIR(&simd_mir2);
        simd_mir2.SetOtherMIR(&simd_mir1);

        IntRange unified_r1(0, 0);
        for (int l1 : Range(test_proxies)) {
          HeapReset hr(lh);
          auto proxy2 = test_proxies[l1];

          FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width()*proxy2->Dimension(), simd_ir_facet_vol1.Size(), lh);
          FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), proxy2->Dimension()*simd_ir_facet_vol1.Size(), &bdbmat1(0,0));
          
          bdbmat1 = 0.0;

          for (int k1 : Range(trial_proxies)) {
            HeapReset hr(lh);
            auto proxy1 = trial_proxies[k1];

            FlatMatrix<SIMD<SCAL>> val(simd_mir1.Size(), 1, lh);
            FlatMatrix<SIMD<SCAL>> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), simd_ir_facet_vol1.Size(), lh);

            for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++) {
              for (size_t j = 0; j < proxy2->Dimension(); j++, kk++) {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = j;

                auto row = proxyvalues.Rows(kk, kk+1);
                cf->Evaluate(simd_mir1, row);
                
                if (lsetintdom->GetDomainType() == IF) // either 2D->0D (no need for weight correction) or 3D->1D ( weights are already corrected)
                  for (int i = 0; i < simd_mir1.Size(); i++){
                    proxyvalues(kk, i) *= ( time_order < 0 ? (simd_ir_scr)[i].Weight() : simd_wei_arr[i]);
                } else { // codim 1
                  for (int i = 0; i < simd_mir1.Size(); i++){
                    proxyvalues(kk,i) *= simd_mir1[i].GetMeasure() * ( time_order < 0 ? (simd_ir_scr)[i].Weight() : simd_wei_arr[i]);
                  }
                }
              }
            }

            IntRange trial_range = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1.GetNDof());
            IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);

            FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width()*proxy1->Dimension(), simd_ir_facet_vol1.Size(), lh);

            if (proxy1->IsOther())
              proxy1->Evaluator()->CalcMatrix(fel2, simd_mir2, bbmat1);
            else
              proxy1->Evaluator()->CalcMatrix(fel1, simd_mir1, bbmat1);
          
            for (auto i : trial_range) {
              for (size_t j = 0; j < proxy2->Dimension(); j++) {
                for (size_t k = 0; k < proxy1->Dimension(); k++) {
                  bdbmat1.Row(i*proxy2->Dimension()+j) += 
                    pw_mult(bbmat1.Row(i*proxy1->Dimension()+k), 
                            proxyvalues.Row(k*proxy2->Dimension()+j));
                }
              }
            }

            if (unified_r1.Next() == 0)
              unified_r1 = r1;
            else
              unified_r1 = IntRange (min2(r1.First(), unified_r1.First()),
                                        max2(r1.Next(), unified_r1.Next()));

          }
 
          IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);

          FlatMatrix<SIMD<SCAL>> bbmat2(elmat.Height()*proxy2->Dimension(), simd_ir_facet_vol1.Size(), lh);
          FlatMatrix<SIMD<SCAL>> hbbmat2(elmat.Height(), proxy2->Dimension()*simd_ir_facet_vol1.Size(), &bbmat2(0,0));

          if (proxy2->IsOther())
            proxy2->Evaluator()->CalcMatrix(fel2, simd_mir2, bbmat2);
          else
            proxy2->Evaluator()->CalcMatrix(fel1, simd_mir1, bbmat2);
        
          AddABt(hbbmat2.Rows(r2), hbdbmat1.Rows(unified_r1), elmat.Rows(r2).Cols(unified_r1));
        }

        return;
      } catch (ExceptionNOSIMD e) {
        cout << e.What() << endl
             << "switching to non-SIMD evaluation (in CalcFacetMatrix)" << endl;
      }
    }

    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<SCAL> val(mir1.Size(), 1,lh);

          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3,SCAL> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
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
          if (lsetintdom->GetDomainType() == IF) // either 2D->0D (no need for weight correction) or 3D->1D ( weights are already corrected)
              for (int i = 0; i < mir1.Size(); i++){
                    //proxyvalues(i,STAR,STAR) *= measure(i) * (*ir_scr)[i].Weight();
                    //proxyvalues(i,STAR,STAR) *= (*ir_facet)[i].Weight();
                    proxyvalues(i,STAR,STAR) *= ( time_order < 0 ? (*ir_scr)[i].Weight() : wei_arr[i]);
                    // proxyvalues(i,STAR,STAR) *= (*ir_scr)[i].Weight(); //The right choice for 2D...
                    //cout << "The two options: " << (*ir_facet)[i].Weight() << "\n" << measure(i) * (*ir_scr)[i].Weight() << endl;
                    //cout << "(*ir_scr)[i].Weight(): " << (*ir_scr)[i].Weight() << endl;
              }

          else // codim 1
          {
             // throw Exception("Foo!");
             for (int i = 0; i < mir1.Size(); i++){
                 //proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * (*ir_scr)[i].Weight();
                 // proxyvalues(i,STAR,STAR) *= measure(i) * ir_scr[i].Weight();
                 proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ( time_order < 0 ? (*ir_scr)[i].Weight() : wei_arr[i]);
                //  if(time_order < 0) proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * (*ir_scr)[i].Weight();
                //  else proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * wei_arr[i];
             }
          }

          IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<SCAL,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<SCAL,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

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
  SymbolicFacetBilinearFormIntegrator2 (shared_ptr<CoefficientFunction> acf)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false)
  {
    simd_evaluate=false;
  }

  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicFacetBilinearFormIntegrator2 ::
  T_CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<SCAL> elmat,
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

    Array<double> ir_st1_wei_arr;

    if (time_order >= 0)
    {
      FlatVector<> st_point(3,lh);
      const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_order);
      ir_st1_wei_arr.SetSize(ir_facet_vol1_tmp.Size()*ir_time.Size());

      auto ir_spacetime1 = new (lh) IntegrationRule (ir_facet_vol1_tmp.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_facet_vol1_tmp.Size(); j++)
        {
          const int ij = i*ir_facet_vol1_tmp.Size()+j;
          ir_st1_wei_arr[ij] = ir_time[i].Weight() * ir_facet_vol1_tmp[j].Weight() ;
          //(*ir_spacetime1)[ij].SetWeight( ir_time[i].Weight() * ir_facet_vol1_tmp[j].Weight() );
          st_point = ir_facet_vol1_tmp[j].Point();
          (*ir_spacetime1)[ij].Point() = st_point;
          (*ir_spacetime1)[ij].SetWeight(ir_time[i](0));
          MarkAsSpaceTimeIntegrationPoint((*ir_spacetime1)[ij]);
        }
      }
      auto ir_spacetime2 = new (lh) IntegrationRule (ir_facet_vol2_tmp.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        for (int j = 0; j < ir_facet_vol2_tmp.Size(); j++)
        {
          const int ij = i*ir_facet_vol2_tmp.Size()+j;
          //(*ir_spacetime2)[ij].SetWeight( ir_time[i].Weight() * ir_facet_vol2_tmp[j].Weight() );
          st_point = ir_facet_vol2_tmp[j].Point();
          (*ir_spacetime2)[ij].Point() = st_point;
          (*ir_spacetime2)[ij].SetWeight(ir_time[i](0));
          MarkAsSpaceTimeIntegrationPoint((*ir_spacetime2)[ij]);
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
    
    BaseMappedIntegrationRule & mir1 = trafo1(*ir_facet_vol1, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(*ir_facet_vol2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<SCAL> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3,SCAL> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

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

          for (int i = 0; i < mir1.Size(); i++){
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
              if(time_order == 0) proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();
              else proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_st1_wei_arr[i];
          }

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*fel1.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<SCAL,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<SCAL_SHAPES,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

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
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)              
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1); // | Lapack;
#else 
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
#endif
            }
        }
  }


  SymbolicFacetPatchBilinearFormIntegrator ::
  SymbolicFacetPatchBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf)
    : SymbolicFacetBilinearFormIntegrator(acf,VOL,false)
  {
    simd_evaluate=false;
  }

  // maps an integration point from inside one element to an integration point of the neighbor element
  // (integration point will be outside), so that the mapped points have the same coordinate
  template<int D>
  void MapPatchIntegrationPoint(IntegrationPoint & from_ip, const ElementTransformation & from_trafo,
                                const ElementTransformation & to_trafo, IntegrationPoint & to_ip,
                                LocalHeap & lh, bool spacetime_mode = false, double from_ip_weight =0.)
  {
    // cout << " ------------------------------------------- " << endl;
    HeapReset hr(lh);

    FlatVector<double> vec(D,lh);
    FlatVector<double> diff(D,lh);
    FlatVector<double> update(D,lh);

    MappedIntegrationPoint<D,D> mip(from_ip, from_trafo);
    const double h = sqrt(mip.GetJacobiDet());

    IntegrationPoint * ip_x0 = new(lh) IntegrationPoint(0,0,0);
    IntegrationPoint * ip_x00 = new(lh) IntegrationPoint(0,0,0);
    vec = mip.GetPoint();
    double w00 = 0;
    double first_diffnorm = 0;

    {
      HeapReset hr(lh);
      auto ip_a0 = new (lh) IntegrationPoint(0,0,0);
      if(spacetime_mode)
      {
        ip_a0->SetWeight(from_ip.Weight());
        MarkAsSpaceTimeIntegrationPoint(*ip_a0);
      }
      auto mip_a0 = new (lh) MappedIntegrationPoint<D,D>(*ip_a0,to_trafo);
      FlatMatrix<double> A(D,D,lh);
      FlatMatrix<double> Ainv(D,D,lh);
      FlatVector<double> f(D,lh);
      f = vec - mip_a0->GetPoint();
      auto ip_ai = new (lh) IntegrationPoint(0.,0.,0.);
      for (int d = 0; d < D ;  d++)
      {
        FlatVector<double> xhat(D,lh);
        for (int di = 0; di < 3;  di++)
        {
          if (di == d)
            ip_ai->Point()[di] = 1;
          else
            ip_ai->Point()[di] = 0;
        }
        if(spacetime_mode)
        {
          ip_ai->SetWeight(from_ip.Weight());
          MarkAsSpaceTimeIntegrationPoint(*ip_ai);
        }
        auto mip_ai = new (lh) MappedIntegrationPoint<D,D>(*ip_ai,to_trafo);
        A.Col(d) = mip_ai->GetPoint() - mip_a0->GetPoint();
      }
      CalcInverse(A,Ainv);
      w00 = abs(Det(A));
      ip_x00->Point().Range(0,D) = Ainv * f;
      ip_x0->Point().Range(0,D) = ip_x00->Point().Range(0,D);
    }

    int its = 0;
    double w = 0;
    while (its==0 || (L2Norm(diff) > globxvar.EPS_FACET_PATCH_INTEGRATOR*h && its < globxvar.NEWTON_ITER_TRESHOLD))
    {
      if(spacetime_mode) { ip_x0->SetWeight(from_ip.Weight()); MarkAsSpaceTimeIntegrationPoint(*ip_x0); }
      MappedIntegrationPoint<D,D> mip_x0(*ip_x0,to_trafo);
      diff = vec - mip_x0.GetPoint();
      if (its==0)
        first_diffnorm = L2Norm(diff);
      update = mip_x0.GetJacobianInverse() * diff;
      ip_x0->Point().Range(0,D) += update;
      its++;
      w = mip_x0.GetMeasure();
    }

    if(its >= globxvar.NEWTON_ITER_TRESHOLD || L2Norm(diff) > globxvar.EPS_FACET_PATCH_INTEGRATOR*h){
      cout << IM(globxvar.NON_CONV_WARN_MSG_LVL ) << "MapPatchIntegrationPoint: Newton did not converge after "
           << its <<" iterations! (" << D <<"D)" << endl;
      cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "taking a low order guess" << endl;
      cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "diff = " << first_diffnorm << endl;
      //globxvar.Output();
      cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "eps_treshold: " << globxvar.EPS_FACET_PATCH_INTEGRATOR << endl;

      to_ip = *ip_x00;
      if(spacetime_mode) to_ip.SetWeight(mip.GetMeasure() * from_ip_weight /w00);
      else to_ip.SetWeight(mip.GetWeight()/w00);
    }
    else if(L2Norm(ip_x0->Point() - ip_x00->Point()) > globxvar.MAX_DIST_NEWTON) {
        cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "Distance warning triggered, dist = " << L2Norm(ip_x0->Point() - ip_x00->Point()) << " its = " << its << endl;
        cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "taking a low order guess" << endl;
        to_ip = *ip_x00;
        if(spacetime_mode) to_ip.SetWeight(mip.GetMeasure() * from_ip_weight /w00);
        else to_ip.SetWeight(mip.GetWeight()/w00);
    }
    else
    {
      to_ip = *ip_x0;
      if(spacetime_mode) to_ip.SetWeight(mip.GetMeasure() * from_ip_weight /w);
      else to_ip.SetWeight(mip.GetWeight()/w);
    }
  }


  template<typename SCAL, typename SCAL_SHAPES>
  void SymbolicFacetPatchBilinearFormIntegrator ::
  T_CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<SCAL> elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;
    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetPatchBFI: LocalFacetNr2==-1");

    bool is_mixedfe1 = typeid(fel1) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe1 = static_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel1_trial = is_mixedfe1 ? mixedfe1->FETrial() : fel1;
    const FiniteElement & fel1_test = is_mixedfe1 ? mixedfe1->FETest() : fel1;
    bool is_mixedfe2 = typeid(fel2) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe2 = static_cast<const MixedFiniteElement*> (&fel2);
    const FiniteElement & fel2_trial = is_mixedfe2 ? mixedfe2->FETrial() : fel2;
    const FiniteElement & fel2_test = is_mixedfe2 ? mixedfe2->FETest() : fel2;

    int D = trafo1.SpaceDim();
    int maxorder = max2(max2 (fel1_trial.Order(), fel1_test.Order()), max2 (fel2_trial.Order(), fel2_test.Order()));

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();

    IntegrationRule ir_vol1(eltype1, 2*maxorder);
    IntegrationRule ir_vol2(eltype2, 2*maxorder);

    // cout << " ir_vol1 = " << ir_vol1 << endl;
    // cout << " ir_vol2 = " << ir_vol2 << endl;
    
    IntegrationRule ir_patch1 (ir_vol1.Size()+ir_vol2.Size(),lh);
    IntegrationRule ir_patch2 (ir_vol1.Size()+ir_vol2.Size(),lh);
    //In the non-space time case, the result of the mapping to the other element does not depend on the time
    //Therefore it is sufficient to do it once here.
    if(time_order == -1){
        for (int l = 0; l < ir_patch1.Size(); l++) {
            if (l<ir_vol1.Size()) {
                ir_patch1[l] = ir_vol1[l];
                if (D==2) MapPatchIntegrationPoint<2>(ir_patch1[l], trafo1, trafo2 ,ir_patch2[l], lh);
                else if(D==1) MapPatchIntegrationPoint<1>(ir_patch1[l], trafo1, trafo2 ,ir_patch2[l], lh);
                else MapPatchIntegrationPoint<3>(ir_patch1[l], trafo1, trafo2 ,ir_patch2[l], lh);
            }
            else {
                ir_patch2[l] = ir_vol2[l - ir_vol1.Size()];
                if (D==2) MapPatchIntegrationPoint<2>(ir_patch2[l], trafo2, trafo1 ,ir_patch1[l], lh);
                else if (D==1) MapPatchIntegrationPoint<1>(ir_patch2[l], trafo2, trafo1 ,ir_patch1[l], lh);
                else MapPatchIntegrationPoint<3>(ir_patch2[l], trafo2, trafo1 ,ir_patch1[l], lh);
            }
            ir_patch1[l].SetNr(l);
            ir_patch2[l].SetNr(l);
        }
    }

    // cout << " ir_patch1 = " << ir_patch1 << endl;
    // cout << " ir_patch2 = " << ir_patch2 << endl;
    
    IntegrationRule * ir1 = nullptr;
    IntegrationRule * ir2 = nullptr;

    Array<double> ir_st1_wei_arr;

    //Here we create approximately (up to the possible change in the element transformation due to the deformation)
    //tensor product quadrature rules for the space-time case
    if (time_order >= 0)
    {
      FlatVector<> st_point(3,lh);
      const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, time_order);

      auto ir_spacetime1 = new (lh) IntegrationRule (ir_patch1.Size()*ir_time.Size(),lh);
      ir_st1_wei_arr.SetSize(ir_spacetime1->Size());
      for (int i = 0; i < ir_time.Size(); i++)
      {
        double tval = ir_time[i](0);
        for (int j = 0; j < ir_patch1.Size(); j++)
        {
          const int ij = i*ir_patch1.Size()+j;

          //Task now: Calculate ir_patch1[j] in the spacetime setting at tval
          if(j< ir_vol1.Size()) ir_patch1[j] = ir_vol1[j];
          else {
            IntegrationPoint tmp = ir_vol2[j - ir_vol1.Size()];
            // ir_patch2[j] = ir_vol2[j - ir_vol1.Size()];
            double physical_weight = tmp.Weight();
            tmp.SetWeight(tval);
            MarkAsSpaceTimeIntegrationPoint(tmp);
            if (D==2) MapPatchIntegrationPoint<2>(tmp, trafo2, trafo1 ,ir_patch1[j], lh, true, physical_weight);
            else if (D==1) MapPatchIntegrationPoint<1>(tmp, trafo2, trafo1 ,ir_patch1[j], lh, true, physical_weight);
            else MapPatchIntegrationPoint<3>(tmp, trafo2, trafo1 ,ir_patch1[j], lh, true, physical_weight);
          }

          ir_st1_wei_arr[ij] = ir_time[i].Weight() * ir_patch1[j].Weight();

          st_point = ir_patch1[j].Point();
          (*ir_spacetime1)[ij].SetFacetNr(-1, VOL);
          (*ir_spacetime1)[ij].Point() = st_point;
          (*ir_spacetime1)[ij].SetWeight( tval);
          (*ir_spacetime1)[ij].SetNr(ij);
          MarkAsSpaceTimeIntegrationPoint((*ir_spacetime1)[ij]);
        }
      }
      auto ir_spacetime2 = new (lh) IntegrationRule (ir_patch2.Size()*ir_time.Size(),lh);
      for (int i = 0; i < ir_time.Size(); i++)
      {
        double tval = ir_time[i](0);
        for (int j = 0; j < ir_patch2.Size(); j++)
        {
          const int ij = i*ir_patch2.Size()+j;

          //Task now: Calculate ir_patch2[j] in the spacetime setting at tval
          if(j< ir_vol1.Size()) {
            IntegrationPoint tmp = ir_vol1[j];
            // ir_patch1[j] = ir_vol1[j];
            double physical_weight = tmp.Weight();
            tmp.SetWeight(tval);
            MarkAsSpaceTimeIntegrationPoint(tmp);
            if (D==2) MapPatchIntegrationPoint<2>(tmp, trafo1, trafo2 ,ir_patch2[j], lh, true, physical_weight);
            else if (D==1) MapPatchIntegrationPoint<1>(tmp, trafo1, trafo2 ,ir_patch2[j], lh, true, physical_weight);
            else MapPatchIntegrationPoint<3>(tmp, trafo1, trafo2 ,ir_patch2[j], lh, true, physical_weight);
          }
          else ir_patch2[j] = ir_vol2[j - ir_vol1.Size()];

          st_point = ir_patch2[j].Point();
          (*ir_spacetime2)[ij].SetFacetNr(-1, VOL);
          (*ir_spacetime2)[ij].Point() = st_point;
          (*ir_spacetime2)[ij].SetWeight(tval);
          (*ir_spacetime2)[ij].SetNr(ij);
          MarkAsSpaceTimeIntegrationPoint((*ir_spacetime2)[ij]);
        }
      }
      ir1 = ir_spacetime1;
      ir2 = ir_spacetime2;
      // cout << " *ir_spacetime1 = " << *ir_spacetime1 << endl;
      // cout << " *ir_spacetime2 = " << *ir_spacetime2 << endl;
    }
    else if (has_tref)
    {
      ir_st1_wei_arr.SetSize(ir_patch1.Size());
      for (int j = 0; j < ir_patch1.Size(); j++)
      {
        ir_st1_wei_arr[j] = ir_patch1[j].Weight();
        ir_patch1[j].SetWeight(tref);
        MarkAsSpaceTimeIntegrationPoint(ir_patch1[j]);
      }
      for (int j = 0; j < ir_patch2.Size(); j++)
      {
        ir_patch2[j].SetWeight(tref);
        MarkAsSpaceTimeIntegrationPoint(ir_patch2[j]);
      }
      ir1 = &ir_patch1;
      ir2 = &ir_patch2;
    }
    else
    {
      ir1 = &ir_patch1;
      ir2 = &ir_patch2;
    }

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    if (simd_evaluate && globxvar.SIMD_EVAL) {
      try {
        SIMD_IntegrationRule simd_ir1(*ir1, lh);
        SIMD_IntegrationRule simd_ir2(*ir2, lh);
        FlatArray<SIMD<double>> simd_ir_st1_wei_arr = CreateSIMD_FlatArray(ir_st1_wei_arr, lh);
        SIMD_BaseMappedIntegrationRule &simd_mir1 = trafo1(simd_ir1, lh);
        SIMD_BaseMappedIntegrationRule &simd_mir2 = trafo2(simd_ir2, lh);

        IntRange unified_r1(0, 0);
        for (int l1 : Range(test_proxies)) {
          HeapReset hr(lh);
          auto proxy2 = test_proxies[l1];

          FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width()*proxy2->Dimension(), simd_ir1.Size(), lh);
          FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), proxy2->Dimension()*simd_ir1.Size(), &bdbmat1(0,0));
          
          bdbmat1 = 0.0;

          for (int k1 : Range(trial_proxies)) {
            HeapReset hr(lh);
            auto proxy1 = trial_proxies[k1];

            FlatMatrix<SIMD<SCAL>> val(simd_mir1.Size(), 1, lh);
            FlatMatrix<SIMD<SCAL>> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), simd_ir1.Size(), lh);

            for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++) {
              for (size_t j = 0; j < proxy2->Dimension(); j++, kk++) {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = j;

                cf->Evaluate(simd_mir1, val);
                auto row = proxyvalues.Row(kk);
                for (auto m : Range(simd_mir1.Size()))
                  row(m) = val(m);
              }
            }

            for (size_t i = 0; i < simd_mir1.Size(); i++) {
              proxyvalues.Col(i) *= simd_mir1[i].GetMeasure()*simd_ir_st1_wei_arr[i];
            }

            IntRange trial_range = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof());
            IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2_trial : fel1_trial);

            FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width()*proxy1->Dimension(), simd_ir1.Size(), lh);

            if (proxy1->IsOther())
              proxy1->Evaluator()->CalcMatrix(fel2_trial, simd_mir2, bbmat1);
            else
              proxy1->Evaluator()->CalcMatrix(fel1_trial, simd_mir1, bbmat1);
          
            for (auto i : trial_range) {
              for (size_t j = 0; j < proxy2->Dimension(); j++) {
                for (size_t k = 0; k < proxy1->Dimension(); k++) {
                  bdbmat1.Row(i*proxy2->Dimension()+j) += 
                    pw_mult(bbmat1.Row(i*proxy1->Dimension()+k), 
                            proxyvalues.Row(k*proxy2->Dimension()+j));
                }
              }
            }

            if (unified_r1.Next() == 0)
              unified_r1 = r1;
            else
              unified_r1 = IntRange (min2(r1.First(), unified_r1.First()),
                                        max2(r1.Next(), unified_r1.Next()));

          }
 
          IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2_test : fel1_test);

          FlatMatrix<SIMD<SCAL>> bbmat2(elmat.Height()*proxy2->Dimension(), simd_ir1.Size(), lh);
          FlatMatrix<SIMD<SCAL>> hbbmat2(elmat.Height(), proxy2->Dimension()*simd_ir1.Size(), &bbmat2(0,0));

          if (proxy2->IsOther())
            proxy2->Evaluator()->CalcMatrix(fel2_test, simd_mir2, bbmat2);
          else
            proxy2->Evaluator()->CalcMatrix(fel1_test, simd_mir1, bbmat2);
        
          AddABt(hbbmat2.Rows(r2), hbdbmat1.Rows(unified_r1), elmat.Rows(r2).Cols(unified_r1));
        }

        return;
      } catch (ExceptionNOSIMD e) {
        cout << IM(6) << e.What() << endl
             << "switching to non-SIMD evaluation (in T_CalcFacetMatrix)" << endl;
      }
    }
    
    // non-SIMD case
    BaseMappedIntegrationRule & mir1 = trafo1(*ir1, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(*ir2, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<SCAL> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3,SCAL> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

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

          for (int i = 0; i < mir1.Size(); i++){
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            // proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();
            if((time_order >=0)||(has_tref)) proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure()*ir_st1_wei_arr[i];
            else proxyvalues(i,STAR,STAR) *= mir1[i].GetWeight();
          }

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*fel1_test.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*fel1_test.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<SCAL,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<SCAL,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2_trial, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1_trial, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2_test, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1_test, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              
              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2_trial : fel1_trial);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2_test : fel1_test);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
    return ;
  }



  void SymbolicCutBilinearFormIntegrator ::
  ApplyElementMatrix (const FiniteElement & fel, 
                      const ElementTransformation & trafo, 
                      const FlatVector<double> elx, 
                      FlatVector<double> ely,
                      void * precomputed,
                      LocalHeap & lh) const
  {
    // THIS IS UGLY CODE DUPLICATION FROM SYMBOLICCUTBFI::T_CalcElementMatrixAdd
    // and SymbolicBilinearFormIntegrator::ApplyElementMatrix (no simd)

    if (element_vb != VOL)
        throw Exception ("Apply for EB not yet implemented");
    
    static bool warned = false;
    if (!warned)
    {
      cout<<IM(3)<<"WARNING: The implementation of ApplyElementMatrix for cut elements is experimental.\n"<<endl;
      warned = true;
    }
    
    HeapReset hr(lh);
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

    LevelsetIntegrationDomain lsetintdom_local(*lsetintdom);    
    if (lsetintdom_local.GetIntegrationOrder() < 0) // integration order shall not be enforced by lsetintdom
      lsetintdom_local.SetIntegrationOrder(intorder);
	
    ProxyUserData ud(trial_proxies.Size(),lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
		
    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    ely = 0;
    if (ir == nullptr)
      return;

		if (simd_evaluate && globxvar.SIMD_EVAL) {
      try {
        SIMD_IntegrationRule & simd_ir = *(new (lh) SIMD_IntegrationRule(*ir, lh));
        FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(wei_arr, lh);
        auto &simd_mir = trafo(simd_ir, lh);

        PrecomputeCacheCF(cache_cfs, simd_mir, lh);

  
			  for (ProxyFunction * proxy : trial_proxies)
			  	ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);

			  for (ProxyFunction * proxy : trial_proxies)
			  	proxy->Evaluator()->Apply(fel_trial, simd_mir, elx, ud.GetAMemory(proxy));
  
			  FlatVector<double> ely1(ely.Size(), lh);
        
        ely = 0;
        FlatMatrix<SIMD<double>> val(simd_mir.Size(), 1,lh);
			  for (auto proxy : test_proxies)
        {
			  	HeapReset hr(lh);
			  	FlatMatrix<SIMD<double>> proxyvalues(simd_mir.Size(), proxy->Dimension(), lh);
			  	for (int k = 0; k < proxy->Dimension(); k++)
			  	{
			  		ud.testfunction = proxy;
			  		ud.test_comp = k;
			  		cf -> Evaluate (simd_mir, val);
			  		proxyvalues.Col(k) = val.Col(0);
			  	}
  
			  	for (int i = 0; i < simd_mir.Size(); i++)
			  		proxyvalues.Row(i) *= simd_mir[i].GetMeasure()*simd_wei_arr[i];
  
			  	proxy->Evaluator()->AddTrans(fel_test, simd_mir, proxyvalues, ely);
			  }
			  return;
		  } catch (ExceptionNOSIMD e) {
        cout << IM(6) << e.What() << endl
             << "switching to non-SIMD evaluation (in ApplyElementMatrix)" << endl;
      }  
    }
		
    ///
    BaseMappedIntegrationRule & mir = trafo(*ir, lh);
    
    
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir->GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      proxy->Evaluator()->Apply(fel_trial, mir, elx, ud.GetMemory(proxy), lh);
    
    // FlatVector<> ely1(ely.Size(), lh);
    FlatVector<double> ely1(ely.Size(), lh);   // can we really skip the <>  ???

    FlatMatrix<double> val(mir.Size(), 1,lh);
    for (auto proxy : test_proxies)
    {
      HeapReset hr(lh);
      FlatMatrix<double> proxyvalues(mir.Size(), proxy->Dimension(), lh);
      for (int k = 0; k < proxy->Dimension(); k++)
      {
        ud.testfunction = proxy;
        ud.test_comp = k;
        cf -> Evaluate (mir, val);
        proxyvalues.Col(k) = val.Col(0);
      }
        
      for (int i = 0; i < mir.Size(); i++)
        proxyvalues.Row(i) *= mir[i].GetMeasure()*wei_arr[i];

      proxy->Evaluator()->ApplyTrans(fel_test, mir, proxyvalues, ely1, lh);
      ely += ely1;
    }

    return;
  }



  void SymbolicCutBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo, 
                               FlatVector<double> elveclin,
                               FlatMatrix<double> elmat,
                               LocalHeap & lh) const
  {
    if (element_vb != VOL)
        throw Exception ("Apply for EB not yet implemented");

    static bool warned = false;
    if (!warned)
    {
      cout<<IM(3)<<"WARNING: The implementation of CalcLinearizedElementMatrix for cut elements is experimental.\n"<<endl;
      warned = true;
    }

    static int timer = NgProfiler::CreateTimer ("symboliccutbfi - calclinearized");
    NgProfiler::RegionTimer reg(timer);
    
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


    LevelsetIntegrationDomain lsetintdom_local(*lsetintdom);    
    if (lsetintdom_local.GetIntegrationOrder() < 0) // integration order shall not be enforced by lsetintdom
      lsetintdom_local.SetIntegrationOrder(intorder);
    
    ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;

    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(lsetintdom_local, trafo, lh);
    elmat = 0;
    if (ir == nullptr)
      return;

    if (simd_evaluate && globxvar.SIMD_EVAL) {
      try {
        SIMD_IntegrationRule & simd_ir = *(new (lh) SIMD_IntegrationRule(*ir, lh));
        FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(wei_arr, lh);
        SIMD_BaseMappedIntegrationRule &simd_mir = trafo(simd_ir, lh);

        PrecomputeCacheCF(cache_cfs, simd_mir, lh);

        for (ProxyFunction * proxy : trial_proxies) {
          ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);
          proxy->Evaluator()->Apply(fel_trial, simd_mir, elveclin, ud.GetAMemory(proxy));
        }

        for (CoefficientFunction * cf : gridfunction_cfs)
            ud.AssignMemory (cf, simd_ir.GetNIP(), cf->Dimension(), lh);
        // essentially an NGSolve copy+past from here on (2023-08-30)
        FlatMatrix<AutoDiff<1,SIMD<double>>> val(1, simd_mir.Size(), lh);
        IntRange unified_r1(0, 0);
        for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);              
          auto proxy2 = test_proxies[l1];
          FlatMatrix<SIMD<double>> bdbmat1(elmat.Width()*proxy2->Dimension(), simd_ir.Size(), lh);
          FlatMatrix<SIMD<double>> hbdbmat1(elmat.Width(), proxy2->Dimension()*simd_ir.Size(),
                                            &bdbmat1(0,0));
          bdbmat1 = 0.0;
          
          for (int k1 : Range(trial_proxies))
            {
              HeapReset hr(lh);
              auto proxy1 = trial_proxies[k1];
              
              FlatMatrix<SIMD<double>> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), simd_ir.Size(), lh);
              
              for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++)
                for (size_t l = 0; l < proxy2->Dimension(); l++, kk++)
                  {
                    ud.trialfunction = proxy1;
                    ud.trial_comp = k;
                    ud.testfunction = proxy2;
                    ud.test_comp = l;
                    // cf -> EvaluateDeriv (mir, val, proxyvalues.Rows(kk,kk+1));
                    cf -> Evaluate (simd_mir, val);
                    auto row = proxyvalues.Row(kk);
                    for (auto j : Range(simd_mir.Size()))
                      row(j) = val(j).DValue(0);
                  }
              
              for (size_t i = 0; i < simd_mir.Size(); i++)
                proxyvalues.Col(i) *= simd_mir[i].GetMeasure() * simd_wei_arr[i];

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
              
              FlatMatrix<SIMD<double>> bbmat1(elmat.Width()*proxy1->Dimension(), simd_ir.Size(), lh);
              
              // bbmat1 = 0.0;
              proxy1->Evaluator()->CalcMatrix(fel_trial, simd_mir, bbmat1);
              for (auto i : r1)
                for (size_t j = 0; j < proxy2->Dimension(); j++)
                  for (size_t k = 0; k < proxy1->Dimension(); k++)
                    {
                      bdbmat1.Row(i*proxy2->Dimension()+j) +=
                        pw_mult (bbmat1.Row(i*proxy1->Dimension()+k),
                                  proxyvalues.Row(k*proxy2->Dimension()+j));
                    }

              if (unified_r1.Next() == 0)
                unified_r1 = r1;
              else
                unified_r1 = IntRange (min2(r1.First(), unified_r1.First()),
                                        max2(r1.Next(), unified_r1.Next()));
            }

          IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
          
          FlatMatrix<SIMD<double>> bbmat2(elmat.Height()*proxy2->Dimension(), simd_ir.Size(), lh);
          FlatMatrix<SIMD<double>> hbbmat2(elmat.Height(), proxy2->Dimension()*simd_ir.Size(),
                                          &bbmat2(0,0));
          // bbmat2 = 0.0;
          proxy2->Evaluator()->CalcMatrix(fel_test, simd_mir, bbmat2);
          
          AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(unified_r1), elmat.Rows(r2).Cols(unified_r1));
        }

        return;
      } catch (ExceptionNOSIMD e) {
        cout << IM(6) << e.What() << endl
             << "switching to non-SIMD evaluation (in CalcLinearizedElementMatrix)" << endl;
      }
    }

    // non-SIMD case
    BaseMappedIntegrationRule & mir = trafo(*ir, lh);
    
    // ud.elx = &elveclin;
    // ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
    {
      ud.AssignMemory (proxy, ir->GetNIP(), proxy->Dimension(), lh);
      proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetMemory(proxy), lh);
    }
    for (CoefficientFunction * cf : gridfunction_cfs)
        ud.AssignMemory (cf, ir->GetNIP(), cf->Dimension(), lh);

    FlatMatrix<AutoDiff<1>> dval(mir.Size(), 1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
      {
        auto proxy1 = trial_proxies[k1];
        auto proxy2 = test_proxies[l1];

        FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
        for (int k = 0; k < proxy1->Dimension(); k++)
          for (int l = 0; l < proxy2->Dimension(); l++)
            // if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k)) // does no work for non-linear 
            if (true)
            {
              ud.trialfunction = proxy1;
              ud.trial_comp = k;
              ud.testfunction = proxy2;
              ud.test_comp = l;
                  
              // cf -> EvaluateDeriv (mir, val, deriv);
              // proxyvalues(STAR,l,k) = deriv.Col(0);
              cf -> Evaluate (mir, dval);
              for (size_t i = 0; i < mir.Size(); i++)
                proxyvalues(i,l,k) = dval(i,0).DValue(0);
            }
            else
              proxyvalues(STAR,l,k) = 0;

        for (int i = 0; i < mir.Size(); i++)
          proxyvalues(i,STAR,STAR) *= mir[i].GetMeasure()*wei_arr[i]; //mir[i].GetWeight();

        FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
        FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

        constexpr size_t BS = 16;
        for (size_t i = 0; i < mir.Size(); i+=BS)
        {
          int rest = min2(size_t(BS), mir.Size()-i);
          HeapReset hr(lh);
          FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

          for (int j = 0; j < rest; j++)
          {
            int ii = i+j;
            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
            proxy1->Evaluator()->CalcMatrix(fel_trial, mir[ii], bmat1, lh);
            proxy2->Evaluator()->CalcMatrix(fel_test, mir[ii], bmat2, lh);
            bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
            bbmat2.Rows(r2) = bmat2;
          }

          IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
          IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
          elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
        }
      }
      
    return;
  }
}
