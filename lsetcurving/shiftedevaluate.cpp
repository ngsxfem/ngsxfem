#define FILE_SHIFTEDEVALUATE_CPP
#include "shiftedevaluate.hpp"
#include <diffop_impl.hpp>
#include "../utils/ngsxstd.hpp"

namespace ngfem
{

  namespace
  {
    // Element-local data of a deformation grid function (back/forth), fetched
    // once per element (independent of the integration point).
    template <int SpaceD>
    struct DeformationData
    {
      bool valid = false;
      const ScalarFiniteElement<SpaceD> * scafe = nullptr;
      FlatMatrixFixWidth<SpaceD> vector {0, (double*)nullptr};
    };

    template <int SpaceD>
    DeformationData<SpaceD>
    GetDeformationData (const GridFunction * gf, ElementId elid, LocalHeap & lh)
    {
      DeformationData<SpaceD> data;
      if (!gf) return data;
      Array<int> dnums;
      gf->GetFESpace()->GetDofNrs(elid, dnums);
      FlatVector<> values(dnums.Size()*SpaceD, lh);
      gf->GetVector().GetIndirect(dnums, values);
      data.vector.Assign(FlatMatrixFixWidth<SpaceD> (dnums.Size(), &values(0)));
      FiniteElement & fe = gf->GetFESpace()->GetFE(elid, lh);
      data.scafe = &dynamic_cast<const ScalarFiniteElement<SpaceD> &> (fe);
      data.valid = true;
      return data;
    }

    // Solve Theta(Phi(x)) = z (with z = Phi_forth(ip)) for the reference point
    // ipx via the fixed point iteration. This is shared between the scalar and
    // the SIMD CalcMatrix/Apply/AddTrans so that both produce identical points.
    template <int SpaceD>
    IntegrationPoint
    SolveShiftedRefPoint (const MappedIntegrationPoint<SpaceD,SpaceD> & mip,
                          const DeformationData<SpaceD> & back,
                          const DeformationData<SpaceD> & forth,
                          LocalHeap & lh)
    {
      HeapReset hr(lh);
      IntegrationPoint ip(mip.IP());

      Vec<SpaceD> z = mip.GetPoint();

      if (forth.valid)
      {
        FlatVector<> shape_forth(forth.vector.Height(), lh);
        forth.scafe->CalcShape(ip, shape_forth);
        z += Trans(forth.vector) * shape_forth;
      }

      // Solve the problem Theta(Phi(x)) = z

      const double h = pow(abs(mip.GetJacobiDet()), 1./SpaceD);
      Vec<SpaceD> diff;
      IntegrationPoint ipx(ip);
      IntegrationPoint ipx0(0,0,0);
      MappedIntegrationPoint<SpaceD,SpaceD> mip_x0(ipx0, mip.GetTransformation());
      Vec<SpaceD> zdiff = z - mip_x0.GetPoint();
      int its = 0;

      if (back.valid)
      {
        FlatVector<> shape_back(back.vector.Height(), lh);
        Vec<SpaceD> dvec_back;

        Vec<SpaceD> ipx_best_so_far;
        double diff_best_so_far;
        bool first = true; int idx_best = 0;

        // Fixed point iteration
        while (its < globxvar.FIXED_POINT_ITER_TRESHOLD)
        {
          back.scafe->CalcShape(ipx, shape_back);
          dvec_back = Trans(back.vector) * shape_back;

          FlatVector<double> fv(SpaceD, &(ipx.Point())(0));

          diff = zdiff - dvec_back - mip.GetJacobian() * fv;
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
      }
      else
      {
        // Fixed point iteration
        while (its < globxvar.FIXED_POINT_ITER_TRESHOLD)
        {
          FlatVector<double> fv(SpaceD, &(ipx.Point())(0));
          diff = zdiff - mip.GetJacobian() * fv;
          if ( L2Norm(diff) < globxvar.EPS_SHIFTED_EVAL*h ) break;
          ipx.Point() = mip.GetJacobianInverse() * zdiff;
          its++;
        }
        if (its == globxvar.FIXED_POINT_ITER_TRESHOLD)
          throw Exception(" shifted eval took FIXED_POINT_ITER_TRESHOLD iterations and didn't (yet?) converge! ");
      }

      return ipx;
    }

    // SIMD version of the fixed point iteration: solves for the shifted
    // reference coordinates of a whole SIMD<IntegrationPoint> group (all lanes
    // at once), mathematically equivalent to the scalar SolveShiftedRefPoint.
    //
    // Like the scalar version, the geometry Jacobian (jac/jacinv) is FROZEN at
    // the original integration point; only the back-deformation displacement
    // dvec_back is recomputed at the current iterate. Lanes converge
    // independently via a frozen-mask (incl. best-so-far fallback / throw).
    template <int D>
    Vec<D,SIMD<double>>
    SolveShiftedRefPointSIMD (const SIMD<IntegrationPoint> & sip,      // original ip group
                              const Vec<D,SIMD<double>> & point0,      // Phi(ip0)
                              const Mat<D,D,SIMD<double>> & jac,
                              const Mat<D,D,SIMD<double>> & jacinv,
                              SIMD<double> jacdet,
                              const Vec<D,SIMD<double>> & phi0,        // Phi(0)
                              const DeformationData<D> & back,
                              const DeformationData<D> & forth,
                              SIMD<double> done_init,                  // 1.0 for padding lanes
                              LocalHeap & lh)
    {
      constexpr size_t W = SIMD<double>::Size();
      Vec<D,SIMD<double>> ip0ref = sip;

      // evaluate a deformation (back/forth) at reference coords -> displacement
      auto eval_dvec = [&] (const Vec<D,SIMD<double>> & ref,
                            const DeformationData<D> & dat) -> Vec<D,SIMD<double>>
      {
        HeapReset hr(lh);
        SIMD<IntegrationPoint> cur = sip;
        for (int d = 0; d < D; d++) cur(d) = ref(d);
        SIMD_IntegrationRule ir1(1, &cur);
        FlatMatrix<SIMD<double>> shape(dat.scafe->GetNDof(), 1, lh);
        dat.scafe->CalcShape(ir1, shape);
        Vec<D,SIMD<double>> dvec(SIMD<double>(0.0));
        for (size_t k = 0; k < dat.vector.Height(); k++)
          for (int d = 0; d < D; d++)
            dvec(d) += dat.vector(k,d) * shape(k,0);
        return dvec;
      };

      Vec<D,SIMD<double>> z = point0;
      if (forth.valid) z += eval_dvec(ip0ref, forth);
      Vec<D,SIMD<double>> zdiff = z - phi0;

      const SIMD<double> h = pow(fabs(jacdet), 1.0/D);
      const SIMD<double> eps_h = globxvar.EPS_SHIFTED_EVAL * h;
      const SIMD<double> eh2 = eps_h * eps_h;

      Vec<D,SIMD<double>> ipx = ip0ref;
      SIMD<double> done = done_init;          // 1.0 => lane frozen
      SIMD<double> best_nrm2 (1e99);
      Vec<D,SIMD<double>> best_ipx = ipx;
      int its = 0;

      while (its < globxvar.FIXED_POINT_ITER_TRESHOLD)
      {
        Vec<D,SIMD<double>> dvec_back (SIMD<double>(0.0));
        if (back.valid) dvec_back = eval_dvec(ipx, back);

        Vec<D,SIMD<double>> rhs = zdiff - dvec_back;
        Vec<D,SIMD<double>> diff = rhs - jac*ipx;
        SIMD<double> nrm2 = diff(0)*diff(0);
        for (int d = 1; d < D; d++) nrm2 += diff(d)*diff(d);

        // best-so-far (used as fallback in the back-case)
        SIMD<mask64> better = nrm2 < best_nrm2;
        best_nrm2 = If(better, nrm2, best_nrm2);
        for (int d = 0; d < D; d++) best_ipx(d) = If(better, ipx(d), best_ipx(d));

        done = If(nrm2 < eh2, SIMD<double>(1.0), done);
        if (HSum(done) >= double(W) - 0.5) break;   // all lanes done

        Vec<D,SIMD<double>> ipx_upd = jacinv * rhs;
        SIMD<mask64> frozen = SIMD<double>(0.5) < done;
        for (int d = 0; d < D; d++) ipx(d) = If(frozen, ipx(d), ipx_upd(d));
        its++;
      }

      SIMD<mask64> notdone = done < SIMD<double>(0.5);
      if (HSum(If(notdone, SIMD<double>(1.0), SIMD<double>(0.0))) > 0.5)
      {
        if (back.valid)
        {
          // bad lane: not converged and best iterate is no good candidate (L2 >= 1)
          SIMD<mask64> bad = notdone && (SIMD<double>(1.0) <= best_nrm2);
          if (HSum(If(bad, SIMD<double>(1.0), SIMD<double>(0.0))) > 0.5)
            throw Exception(" shifted eval took FIXED_POINT_ITER_TRESHOLD = "+to_string(globxvar.FIXED_POINT_ITER_TRESHOLD)+" iterations and didn't (yet?) converge! In addition, the best interation step is no good fallback candidate.");
          cout << IM(globxvar.NON_CONV_WARN_MSG_LVL) << "In Shifted_eval (SIMD): Not converged, using best iterate as fallback candidate." << endl;
          for (int d = 0; d < D; d++) ipx(d) = If(notdone, best_ipx(d), ipx(d));
        }
        else
          throw Exception(" shifted eval took FIXED_POINT_ITER_TRESHOLD iterations and didn't (yet?) converge! ");
      }

      return ipx;
    }
  }

  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              BareSliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {
    const MappedIntegrationPoint<SpaceD,SpaceD> & mip =
      static_cast<const MappedIntegrationPoint<SpaceD,SpaceD>&> (bmip);

    auto elid = mip.GetTransformation().GetElementId();
    auto data_back  = GetDeformationData<SpaceD> (back.get(), elid, lh);
    auto data_forth = GetDeformationData<SpaceD> (forth.get(), elid, lh);

    IntegrationPoint ipx = SolveShiftedRefPoint<SpaceD> (mip, data_back, data_forth, lh);

    MappedIntegrationPoint<SpaceD, SpaceD> mipx(ipx, mip.GetTransformation());
    evaluator->CalcMatrix(bfel, mipx, mat, lh);
  }

  template <int SpaceD>
  SIMD_BaseMappedIntegrationRule &
  DiffOpShiftedEval<SpaceD> ::
  CreateShiftedMIR (const SIMD_BaseMappedIntegrationRule & bmir,
                    LocalHeap & lh) const
  {
    const ElementTransformation & trafo = bmir.GetTransformation();
    auto elid = trafo.GetElementId();
    auto data_back  = GetDeformationData<SpaceD> (back.get(), elid, lh);
    auto data_forth = GetDeformationData<SpaceD> (forth.get(), elid, lh);

    const SIMD_IntegrationRule & ir = bmir.IR();
    const size_t n = ir.Size();
    const size_t nip = ir.GetNIP();
    constexpr size_t W = SIMD<IntegrationPoint>::Size();

    // Phi(0): physical image of the reference origin (same for all lanes).
    IntegrationPoint ip0(0,0,0);
    MappedIntegrationPoint<SpaceD,SpaceD> mip0(ip0, trafo);
    Vec<SpaceD,SIMD<double>> phi0;
    for (int d = 0; d < SpaceD; d++) phi0(d) = SIMD<double>(mip0.GetPoint()(d));

    const SIMD<double> lanes ([] (int l) { return double(l); });

    SIMD<IntegrationPoint> * mem = new (lh) SIMD<IntegrationPoint>[n];

    for (size_t i = 0; i < n; i++)
    {
      const SIMD<IntegrationPoint> & sip = ir[i];
      const auto & smip =
        static_cast<const SIMD<MappedIntegrationPoint<SpaceD,SpaceD>>&> (bmir[i]);

      // padding lanes (index >= nip) are marked done so they neither block the
      // iteration nor trigger spurious non-convergence.
      const size_t remain = (i*W < nip) ? min(W, nip - i*W) : size_t(0);
      SIMD<double> done_init = If(SIMD<double>(double(remain)) <= lanes,
                                  SIMD<double>(1.0), SIMD<double>(0.0));

      Vec<SpaceD,SIMD<double>> ipx =
        SolveShiftedRefPointSIMD<SpaceD> (sip, smip.GetPoint(), smip.GetJacobian(),
                                          smip.GetJacobianInverse(), smip.GetJacobiDet(),
                                          phi0, data_back, data_forth, done_init, lh);

      SIMD<IntegrationPoint> res = sip;   // keep weight / facetnr / unused coord
      for (int d = 0; d < SpaceD; d++) res(d) = ipx(d);
      mem[i] = res;
    }

    SIMD_IntegrationRule shifted_ir(n, mem);
    shifted_ir.SetNIP(nip);
    return trafo(shifted_ir, lh);
  }

  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  CalcMatrix (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> mat) const
  {
    static thread_local LocalHeap lh(10*1000*1000, "shiftedeval_simd_calcmatrix", true);
    HeapReset hr(lh);
    auto & smir = CreateShiftedMIR(bmir, lh);
    evaluator->CalcMatrix(bfel, smir, mat);
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
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x,
         BareSliceMatrix<SIMD<double>> flux) const
  {
    static thread_local LocalHeap lh(10*1000*1000, "shiftedeval_simd_apply", true);
    HeapReset hr(lh);
    auto & smir = CreateShiftedMIR(bmir, lh);
    evaluator->Apply(bfel, smir, x, flux);
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

  template <int SpaceD>
  void DiffOpShiftedEval<SpaceD> ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    static thread_local LocalHeap lh(10*1000*1000, "shiftedeval_simd_addtrans", true);
    HeapReset hr(lh);
    auto & smir = CreateShiftedMIR(bmir, lh);
    evaluator->AddTrans(bfel, smir, flux, x);
  }

  template class DiffOpShiftedEval<1>;
  template class DiffOpShiftedEval<2>;
  template class DiffOpShiftedEval<3>;

}


