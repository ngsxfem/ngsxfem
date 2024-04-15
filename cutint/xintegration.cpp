#include "xintegration.hpp"
#include "mlsetintegration.hpp"
#include "straightcutrule.hpp"
#include "spacetimecutrule.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"


namespace xintegration
{
  tuple<const IntegrationRule *, Array<double>> CreateCutIntegrationRule(const LevelsetIntegrationDomain & lsetintdom,
                                                                         const ElementTransformation & trafo,
                                                                         LocalHeap & lh)
  {
    static Timer timer("CreateCutIntegrationRule");
    RegionTimer reg (timer);
    if (lsetintdom.IsMultiLevelsetDomain())
    {
      if (lsetintdom.HasReferenceTime())
        throw Exception("No tref-fixing for mlset yet.");
      return CreateMultiLevelsetCutIntegrationRule(lsetintdom, trafo, lh);
    }
    else
    {
      auto gflset = lsetintdom.GetLevelsetGF();
      auto cflset = lsetintdom.GetLevelsetCF();
      int intorder = lsetintdom.GetIntegrationOrder();
      int time_intorder = lsetintdom.GetTimeIntegrationOrder();
      DOMAIN_TYPE dt = lsetintdom.GetDomainType();
      SWAP_DIMENSIONS_POLICY quad_dir_policy = lsetintdom.GetSwapDimensionPolicy();
      int subdivlvl = lsetintdom.GetNSubdivisionLevels();
      
      if (gflset != nullptr)
      {
        Array<DofId> dnums(0,lh);
        gflset->GetFESpace()->GetDofNrs(trafo.GetElementId(),dnums);
        FlatVector<> elvec(dnums.Size(),lh);
        gflset->GetVector().GetIndirect(dnums,elvec);
        if (time_intorder >= 0) {
          if (lsetintdom.HasReferenceTime())
            throw Exception("space-time quadrature rule shall not have a fixed reference time.");
          FESpace* raw_FE = (gflset->GetFESpace()).get();
          SpaceTimeFESpace * st_FE = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
          ScalarFiniteElement<1>* fe_time = nullptr;
          if (!st_FE)
          {
            static bool warned = false;
            if (!warned)
            { 
              warned = true;
              cout << IM(2) << "WARNING: You are using a space-time integration rule with an only spatial gridfunction (not a space time FE). This gridfunction will be treated as a constant-in-time function.\n" << endl;
            }
            fe_time = new (lh) L2HighOrderFE<ET_SEGM>(0);
          }
          else
            fe_time = dynamic_pointer_cast<ScalarFiniteElement<1>>(st_FE->GetTimeFE()).get();
          return SpaceTimeCutIntegrationRule(elvec, trafo, fe_time, dt, time_intorder, intorder, quad_dir_policy, lh);
        } else {
          const IntegrationRule * ir = StraightCutIntegrationRule(elvec, trafo, dt, intorder, quad_dir_policy, lh);
          if (ir != nullptr) 
          {
            if (lsetintdom.HasReferenceTime())
            {
              IntegrationRule * ir2 = new(lh) IntegrationRule(ir->Size(),lh);
              Array<double> wei_arr (ir->Size());
              for (int i = 0; i < ir->Size(); i ++)
              for(int i=0; i< ir->Size(); i++)
              { 
                wei_arr [i] = (*ir)[i].Weight();
                (*ir2)[i] = (*ir)[i];
                (*ir2)[i].SetWeight(lsetintdom.ReferenceTime());
                MarkAsSpaceTimeIntegrationPoint((*ir2)[i]);                
              }
              return make_tuple(ir2, wei_arr);
            }
            else 
            {
              Array<double> wei_arr (ir->Size());
              for(int i=0; i< ir->Size(); i++) 
                wei_arr [i] = (*ir)[i].Weight();
              return make_tuple(ir, wei_arr);
            }
          }
          else return make_tuple(nullptr, Array<double>());
        }
      }
      else if (cflset != nullptr)
      {
        if (lsetintdom.HasReferenceTime())
          throw Exception("No tref-fixing for old quadrature rules (yet).");
        if (time_intorder < 0) {
          const IntegrationRule * ir = CutIntegrationRule(cflset, trafo, dt, intorder, subdivlvl, lh);
          if(ir != nullptr) {
            Array<double> wei_arr (ir->Size());
            for(int i=0; i< ir->Size(); i++) wei_arr [i] = (*ir)[i].Weight();
            return make_tuple(ir, wei_arr);
          }
          else return make_tuple(nullptr, Array<double>());
        }
        else throw Exception("Space-time requires the levelset as a GridFunction!");
      }
      else throw Exception("Only null information provided, null integration rule served!");
      
    }
  }

  // old style (before introduction of LevelsetIntegrationDomain):
  tuple<const IntegrationRule *, Array<double>> CreateCutIntegrationRule(shared_ptr<CoefficientFunction> cflset,
                                                   shared_ptr<GridFunction> gflset,
                                                   const ElementTransformation & trafo,
                                                   DOMAIN_TYPE dt,
                                                   int intorder,
                                                   int time_intorder,
                                                   LocalHeap & lh,
                                                   int subdivlvl,
                                                   SWAP_DIMENSIONS_POLICY quad_dir_policy)
  {
    LevelsetIntegrationDomain lsetintdom(cflset, gflset,dt,intorder,time_intorder,subdivlvl,quad_dir_policy);
    return CreateCutIntegrationRule(lsetintdom, trafo, lh);
  }
  
  template<class SCAL> 
  FlatArray<SIMD<SCAL>> CreateSIMD_FlatArray(FlatArray<SCAL> ns_arr, LocalHeap & lh)
  {
    const int simd_blocks = (ns_arr.Size() + SIMD<IntegrationPoint>::Size() - 1) / SIMD<IntegrationPoint>::Size();
    FlatArray<SIMD<SCAL>> simd_arr(simd_blocks,lh);

    for (int i = 0; i < simd_blocks; i++){
      simd_arr[i] = [&] (int j)
      {
        const int nr = i * SIMD<IntegrationPoint>::Size() + j;
        if (nr < ns_arr.Size())
          return ns_arr[nr];
        return SCAL(0.0);
      };
    }
    return simd_arr;
  }


  template<int SD>
  PointContainer<SD>::PointContainer()
  {
#ifdef DEBUG
    k=0;
#endif
    pset.clear();
  };


  template<int SD>
  const Vec<SD>* PointContainer<SD>::operator()(const Vec<SD> & p)
  {

    // RegionTimer reg (timer);

    typename SetOfPoints::iterator it;
    it = pset.find(p);
    if (it == pset.end()) return &(*pset.insert(p).first);
    else
    {
#ifdef DEBUG
      k++;
#endif
      return &(*it);
    }
  }

  template<int SD>
  void PointContainer<SD>::Report(std::ostream & out) const
  {
    out << " PointContainer stored " << pset.size() << " points.\n";
#ifdef DEBUG
    out << " PointContainer rejected " << k << " points.\n";
#endif
  }

  template class PointContainer<2>;
  template class PointContainer<3>;
  template class PointContainer<4>;

  // ostream & operator<<(ostream & s, DOMAIN_TYPE dt)
  // {
  //   switch (dt)
  // {
  // case NEG:
  //     s << "NEG";
  //     break;
  // case POS:
  //     s << "POS";
  //     break;
  // case IF:
  //     s << "IF";
  //     break;
  // default:
  //     ;
  // }
  //   return s;
  // }


  template<>
  void FillSimplexWithRule<4> (const Simplex<4> & s, QuadratureRule<4> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 4");
  }

  template<>
  void FillSimplexWithRule<3> (const Simplex<3> & s, QuadratureRule<3> & quaddom, int intorder)
  {
    const double trafofac = Measure<3,3>(s.p) * 6.0;
    const IntegrationRule & ir = SelectIntegrationRule (ET_TET, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 3 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s.p[0]);
      for (int m = 0; m < 3 ;++m)
        point += ir[k](m) * (*s.p[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }

  template<>
  void FillSimplexWithRule<2> (const Simplex<2> & s, QuadratureRule<2> & quaddom, int intorder)
  {
    const double trafofac = Measure<2,2>(s.p) * 2.0;
    const IntegrationRule & ir = SelectIntegrationRule (ET_TRIG, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<2> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 2 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s.p[0]);
      for (int m = 0; m < 2 ;++m)
        point += ir[k](m) * (*s.p[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }


  template<>
  void FillSimplexWithRule<1> (const Simplex<1> & s, QuadratureRule<1> & quaddom, int intorder)
  {
    const double trafofac = Measure<1,1>(s.p);
    const IntegrationRule & ir = SelectIntegrationRule (ET_SEGM, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<1> point = ir[k](0) * (*s.p[1]) + (1.0 - ir[k](0)) * (*s.p[0]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }


  template<>
  void FillSimplexCoDim1WithRule<4> (const Array< const Vec<4> *> & s, const Vec<4> & pospoint,
                                     QuadratureRuleCoDim1<4> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 4");
  }

  template<>
  void FillSimplexCoDim1WithRule<3> (const Array< const Vec<3> *> & s, const Vec<3> & pospoint,
                                     QuadratureRuleCoDim1<3> & quaddom, int intorder)
  {
    Vec<3> a = *s[1] - *s[0];
    Vec<3> b = *s[2] - *s[0];
    Vec<3> c = Cross(a,b);
    const double trafofac = L2Norm(c);
    const double rel = max(L2Norm(a),L2Norm(b));

    if (trafofac < 1e-14 * rel)
      return;

    c /= trafofac;

    //sign check and probably switch
    Vec<3> d = pospoint - *s[0];
    if (InnerProduct(d,c) >= 0)
      c *= -1.0;


    const IntegrationRule & ir = SelectIntegrationRule (ET_TRIG, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 2 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s[0]);
      for (int m = 0; m < 2 ;++m)
        point += ir[k](m) * (*s[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
      quaddom.normals.Append(c);
    }
  }

  template<>
  void FillSimplexCoDim1WithRule<2> (const Array< const Vec<2> *> & s, const Vec<2> & pospoint,
                                     QuadratureRuleCoDim1<2> & quaddom, int intorder)
  {
    Vec<2> a = *s[1] - *s[0];
    Vec<2> n(-a(1),a(0));
    const double trafofac = L2Norm(a);

    const double rel = L2Norm(a);
    if (trafofac < 1e-14 * rel)
      return;

    n /= trafofac;

    //sign check and probably switch
    Vec<2> d = pospoint - *s[0];
    if (InnerProduct(d,n) >= 0)
      n *= -1.0;

    const IntegrationRule & ir = SelectIntegrationRule (ET_SEGM, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<2> point(0.0);
      point = (1.0-ir[k](0)) * (*s[0]);
      point += ir[k](0) * (*s[1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
      quaddom.normals.Append(n);
    }
  }



  double XLocalGeometryInformation::EvaluateLsetAtPoint( const IntegrationPoint & ip, double time) const
  {
    throw Exception("base class member function XLocalGeometryInformation::EvaluateLsetAtPoint called!");
    return 0.;
  }

  DOMAIN_TYPE XLocalGeometryInformation::MakeQuadRule() const
  {
    throw Exception("base class member function XLocalGeometryInformation::MakeQuadRule called!");
    return IF;
  }


  shared_ptr<XLocalGeometryInformation> XLocalGeometryInformation::Create(ELEMENT_TYPE ET_SPACE,
                                                                          ELEMENT_TYPE ET_TIME,
                                                                          const ScalarFieldEvaluator & a_lset,
                                                                          CompositeQuadratureRule<1> & a_compquadrule1,
                                                                          LocalHeap & a_lh,
                                                                          int a_int_order_space, int a_int_order_time,
                                                                          int a_ref_level_space, int a_ref_level_time)
  {
    return  XLocalGeometryInformation::Create(ET_SPACE,ET_TIME, a_lset,
                                              &a_compquadrule1, NULL, NULL, NULL,
                                              a_lh, a_int_order_space, a_int_order_time,
                                              a_ref_level_space, a_ref_level_time);
  }

  shared_ptr<XLocalGeometryInformation>XLocalGeometryInformation::Create(ELEMENT_TYPE ET_SPACE,
                                                                         ELEMENT_TYPE ET_TIME,
                                                                         const ScalarFieldEvaluator & a_lset,
                                                                         CompositeQuadratureRule<2> & a_compquadrule2,
                                                                         LocalHeap & a_lh,
                                                                         int a_int_order_space, int a_int_order_time,
                                                                         int a_ref_level_space, int a_ref_level_time)
  {
    return XLocalGeometryInformation::Create(ET_SPACE,ET_TIME, a_lset,
                                             NULL, &a_compquadrule2, NULL, NULL,
                                             a_lh, a_int_order_space, a_int_order_time,
                                             a_ref_level_space, a_ref_level_time);
  }

  shared_ptr<XLocalGeometryInformation> XLocalGeometryInformation::Create(ELEMENT_TYPE ET_SPACE,
                                                                          ELEMENT_TYPE ET_TIME,
                                                                          const ScalarFieldEvaluator & a_lset,
                                                                          CompositeQuadratureRule<3> & a_compquadrule3,
                                                                          LocalHeap & a_lh,
                                                                          int a_int_order_space, int a_int_order_time,
                                                                          int a_ref_level_space, int a_ref_level_time)
  {
    return XLocalGeometryInformation::Create(ET_SPACE,ET_TIME, a_lset,
                                             NULL, NULL, &a_compquadrule3, NULL,
                                             a_lh, a_int_order_space, a_int_order_time,
                                             a_ref_level_space, a_ref_level_time);
  }

  shared_ptr<XLocalGeometryInformation> XLocalGeometryInformation::Create(ELEMENT_TYPE ET_SPACE,
                                                                          ELEMENT_TYPE ET_TIME,
                                                                          const ScalarFieldEvaluator & a_lset,
                                                                          CompositeQuadratureRule<4> & a_compquadrule4,
                                                                          LocalHeap & a_lh,
                                                                          int a_int_order_space, int a_int_order_time,
                                                                          int a_ref_level_space, int a_ref_level_time)
  {
    return XLocalGeometryInformation::Create(ET_SPACE,ET_TIME, a_lset,
                                             NULL, NULL, NULL, &a_compquadrule4,
                                             a_lh, a_int_order_space, a_int_order_time,
                                             a_ref_level_space, a_ref_level_time);
  }

  shared_ptr<XLocalGeometryInformation> XLocalGeometryInformation::Create(ELEMENT_TYPE ET_SPACE,
                                                                ELEMENT_TYPE ET_TIME,
                                                                const ScalarFieldEvaluator & a_lset,
                                                                CompositeQuadratureRule<1> * a_compquadrule1,
                                                                CompositeQuadratureRule<2> * a_compquadrule2,
                                                                CompositeQuadratureRule<3> * a_compquadrule3,
                                                                CompositeQuadratureRule<4> * a_compquadrule4,
                                                                LocalHeap & a_lh,
                                                                int a_int_order_space, int a_int_order_time,
                                                                int a_ref_level_space, int a_ref_level_time)
  {
    if (ET_TIME == ET_POINT)
    {
      switch (ET_SPACE)
      {
      case ET_SEGM:
        return
          make_shared<NumericalIntegrationStrategy<ET_SEGM,ET_POINT>>
          (a_lset, *a_compquadrule1, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      case ET_TRIG:
        return
          make_shared<NumericalIntegrationStrategy<ET_TRIG,ET_POINT>>
          (a_lset, *a_compquadrule2, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      case ET_TET:
        return
          make_shared<NumericalIntegrationStrategy<ET_TET,ET_POINT>>
          (a_lset, *a_compquadrule3, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      default:
        throw Exception(" XLocalGeometryInformation * Create | ELEMENT_TYPE is not treated ");
        break;
      }
    }
    else // ET_SEGM
    {
      switch (ET_SPACE)
      {
      case ET_SEGM:
        return
          make_shared<NumericalIntegrationStrategy<ET_SEGM,ET_SEGM>>
          (a_lset, *a_compquadrule2, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      case ET_TRIG:
        return
          make_shared<NumericalIntegrationStrategy<ET_TRIG,ET_SEGM>>
          (a_lset, *a_compquadrule3, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      case ET_TET:
        return
          make_shared<NumericalIntegrationStrategy<ET_TET,ET_SEGM>>
          (a_lset, *a_compquadrule4, a_lh,
           a_int_order_space, a_int_order_time,
           a_ref_level_space, a_ref_level_time);
      default:
        throw Exception(" XLocalGeometryInformation * Create | ELEMENT_TYPE is not treated ");
        break;
      }
    }
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: NumericalIntegrationStrategy(const NumericalIntegrationStrategy & a, int reduce_ref_space, int reduce_ref_time)
    : XLocalGeometryInformation(a.lset), pc(a.pc),
      simplex_array_neg(a.simplex_array_neg),
      simplex_array_pos(a.simplex_array_pos),
      ref_level_space(a.ref_level_space-reduce_ref_space), ref_level_time(a.ref_level_time-reduce_ref_time),
      int_order_space(a.int_order_space), int_order_time(a.int_order_time),
      lh(a.lh), compquadrule(a.compquadrule)
  {
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: NumericalIntegrationStrategy(const ScalarFieldEvaluator & a_lset,
                                  PointContainer<SD> & a_pc,
                                  CompositeQuadratureRule<SD> & a_compquadrule,
                                  LocalHeap & a_lh,
                                  int a_int_order_space, int a_int_order_time,
                                  int a_ref_level_space, int a_ref_level_time)
    : XLocalGeometryInformation(&a_lset), pc(a_pc),
      ref_level_space(a_ref_level_space), ref_level_time(a_ref_level_time),
      int_order_space(a_int_order_space), int_order_time(a_int_order_time),
      lh(a_lh), compquadrule(a_compquadrule)
  {
    SetVerticesSpace();
    SetVerticesTime();
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: NumericalIntegrationStrategy(const ScalarFieldEvaluator & a_lset,
                                  CompositeQuadratureRule<SD> & a_compquadrule,
                                  LocalHeap & a_lh,
                                  int a_int_order_space, int a_int_order_time,
                                  int a_ref_level_space, int a_ref_level_time)
    : XLocalGeometryInformation(&a_lset), pc(*(new PointContainer<SD>())),
      ref_level_space(a_ref_level_space), ref_level_time(a_ref_level_time),
      int_order_space(a_int_order_space), int_order_time(a_int_order_time),
    lh(a_lh), compquadrule(a_compquadrule), ownpc(true)
  {
    SetVerticesSpace();
    SetVerticesTime();
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesSpace()
  {
    const POINT3D * verts = ElementTopology::GetVertices(ET_SPACE);
    const int nv = ElementTopology::GetNVertices(ET_SPACE);

    verts_space.SetSize(nv);
    for (int i = 0; i < nv; ++i)
      for (int d = 0; d < D; ++d)
        verts_space[i][d] = verts[i][d];
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesTime(const Array<double> & verts)
  {
    verts_time.SetSize(verts.Size());
    for (int i = 0; i < verts.Size(); ++i)
      verts_time[i] = verts[i];
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesTime()
  {
    const int np1dt = pow(2,ref_level_time);
    const double ht = 1.0 / np1dt;
    verts_time.SetSize(np1dt+1);
    for (int i = 0; i < np1dt + 1; ++i)
    {
      verts_time[i] = i * ht;
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesTimeFromUpperHalf(const Array< double >& verts_t)
  {
    const int newsize = (verts_t.Size()+1)/2;
    const int offset = (verts_t.Size()-1)/2;
    verts_time.SetSize(newsize);
    for (int i = 0; i < newsize; ++i)
    {
      verts_time[i] = verts_t[offset+i];
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesTimeFromLowerHalf(const Array< double >& verts_t)
  {
    const int newsize = (verts_t.Size()+1)/2;
    verts_time.SetSize(newsize);
    for (int i = 0; i < newsize; ++i)
    {
      verts_time[i] = verts_t[i];
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: SetVerticesSpace(const Array<Vec<D> > & verts)
  {
    verts_space.SetSize(verts.Size());
    for (int i = 0; i < verts.Size(); ++i)
      verts_space[i] = verts[i];
  }

  // Check prism for cut
  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: CheckIfCut() const
  {
    //enum { D = ET_trait<ET_SPACE>::DIM }; // spatial dimension
    //enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM}; // total dimension (space+time)

    //static Timer timer ("NumIntStrategy::CheckifCut (the prism check)");
    // ThreadRegionTimer reg (timer, TaskManager::GetThreadId());
    // RegionTimer reg (timer);

    bool haspos = false;
    bool hasneg = false;

    int np1ds = pow(2,ref_level_space);
    int np1dt = pow(2,ref_level_time);

    double dx_scalar = 1.0 / np1ds;

    Array<double> time (0);
    // double ht = 1.0;

    switch (ET_SPACE)
    {
    case ET_SEGM:
    case ET_TRIG:
    case ET_TET:
    {
      // int sum = 0;
      IVec< D > I;
      Vec< SD > position;
      for (int i = 0; i < D; ++i)
        I[i] = 0;

      // cout << " index = ";
      // for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
      //   cout << I[i] << ", \t";
      // cout << endl;

      bool finish = false;
      while (finish == false)
      {
        // do something

        //calculate all points corresponding to the current space position
        // loop over time points
        for (int i = 0; i < np1dt + 1; ++i)
        {
          for (int d = 0; d < D; ++d)
            position[d] = verts_space[0][d];

          for (int j = 0; j < D; ++j)
          {
            for (int d = 0; d < D; ++d)
            {
               position[d] += I[j] * dx_scalar * (verts_space[j+1][d] - verts_space[0][d]);
            }
          }
          // for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
          //   cout << position[d] << ",\t";

          if (ET_TIME == ET_SEGM)
          {
            position[ET_trait<ET_SPACE>::DIM] = verts_time[i];
            // cout << position[ET_trait<ET_SPACE>::DIM] << ",\t";
          }
          const ngfem::ScalarFieldEvaluator & eval (*lset);
          const double lsetval = eval(position);

          if (lsetval > distance_threshold)
            return POS;

          if (lsetval < -distance_threshold)
            return NEG;

          if (lsetval >= 0.0)
            haspos = true;
          else
            hasneg = true;
          // cout << " :: " << lsetval << ",\t";
          // cout << endl;

          if(haspos && hasneg)
          {
            // finish = true;
            // break;
            // cout << " IF " << endl;
            return IF;
          }
        }

        I[0]++;
        int sum = 0;
        for (int checkdim = 0; checkdim < D; ++checkdim)
        {
          sum = 0;
          for (int i = 0; i < D; ++i)
            sum += I[i];

          if (sum >= np1ds + 1)
          {
            if ( checkdim == D - 1)
            {
              finish = true;
              break;
            }
            else
            {
              I[checkdim] = 0;
              I[checkdim+1]++;
            }
          }
          else
            break;
        }
        // if (!finish)
        // {
        //   cout << " index = ";
        //   for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
        //     cout << I[i] << ", \t";
        //   cout << endl;
        // }
      }
      if (haspos)
        return POS;
      else
        return NEG;
    }
    break;
    default:
      throw Exception(" this CheckIfCut::ELEMENT_TYPE is not treated... yet ");
      break;
    }

  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE NumericalIntegrationStrategy<ET_SPACE,ET_TIME>
  :: MakeQuadRule() const
  {

    //enum { D = ET_trait<ET_SPACE>::DIM }; // spatial dimension
    //enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM}; // total dimension (space+time)

    // check with the help of regularly distributed points if current
    // space(-time) geometry is cut (has different sign in lset-value)
    DOMAIN_TYPE dt_self = CheckIfCut();

    if (dt_self == IF)
    {
      // cout << " cut " << endl;

      bool refine_time = ET_TIME == ET_SEGM && ref_level_time > 0;
      bool refine_space = ref_level_space > 0;

      // only refinement in space or time for the recursive call
      // decide which refinement should come first here:
      if (refine_time && refine_space)
      {
        if (false && ref_level_time >= ref_level_space) //better suited for high (relative) velocities...
        {
          refine_time = true;
          refine_space = false;
        }
        else //better suited for slow (relative) velocities...
        {
          refine_time = false;
          refine_space = true;
        }
      }

      // divide space-time prism into upper and lower half
      if (refine_time)
      {
        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_upper (*this, 0, 1);
        numint_upper.SetVerticesSpace(verts_space);
        numint_upper.SetVerticesTimeFromUpperHalf(verts_time);
        numint_upper.SetDistanceThreshold(distance_threshold);
        numint_upper.MakeQuadRule();  // recursive call!

        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_lower (*this, 0, 1);
        numint_lower.SetVerticesSpace(verts_space);
        numint_lower.SetVerticesTimeFromLowerHalf(verts_time);
        numint_lower.SetDistanceThreshold(distance_threshold);
        numint_lower.MakeQuadRule();  // recursive call!
      }

      // divide space-time prism into 4/8 prism of same height
      if (refine_space)
      {
        if ( ET_SPACE == ET_TRIG)
        {
          // barycentric coordinates for new points
          const double baryc[6][3] = { { 0.0, 0.0, 1.0},
                                        { 0.5, 0.0, 0.5},
                                        { 1.0, 0.0, 0.0},
                                        { 0.0, 0.5, 0.5},
                                        { 0.5, 0.5, 0.0},
                                        { 0.0, 1.0, 0.0}};

          // new triangles as connectivity information of the vertices baryc above
          const int trigs[4][3] = { { 0, 1, 3},
                                     { 1, 2, 4},
                                     { 1, 3, 4},
                                     { 3, 4, 5}};

          for (int i = 0; i < 4; ++i) // triangles
          {
            NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_i (*this, 1, 0);
            numint_i.SetVerticesTime(verts_time);
            Array< Vec<D> > newverts(3);
            for (int j = 0; j < 3; ++j) //vertices
            {
              newverts[j] = Vec<D>(0.0);
              for (int d = 0; d < 3; ++d)
                newverts[j] += baryc[trigs[i][j]][d] * verts_space[d];
            }
            numint_i.SetVerticesSpace(newverts);
            if (ET_TIME == ET_POINT)
              numint_i.SetDistanceThreshold(0.5*distance_threshold);
            else
              numint_i.SetDistanceThreshold(distance_threshold);
            numint_i.MakeQuadRule(); // recursive call!
          }

        }
        else if ( ET_SPACE == ET_SEGM)
        {
            // barycentric coordinates for new points
            const double baryc[3][2] = { { 0.0, 1.0},
                                         { 0.5, 0.5},
                                         { 1.0, 0.0}};

            // new segms as connectivity information of the vertices baryc above
            const int segm[2][2] = { { 0, 1},
                                     { 1, 2}};

            for (int i = 0; i < 2; ++i) // segms
            {
              NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_i (*this, 1, 0);
              numint_i.SetVerticesTime(verts_time);
              Array< Vec<D> > newverts(2);
              for (int j = 0; j < 2; ++j) //vertices
              {
                newverts[j] = Vec<D>(0.0);
                for (int d = 0; d < 2; ++d)
                  newverts[j] += baryc[segm[i][j]][d] * verts_space[d];
              }
              numint_i.SetVerticesSpace(newverts);
              numint_i.SetDistanceThreshold(0.5*distance_threshold);
              numint_i.MakeQuadRule(); // recursive call!
            }
        }
        else if ( ET_SPACE == ET_TET)
        {
          const double baryc[10][4] = { { 0.0, 0.0, 0.0, 1.0},
                                        { 0.5, 0.0, 0.0, 0.5},
                                        { 1.0, 0.0, 0.0, 0.0},
                                        { 0.0, 0.5, 0.0, 0.5},
                                        { 0.5, 0.5, 0.0, 0.0},
                                        { 0.0, 1.0, 0.0, 0.0},
                                        { 0.0, 0.0, 0.5, 0.5},
                                        { 0.5, 0.0, 0.5, 0.0},
                                        { 0.0, 0.5, 0.5, 0.0},
                                        { 0.0, 0.0, 1.0, 0.0}};

          // new triangles as connectivity information of the vertices baryc above
          const int tets[8][4] = { { 1, 2, 4, 7},  //corner x
                                   { 3, 4, 5, 8},  //corner y
                                   { 6, 7, 8, 9},  //corner z
                                   { 3, 4, 8, 6},  //prism part 1
                                   { 3, 4, 1, 6},  //prism part 2
                                   { 0, 1, 3, 6},  //prism part 3
                                   { 8, 6, 1, 7},  //pyramid part 1
                                   { 8, 4, 1, 7}}; //pyramid part 1

          for (int i = 0; i < 8; ++i) // tets
          {
            NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_i (*this, 1, 0);
            numint_i.SetVerticesTime(verts_time);
            Array< Vec<D> > newverts(4);
            for (int j = 0; j < 4; ++j) //vertices
            {
              newverts[j] = Vec<D>(0.0);
              for (int d = 0; d < 4; ++d)
                newverts[j] += baryc[tets[i][j]][d] * verts_space[d];
            }
            numint_i.SetVerticesSpace(newverts);
            if (ET_TIME == ET_POINT)
              numint_i.SetDistanceThreshold(0.5*distance_threshold);
            else
              numint_i.SetDistanceThreshold(distance_threshold);
            numint_i.MakeQuadRule(); // recursive call!
          }

        }
        else
            throw Exception(" refine_space in 3D in NumInt::MakeQuad not yet implemented");
      }

      if (!refine_space && !refine_time) // already on finest level: deal with cut situation
      {
        // RegionTimer reg (timer);
        Array<Simplex<SD> *> simplices;
        const int nvt = ET_TIME == ET_SEGM ? 2 : 1;
        const int nvs = verts_space.Size();
        Array<const Vec<SD> * > verts(nvs * nvt);
        for (int K = 0; K < nvt; ++K)
          for (int i = 0; i < nvs; ++i)
          {
            Vec<SD> newpoint;
            for (int j = 0; j < D; ++j)
            {
              newpoint[j] = verts_space[i][j];
              if (ET_TIME == ET_SEGM)
                newpoint[D] = K == 0 ? verts_time[0] : verts_time[verts_time.Size()-1];
            }
            verts[i+K*nvs] = pc(newpoint);
            // cout << "verts["<<i+K*nvs<<"]:" << newpoint << endl;
          }

        // spacetime: prism to simplices
        // only space: set simplix
        // in both cases we have a list of SD-dimensional simplices
        if (ET_TIME==ET_POINT)
        {
          simplices.SetSize(1);
          simplices[0] = new Simplex<SD>(verts);
        }
        else
        {
          DecomposePrismIntoSimplices<SD>(verts, simplices, pc, lh);
        }

        // cout << "simplices:\n";
        // for (int i = 0; i < simplices.Size(); ++i)
        // {
        //   cout << *simplices[i] << endl;
        // }

        const ScalarFieldEvaluator & eval (*lset);

        for (int i = 0; i < simplices.Size(); ++i)
        {
          // Check for each simplex if it is cut.
          // If yes call decomposition strategy for according dimension
          // If no  direction fill the composition rule accordingly
          DOMAIN_TYPE dt_simplex = simplices[i]->CheckIfCut(eval);
          if (dt_simplex == IF)
          {
            MakeQuadRuleOnCutSimplex<SD>(*simplices[i], *this);
          }
          else
          {
            FillSimplexWithRule<SD>(*simplices[i],
                                    compquadrule.GetRule(dt_simplex),
                                    GetIntegrationOrderMax());
            if (SD==2 && simplex_array_neg)
            {
              if (dt_simplex == NEG)
                simplex_array_neg->Append(new Simplex<SD> (*simplices[i]));
              else
                simplex_array_pos->Append(new Simplex<SD> (*simplices[i]));
            }
          }
          delete simplices[i];
        }
      }
      quaded = true;
      return IF;
    }
    else // no cut
    {
      //static Timer timer ("MakeQuadRule::FillUnCutSimplex");
      //RegionTimer reg (timer);
      // ThreadRegionTimer reg (timer, TaskManager::GetThreadId());

      double trafofac = 1.0;
      if (D==2)
      {
        Vec<D> a = verts_space[1] - verts_space[0];
        Vec<D> b = verts_space[2] - verts_space[0];
        trafofac = abs(a(0) * b(1) - a(1) * b(0));
        if (SD==2 && simplex_array_neg)
        {
          Array<const Vec<SD> *> simpl_verts(3);
          simpl_verts[0] = pc(verts_space[0]);
          simpl_verts[1] = pc(verts_space[1]);
          simpl_verts[2] = pc(verts_space[2]);
          if (dt_self == NEG)
            simplex_array_neg->Append(new Simplex<SD> (simpl_verts));
          else
            simplex_array_pos->Append(new Simplex<SD> (simpl_verts));
        }
      }
      else if (D==3)
      {
        Vec<D> a = verts_space[1] - verts_space[0];
        Vec<D> b = verts_space[2] - verts_space[0];
        Vec<D> c = verts_space[3] - verts_space[0];
        trafofac = abs(Determinant(a,b,c));
      }
      else // D==1
      {
        Vec<D> a = verts_space[1] - verts_space[0];
        trafofac = abs(a(0));
      }

      const double dt = verts_time[verts_time.Size()-1] - verts_time[0];
      const double t0 = verts_time[0];


      const IntegrationRule & ir_time = SelectIntegrationRule (ET_TIME, int_order_time);
      const IntegrationRule & ir_space = SelectIntegrationRule (ET_SPACE, int_order_space);


      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        double current = t0 + ir_time[l](0) * dt;
        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          Vec<D> point(0.0);
          double originweight = 1.0;
          for (int m = 0; m < D ;++m)
            originweight -= ir_space[k](m);
          point = originweight * verts_space[0];
          for (int m = 0; m < D ;++m)
            point += ir_space[k](m) * verts_space[m+1];
          const double weight = ir_time[l].Weight() * dt * ir_space[k].Weight() * trafofac;
          Vec<SD> ipoint(0.0);

          for (int m = 0; m < D ;++m)
            ipoint(m) = point(m);
          if (ET_trait<ET_TIME>::DIM > 0)
            ipoint(SD-1) = current;

          compquadrule.GetRule(dt_self).points.Append(ipoint);
          compquadrule.GetRule(dt_self).weights.Append(weight);

        }
      }
      quaded = true;
      return dt_self;
    }
  }

  template class NumericalIntegrationStrategy<ET_SEGM, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_SEGM, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;


  namespace DecompositionRules
  {
    template<int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <D> & s,
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    {
      cout << IM(1) << " ET_SPACE = " << ET_SPACE << ", ET_TIME = " << ET_TIME << endl;
      throw Exception("CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad --- no implementation for these Element Types");
    }


    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<3,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <3> & s,
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    {

      // {
      //   std::cout << "\n\n\n repeat inner <D=3> " << std::endl;

      //   for (int j = 0; j < 3+1; ++j)
      //   {
      //     Vec<3> test = *(s.p[j]);
      //     // test[2] = 0.0;
      //     std::cout << " (*numint.lset)(*(s.p[" << j << "])) = " << (*numint.lset)(test) << std::endl;
      //   }
      //   std::cout << " repeat again inner <D=3> \n \n " << std::endl;

      //   for (int j = 0; j < 3+1; ++j)
      //   {
      //     Vec<3> test = *(s.p[j]);
      //     // test[2] = 0.0;
      //     std::cout << " (*numint.lset)(*(s.p[" << j << "])) = " << (*numint.lset)(test) << std::endl;
      //   }
      //   std::cout << " \n\n\n " << std::endl;

      // }

      enum { SD = 3};

      Array< const Vec<SD> * > cutpoints(4);
      Array< const Vec<SD> * > pospoints(8);
      Array< const Vec<SD> * > negpoints(8);

      int ncutpoints = 0;
      int npospoints = 0;
      int nnegpoints = 0;

      Array<int> posvidx(0);
      Array<int> negvidx(0);

      // vertex idx connected to cut idx (just in case of 4 cut positions)
      // connectivity information of cuts
      Array<int> v2cut_1(4);
      Array<int> v2cut_2(4);
      v2cut_1 = -1;
      v2cut_2 = -1;

      const int edge[6][2] = { {0, 1},
                               {0, 2},
                               {0, 3},
                               {1, 2},
                               {1, 3},
                               {2, 3}};

      double vvals[4];
      bool zero[4];


      // cout << "\n\n\n vvals = \n";
      // for (int l = 0; l < 4; ++l)
      //   cout << l << ":" << (*numint.lset)(*(s.p[l])) << endl;

      for (int j = 0; j < 4; ++j)
      {
        zero[j] = false;
        vvals[j] = (*numint.lset)(*(s.p[j]));
        if (vvals[j] > 0)
        {
          pospoints[npospoints++] = numint.pc(*(s.p[j]));
          posvidx.Append(j);
        }
        else if (vvals[j] < 0)
        {
          negpoints[nnegpoints++] = numint.pc(*(s.p[j]));
          negvidx.Append(j);
        }
        else // (vvals[j] == 0.0)
        {
          pospoints[npospoints++] = numint.pc(*(s.p[j]));
          posvidx.Append(j);
          zero[j] = true;
        }
      }

      // cout << " vvals = \n";
      // for (int l = 0; l < 4; ++l)
      //   cout << l << ":" << vvals[l] << endl;

      int cntcuts = 0;
      for (int j = 0; j < 6; ++j)
      {
        const int lv = edge[j][0];
        const int rv = edge[j][1];
        const double valleft = vvals[lv];
        const double valright = vvals[rv];
        bool hascut = (valleft * valright < 0);
        if (zero[lv] && valright < 0)
          hascut = true;

        if (zero[rv] && valleft < 0)
          hascut = true;

        if (hascut)
        {
          const double cutpos = valleft / (valleft - valright);
          // std::cout << " cutpos = " << cutpos << std::endl;
          Vec<SD> p = (1-cutpos) * *(s.p[lv]) + cutpos * *(s.p[rv]) ;
          cutpoints[ncutpoints] = numint.pc(p);
          ncutpoints++;
          // collect connectivity of cut and vertices
          if (v2cut_1[lv] == -1)
            v2cut_1[lv] = cntcuts;
          else
            v2cut_2[lv] = cntcuts;
          if (v2cut_1[rv] == -1)
            v2cut_1[rv] = cntcuts;
          else
            v2cut_2[rv] = cntcuts;
          cntcuts ++;
        }
      }


      if (ncutpoints == 3) // three intersections: prism + tetra
      {
        Array< const Vec<SD> *> & minorgroup ( nnegpoints > npospoints ?
                                               pospoints : negpoints);
        int & nminorgroup ( nnegpoints > npospoints ? npospoints : nnegpoints);
        DOMAIN_TYPE dt_minor = nnegpoints > npospoints ? POS : NEG;
        Array< const Vec<SD> *> & majorgroup ( nnegpoints <= npospoints ?
                                               pospoints : negpoints);
        int & nmajorgroup ( nnegpoints <= npospoints ? npospoints : nnegpoints);
        DOMAIN_TYPE dt_major = nnegpoints <= npospoints ? POS : NEG;

        Array<int> & majvidx( negvidx.Size() > posvidx.Size() ? negvidx : posvidx);

        for (int k = 0; k < 3; ++k)
          minorgroup[nminorgroup++] = cutpoints[k];
        // minorgroup is a simplex of type dt_minor
        FillSimplexWithRule<SD>(minorgroup,
                                numint.compquadrule.GetRule(dt_minor),
                                numint.GetIntegrationOrderMax());

        Array< Simplex<SD> * > innersimplices(0);
        for (int k = 0; k < 3; ++k)
        {
          int corresponding_cut = v2cut_1[majvidx[k]];
          majorgroup[nmajorgroup++] = cutpoints[corresponding_cut];
        }
        DecomposePrismIntoSimplices<SD>(majorgroup, innersimplices, numint.pc, numint.lh);
        for (int l = 0; l < innersimplices.Size(); ++l)
        {
          FillSimplexWithRule<SD>(innersimplices[l]->p,
                                  numint.compquadrule.GetRule(dt_major),
                                  numint.GetIntegrationOrderMax());
          delete innersimplices[l];
        }

        // and the interface:
        FillSimplexCoDim1WithRule<SD> ( cutpoints, *pospoints[0],
                                        numint.compquadrule.GetInterfaceRule(),
                                        numint.GetIntegrationOrderMax());


      }
      else if (ncutpoints == 4) // four intersections: prism + prism
      {
        //pos domain
        {
          Array< const Vec<SD> *> posprism(6);
          posprism[0] = pospoints[0];
          const int idxn = posvidx[0];
          const int cut1 = v2cut_1[idxn];
          const int cut2 = v2cut_2[idxn];
          posprism[1] = cutpoints[cut1];
          posprism[2] = cutpoints[cut2];
          posprism[3] = pospoints[1];
          int cut3 = -1;
          for (int l = 0; l < 4; ++l)
            if (cut1 != l && cut2 != l)
              cut3 = l;
          int cut4 = 6 - cut3 - cut2 - cut1;

          // possibly switch orientation of cut3 / cut4
          if ((v2cut_1[negvidx[0]] == cut1 && v2cut_2[negvidx[0]] == cut4)
              || (v2cut_1[negvidx[0]] == cut2 && v2cut_2[negvidx[0]] == cut3))
          {
            const int cutt = cut3;
            cut3 = cut4;
            cut4 = cutt;
          }
          posprism[4] = cutpoints[cut3];
          posprism[5] = cutpoints[cut4];

          // cout << " cutpoints vertices: " << endl;
          // for (int l = 0; l < 4; ++l)
          //   cout << *cutpoints[l] << endl;

          // cout << " pospoints vertices: " << endl;
          // for (int l = 0; l < 2; ++l)
          //   cout << *pospoints[l] << endl;

          // cout << " posprism vertices: " << endl;
          // for (int l = 0; l < 6; ++l)
          //   cout << *posprism[l] << endl;

          Array< Simplex<SD> * > innersimplices(0);
          DecomposePrismIntoSimplices<SD>(posprism, innersimplices, numint.pc, numint.lh);
          for (int l = 0; l < innersimplices.Size(); ++l)
          {
            // std::cout << " *innersimplices[l] = " << *innersimplices[l] << std::endl;
            FillSimplexWithRule<SD>(innersimplices[l]->p,
                                    numint.compquadrule.GetRule(POS),
                                    numint.GetIntegrationOrderMax());
            delete innersimplices[l];
          }
        }
        //neg domain
        {
          Array< const Vec<SD> *> negprism(6);
          negprism[0] = negpoints[0];
          const int idxn = negvidx[0];
          const int cut1 = v2cut_1[idxn];
          const int cut2 = v2cut_2[idxn];
          negprism[1] = cutpoints[cut1];
          negprism[2] = cutpoints[cut2];
          negprism[3] = negpoints[1];
          int cut3 = -1;
          for (int l = 0; l < 4; ++l)
            if (cut1 != l && cut2 != l)
              cut3 = l;
          int cut4 = 6 - cut3 - cut2 - cut1;

          // possibly switch orientation of cut3 / cut4
          if ((v2cut_1[posvidx[0]] == cut1 && v2cut_2[posvidx[0]] == cut4)
              || (v2cut_1[posvidx[0]] == cut2 && v2cut_2[posvidx[0]] == cut3))
          {
            const int cutt = cut3;
            cut3 = cut4;
            cut4 = cutt;
          }
          negprism[4] = cutpoints[cut3];
          negprism[5] = cutpoints[cut4];

          Array< Simplex<SD> * > innersimplices(0);
          DecomposePrismIntoSimplices<SD>(negprism, innersimplices, numint.pc, numint.lh);
          for (int l = 0; l < innersimplices.Size(); ++l)
          {
            FillSimplexWithRule<SD>(innersimplices[l]->p,
                                    numint.compquadrule.GetRule(NEG),
                                    numint.GetIntegrationOrderMax());
            delete innersimplices[l];
          }
        }
        //interface
        {
          int diag1, diag2;
          int ndiag1, ndiag2;
          if (v2cut_1[negvidx[0]] == v2cut_1[posvidx[0]] || v2cut_2[negvidx[0]] == v2cut_1[posvidx[0]])
          {
            diag1 = v2cut_1[posvidx[0]];
            ndiag1 = v2cut_2[posvidx[0]];
          }
          else
          {
            diag1 = v2cut_2[posvidx[0]];
            ndiag1 = v2cut_1[posvidx[0]];
          }

          if (v2cut_1[negvidx[1]] == v2cut_1[posvidx[1]] || v2cut_2[negvidx[1]] == v2cut_1[posvidx[1]])
          {
            diag2 = v2cut_1[posvidx[1]];
            ndiag2 = v2cut_2[posvidx[1]];
          }
          else
          {
            diag2 = v2cut_2[posvidx[1]];
            ndiag2 = v2cut_1[posvidx[1]];
          }

          Array< const Vec<SD> * > trig1(3);
          Array< const Vec<SD> * > trig2(3);

          trig1[0] = cutpoints[diag1];
          trig1[1] = cutpoints[ndiag1];
          trig1[2] = cutpoints[diag2];

          trig2[0] = cutpoints[diag2];
          trig2[1] = cutpoints[ndiag2];
          trig2[2] = cutpoints[diag1];

          FillSimplexCoDim1WithRule<SD> ( trig1, *pospoints[0],
                                          numint.compquadrule.GetInterfaceRule(),
                                          numint.GetIntegrationOrderMax());
          FillSimplexCoDim1WithRule<SD> ( trig2, *pospoints[0],
                                          numint.compquadrule.GetInterfaceRule(),
                                          numint.GetIntegrationOrderMax());
        }

      } // end of 3 or 4 cutpoints
      else
      {
        cout << "ncutpoints = " << ncutpoints << endl;
        throw Exception(" did not expect this.. -3-");
      }
    }

    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<2,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <2> & s,
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    {
      enum { SD = 2};

      Array< const Vec<SD> * > cutpoints(0);
      Array< const Vec<SD> * > pospoints(0);
      Array< const Vec<SD> * > negpoints(0);

      Array<int> posvidx(0);
      Array<int> negvidx(0);

      // vertex idx connected to cut idx (just in case of 4 cut positions)
      // connectivity information of cuts
      Array<int> v2cut_1(4);
      Array<int> v2cut_2(4);
      v2cut_1 = -1;
      v2cut_2 = -1;

      const int edge[3][2] = { {0, 1},
                               {0, 2},
                               {1, 2}};

      double vvals[3];
      bool zero[3];

      for (int j = 0; j < 3; ++j)
      {
        zero[j] = false;
        vvals[j] = (*numint.lset)(*(s.p[j]));
        if (vvals[j] > 0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
        }
        else if (vvals[j] < 0)
        {
          negpoints.Append(numint.pc(*(s.p[j])));
          negvidx.Append(j);
        }
        else // (vvals[j] == 0.0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
          zero[j] = true;
        }
      }

      // cout << " Avvals = \n";
      // for (int i = 0; i < 3; ++i)
      //   cout << i << ":" << vvals[i] << endl;

      int cntcuts = 0;
      for (int j = 0; j < 3; ++j) //edges
      {
        const int lv = edge[j][0];
        const int rv = edge[j][1];
        const double valleft = vvals[lv];
        const double valright = vvals[rv];
        bool hascut = (valleft * valright < 0);
        if (zero[lv] && valright < 0)
          hascut = true;

        if (zero[rv] && valleft < 0)
          hascut = true;

        if (hascut)
        {
          const double cutpos = valleft / (valleft - valright);
          // std::cout << " cutpos = " << cutpos << std::endl;
          Vec<SD> p = (1-cutpos) * *(s.p[lv]) + cutpos * *(s.p[rv]) ;
          cutpoints.Append(numint.pc(p));
          // collect connectivity of cut and vertices
          if (v2cut_1[lv] == -1)
            v2cut_1[lv] = cntcuts;
          else
            v2cut_2[lv] = cntcuts;
          if (v2cut_1[rv] == -1)
            v2cut_1[rv] = cntcuts;
          else
            v2cut_2[rv] = cntcuts;
          cntcuts ++;
        }
      }

      // std::cout << " cutpoints[0] = " << *cutpoints[0] << std::endl;
      // std::cout << " cutpoints[1] = " << *cutpoints[1] << std::endl;

      if (cutpoints.Size() == 2) // three intersections: prism + tetra
      {

        Array< const Vec<SD> *> & minorgroup ( negpoints.Size() > pospoints.Size() ?
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_minor = negpoints.Size() > pospoints.Size() ? POS : NEG;
        Array< const Vec<SD> *> & majorgroup ( negpoints.Size() <= pospoints.Size() ?
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_major = negpoints.Size() <= pospoints.Size() ? POS : NEG;

        Array<int> & majvidx( negvidx.Size() > posvidx.Size() ? negvidx : posvidx);

        for (int k = 0; k < 2; ++k)
          minorgroup.Append(cutpoints[k]);
        // minorgroup is a simplex of type dt_minor
        FillSimplexWithRule<SD>(minorgroup,
                                numint.compquadrule.GetRule(dt_minor),
                                numint.GetIntegrationOrderMax());

        // for result visualization
        if (numint.simplex_array_neg && (dt_minor == NEG))
            numint.simplex_array_neg->Append(new Simplex<SD> (minorgroup));
        if (numint.simplex_array_pos && (dt_minor == POS))
            numint.simplex_array_pos->Append(new Simplex<SD> (minorgroup));

        Array< Simplex<SD> * > innersimplices(0);
        for (int k = 0; k < 2; ++k)
        {
          int corresponding_cut = v2cut_1[majvidx[k]];
          majorgroup.Append(cutpoints[corresponding_cut]);
        }

        DecomposePrismIntoSimplices<SD>(majorgroup, innersimplices, numint.pc, numint.lh);
        for (int l = 0; l < innersimplices.Size(); ++l)
        {
          FillSimplexWithRule<SD>(innersimplices[l]->p,
                                  numint.compquadrule.GetRule(dt_major),
                                  numint.GetIntegrationOrderMax());

          // for result visualization
          if (numint.simplex_array_neg && (dt_minor == POS))
            numint.simplex_array_neg->Append(new Simplex<SD> (innersimplices[l]->p));
          if (numint.simplex_array_pos && (dt_minor == NEG))
            numint.simplex_array_pos->Append(new Simplex<SD> (innersimplices[l]->p));

          delete innersimplices[l];
        }

        // and the interface:
        FillSimplexCoDim1WithRule<SD> ( cutpoints, *pospoints[0],
                                        numint.compquadrule.GetInterfaceRule(),
                                        numint.GetIntegrationOrderMax());


      }
      else
      {
        cout << IM(1) << "cutpoints.Size() = " << cutpoints.Size() << endl;
        cout << IM(1) << " Avvals = \n";
        for (int i = 0; i < 3; ++i)
          cout << IM(1) << i << ":" << vvals[i] << endl;
        throw Exception(" did not expect this.. -2-");
      }
    }


    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<1,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <1> & s,
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    {
      enum { SD = 1};

      // static Timer timer ("CutSimplex<1>::MakeQuad");
      // ThreadRegionTimer reg (timer, TaskManager::GetThreadId());
      // RegionTimer reg (timer);

      const Vec<1> & left = *(s.p[0]);
      const Vec<1> & right = *(s.p[1]);

      double valleft = (*numint.lset)(left);
      double valright = (*numint.lset)(right);

      const double cutpos = valleft / (valleft - valright);
      Vec<SD> mid = (1-cutpos) * left + cutpos * right ;
      const Vec<SD> * p = numint.pc(mid);

      Array < const Vec<SD> * > leftint(2);
      leftint[0] = s.p[0]; leftint[1] = p;

      Array < const Vec<SD> * > rightint(2);
      rightint[0] = p; rightint[1] = s.p[1];

      Simplex<1> leftsimplex (leftint);
      Simplex<1> rightsimplex (rightint);

      DOMAIN_TYPE dt_left, dt_right;

      if (valleft > 0)
      {
        dt_left = POS;
        dt_right = NEG;
      }
      else
      {
        dt_left = NEG;
        dt_right = POS;
      }


      FillSimplexWithRule<SD>(leftsimplex,
                              numint.compquadrule.GetRule(dt_left),
                              numint.GetIntegrationOrderMax());

      FillSimplexWithRule<SD>(rightsimplex,
                              numint.compquadrule.GetRule(dt_right),
                              numint.GetIntegrationOrderMax());

    }


  } // end of namespace DecompositionRules



  template <int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void MakeQuadRuleOnCutSimplex(const Simplex <D> & s,
                                const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {
    // RegionTimer reg (timer);
    // std::cout << " from here MakeQuadRuleOnCutSimplex "<< D << " " << ET_SPACE << " " << ET_TIME  << std::endl;
    // std::cout << " simplex s = " << s << std::endl;
    // if (D==3)
    // {

    //   for (int j = 0; j < D+1; ++j)
    //     std::cout << " (*numint.lset)(*(s.p[" << j << "])) = " << (*numint.lset)(*(s.p[j])) << std::endl;
    //   std::cout << " repeat ... " << std::endl;

    //   for (int j = 0; j < D+1; ++j)
    //     std::cout << " (*numint.lset)(*(s.p[" << j << "])) = " << (*numint.lset)(*(s.p[j])) << std::endl;
    // }
    DecompositionRules::CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad(s,numint);
  }


  // integration rules that are returned assume that a scaling with mip.GetMeasure() gives the
  // correct weight on the "physical" domain (note that this is not a natural choicefor interface integrals)
  const IntegrationRule * CutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                             const ElementTransformation & trafo,
                                             DOMAIN_TYPE dt,
                                             int intorder,
                                             int subdivlvl,
                                             LocalHeap & lh)
  {
    static int timer = NgProfiler::CreateTimer ("OldCutIntegrationRule"); NgProfiler::RegionTimer reg (timer);

    int DIM = trafo.SpaceDim();
    auto lset_eval
      = ScalarFieldEvaluator::Create(DIM,*cf_lset,trafo,lh);

    if (trafo.VB() == BND)
      DIM--;
      // tstart.Stop();

    auto et = trafo.GetElementType();

    shared_ptr<XLocalGeometryInformation> xgeom = nullptr;

    CompositeQuadratureRule<1> cquad1d;
    CompositeQuadratureRule<2> cquad2d;
    CompositeQuadratureRule<3> cquad3d;
    if (DIM == 1)
      xgeom = XLocalGeometryInformation::Create(et, ET_POINT,
                                                *lset_eval, cquad1d, lh,
                                                intorder, 0, subdivlvl, 0);
    else if (DIM == 2)
      xgeom = XLocalGeometryInformation::Create(et, ET_POINT,
                                                *lset_eval, cquad2d, lh,
                                                intorder, 0, subdivlvl, 0);
    else
      xgeom = XLocalGeometryInformation::Create(et, ET_POINT,
                                                *lset_eval, cquad3d, lh,
                                                intorder, 0, subdivlvl, 0);
    DOMAIN_TYPE element_domain = xgeom->MakeQuadRule();

    const IntegrationRule* ir = nullptr;

    if (element_domain == IF) // there is a cut on the current element
    {
      if (dt == IF)
      {
        if (DIM == 1)
        {
          throw Exception("no interface quad rule for 1D for now...");
        }
        else if (DIM == 2)
        {
          const QuadratureRuleCoDim1<2> & interface_quad(cquad2d.GetInterfaceRule());
          IntegrationRule * ir_interface  = new (lh) IntegrationRule(interface_quad.Size(),lh);
          for (int i = 0; i < interface_quad.Size(); ++i)
          {
            IntegrationPoint ip(&interface_quad.points[i](0),interface_quad.weights[i]);
            MappedIntegrationPoint<2,2> mip(ip,trafo);

            Mat<2,2> Finv = mip.GetJacobianInverse();
            const double absdet = mip.GetMeasure();

            Vec<2> nref = interface_quad.normals[i];
            Vec<2> normal = absdet * Trans(Finv) * nref ;
            double len = L2Norm(normal);
            // const double weight = interface_quad.weights[i] * len;

            (*ir_interface)[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i] * len / mip.GetMeasure());
            ir = ir_interface;
          }
        }
        else
        {
          const QuadratureRuleCoDim1<3> & interface_quad(cquad3d.GetInterfaceRule());
          IntegrationRule * ir_interface  = new (lh) IntegrationRule(interface_quad.Size(),lh);
          for (int i = 0; i < interface_quad.Size(); ++i)
          {
            IntegrationPoint ip(&interface_quad.points[i](0),interface_quad.weights[i]);
            MappedIntegrationPoint<3,3> mip(ip,trafo);

            Mat<3,3> Finv = mip.GetJacobianInverse();
            const double absdet = mip.GetMeasure();

            Vec<3> nref = interface_quad.normals[i];
            Vec<3> normal = absdet * Trans(Finv) * nref ;
            double len = L2Norm(normal);
            // const double weight = interface_quad.weights[i] * len;

            (*ir_interface)[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i] * len / mip.GetMeasure());
            ir = ir_interface;
          }
        }
      }
      else
      {
        if (DIM == 1)
        {
          const QuadratureRule<1> & domain_quad = cquad1d.GetRule(dt);
          auto ir_domain = new (lh) IntegrationRule (domain_quad.Size(),lh);
          for (int i = 0; i < ir_domain->Size(); ++i)
            (*ir_domain)[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
          ir = ir_domain;
        }
        else if (DIM == 2)
        {
          const QuadratureRule<2> & domain_quad = cquad2d.GetRule(dt);
          auto ir_domain = new (lh) IntegrationRule (domain_quad.Size(),lh);
          for (int i = 0; i < ir_domain->Size(); ++i)
            (*ir_domain)[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
          ir = ir_domain;
        }
        else
        {
          const QuadratureRule<3> & domain_quad = cquad3d.GetRule(dt);
          auto ir_domain = new (lh) IntegrationRule (domain_quad.Size(),lh);
          for (int i = 0; i < ir_domain->Size(); ++i)
            (*ir_domain)[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
          ir = ir_domain;
        }
      }
    }
    else
    {
      if (element_domain != dt) //no integration on this element
        return nullptr;
      ir = & (SelectIntegrationRule (trafo.GetElementType(), intorder));
    }

    return ir;
  }

  template FlatArray<SIMD<double>> CreateSIMD_FlatArray(FlatArray<double> ns_arr, LocalHeap & lh);
  //template FlatArray<SIMD<Complex>> CreateSIMD_FlatArray(FlatArray<Complex> ns_arr, LocalHeap & lh);


} // end of namespace
