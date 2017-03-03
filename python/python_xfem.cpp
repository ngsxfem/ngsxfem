//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../utils/bitarraycf.hpp"
#include "../xfem/cutinfo.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/ghostpenalty.hpp"
#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../cutint/straightcutrule.hpp"
// #include "../utils/error.hpp"

//using namespace ngcomp;

void ExportNgsx(py::module &m)
{
  



  typedef PyWrapper<FESpace> PyFES;
  typedef PyWrapper<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef PyWrapper<GF> PyGF;
  typedef PyWrapper<shared_ptr<BitArray>> PyBA;

  py::enum_<DOMAIN_TYPE>(m, "DOMAIN_TYPE")
  .value("POS", POS)
  .value("NEG", NEG)
  .value("IF", IF)
  .export_values()
  ;

  // typedef PyWrapperDerived<CompoundFESpace, FESpace> PyCompFES;

  typedef PyWrapperDerived<XFESpace, FESpace> PyXFES;
  typedef PyWrapperDerived<XStdFESpace, FESpace> PyXStdFES;

  m.def("XToNegPos", FunctionPointer( [] (PyGF gfx, PyGF gfnegpos) {
    XFESpace::XToNegPos(gfx.Get(),gfnegpos.Get());
  } ) );

  py::class_<CutInformation>
    (m, "CutInfo")
  .def("__init__", FunctionPointer( [] (CutInformation *instance,
                                        shared_ptr<MeshAccess> ma,
                                        py::object lset,
                                        int heapsize)
  {
    new (instance) CutInformation (ma);
    if (py::extract<PyCF> (lset).check())
    {
      PyCF cflset = py::extract<PyCF>(lset)();
      LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
      instance->Update(cflset.Get(),lh);
    }
  }),
       py::arg("mesh"),
       py::arg("levelset") = DummyArgument(),
       py::arg("heapsize") = 1000000
       )
  .def("Update", FunctionPointer ([](CutInformation & self,
                                     PyCF lset,
                                     int heapsize)
  {
    LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
    self.Update(lset.Get(),lh);
  }),
       py::arg("levelset"),
       py::arg("heapsize") = 1000000
       )
  .def("Mesh", FunctionPointer ([](CutInformation & self)
  {
    return self.GetMesh();
  })
       )
  .def("GetElementsOfType", FunctionPointer ([](CutInformation & self,
                                                DOMAIN_TYPE dt,
                                                VorB vb)
  {
    return self.GetElementsOfDomainType(dt,vb);
  }),
       py::arg("domain_type") = IF,
       py::arg("VOL_or_BND") = VOL
       )
  .def("GetFacetsOfType", FunctionPointer ([](CutInformation & self,
                                              DOMAIN_TYPE dt)
  {
    return self.GetFacetsOfDomainType(dt);
  }),
       py::arg("domain_type") = IF)

  .def("GetCutRatios", FunctionPointer ([](CutInformation & self,
                                           VorB vb)
  {
    return self.GetCutRatios(vb);
  }),
       py::arg("VOL_or_BND") = VOL)
  ;


  m.def("GetFacetsWithNeighborTypes",
        FunctionPointer( [] (shared_ptr<MeshAccess> ma,
                             shared_ptr<BitArray> a,
                             bool bv_a,
                             bool bv_b,
                             bool use_and,
                             py::object bb,
                             int heapsize)
  {
    LocalHeap lh (heapsize, "FacetsWithNeighborTypes-heap", true);
    shared_ptr<BitArray> b = nullptr;
    if (py::extract<PyBA> (bb).check())
      b = py::extract<PyBA>(bb)();
    else
      b = a;
    return GetFacetsWithNeighborTypes(ma,a,b,bv_a,bv_b,use_and,lh);
  } ),
        py::arg("mesh"),
        py::arg("a"),
        py::arg("bnd_val_a") = true,
        py::arg("bnd_val_b") = true,
        py::arg("use_and") = true,
        py::arg("b") = DummyArgument(),
        py::arg("heapsize") = 1000000
        );

  m.def("GetElementsWithNeighborFacets",
        FunctionPointer( [] (shared_ptr<MeshAccess> ma,
                             shared_ptr<BitArray> a,
                             int heapsize)
  {
    LocalHeap lh (heapsize, "GetElementsWithNeighborFacets-heap", true);
    return GetElementsWithNeighborFacets(ma,a,lh);
  } ),
        py::arg("mesh"),
        py::arg("a"),
        py::arg("heapsize") = 1000000
        );


  // Export RandomCoefficientFunction to python (name "RandomCF")
  typedef PyWrapperDerived<BitArrayCoefficientFunction,CoefficientFunction> PyBACF;
  py::class_<PyBACF, PyCF>
    (m, "BitArrayCF")
  .def("__init__",
       [](PyBACF *instance, shared_ptr<BitArray> ba)
  {
    new (instance) PyBACF(make_shared<BitArrayCoefficientFunction> (ba));
  },
       py::arg("bitarray")
       );

  py::class_<PyXFES, PyFES>
    (m, "XFESpace")
  .def("__init__", FunctionPointer( [] (PyXFES *instance,
                                        PyFES basefes,
                                        PyCF lset,
                                        py::dict bpflags,
                                        int heapsize)
  {
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<XFESpace> ret = nullptr;
    shared_ptr<MeshAccess> ma = basefes.Get()->GetMeshAccess();
    if (ma->GetDimension()==2)
      ret = make_shared<T_XFESpace<2>> (ma, basefes.Get(), lset.Get(), flags);
    else
      ret = make_shared<T_XFESpace<3>> (ma, basefes.Get(), lset.Get(), flags);
    LocalHeap lh (heapsize, "XFESpace::Update-heap", true);
    ret->Update(lh);
    new (instance) PyFES(ret);
  }),
       py::arg("basefes"),
       py::arg("lset"),
       py::arg("flags") = py::dict(),
       py::arg("heapsize") = 1000000)

  .def("SetLevelSet", FunctionPointer ([](PyXFES self, PyCF cf)
  {
    self.Get()->SetLevelSet(cf.Get());
  }),
       "Update information on level set function")
  .def("SetLevelSet", FunctionPointer ([](PyXFES self, PyGF gf)
  {
    self.Get()->SetLevelSet(gf.Get());
  }),
       "Update information on level set function")
  .def("SetBaseFESpace", FunctionPointer ([](PyXFES self, PyFES fes)
  {
    self.Get()->SetBaseFESpace(fes.Get());
  }),
       "Update information on base FESpace")
  .def("BaseDofOfXDof", FunctionPointer ([](PyXFES self, int i)
  {
    return self.Get()->GetBaseDofOfXDof(i);
  }),
       "get corresponding dof of base FESpace")
  .def("GetNVertexDofs", FunctionPointer ([](PyXFES self)
  {
    return self.Get()->GetNVertexDof();
  }),
       "get number of x dofs at vertices")
  .def("CutElements", FunctionPointer ([](PyXFES self)
  {
    return self.Get()->CutElements();
  }),
       "get BitArray of cut elements")
  .def("CutSurfaceElements", FunctionPointer ([](PyXFES self)
  {
    return self.Get()->CutSurfaceElements();
  }),
       "get BitArray of cut surface elements")
  .def("GetDomainOfDof", FunctionPointer ([](PyXFES self, int i)
  {
    return self.Get()->GetDomainOfDof(i);
  }),
       "get domain_type of degree of freedom")
  .def("GetDomainOfElement", FunctionPointer ([](PyXFES self, int i)
  {
    return self.Get()->GetDomainOfElement(i);
  }),
       "get domain_type of element")
  .def("GetDomainNrs",  FunctionPointer( [] (PyXFES self, int elnr) {
    Array<DOMAIN_TYPE> domnums;
    self.Get()->GetDomainNrs( elnr, domnums );
    return domnums;
  }))
  ;

  m.def("InterpolateToP1", FunctionPointer( [] (PyGF gf_ho, PyGF gf_p1, int heapsize)
  {
    InterpolateP1 interpol(gf_ho.Get(), gf_p1.Get());
    LocalHeap lh (heapsize, "InterpolateP1-Heap");
    interpol.Do(lh);
  } ),
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,py::arg("heapsize")=1000000)
  ;

  m.def("InterpolateToP1", FunctionPointer( [] (PyCF coef, PyGF gf_p1, int heapsize)
  {
    InterpolateP1 interpol(coef.Get(), gf_p1.Get());
    LocalHeap lh (heapsize, "InterpolateP1-Heap");
    interpol.Do(lh);
  } ),
        py::arg("coef"),py::arg("gf"),py::arg("heapsize")=1000000)
  ;

  py::class_<StatisticContainer, shared_ptr<StatisticContainer>>(m, "StatisticContainer")
  .def(py::init<>())
  .def("Print", FunctionPointer ([](StatisticContainer & self, string label, string select)
  {
    if (select == "L1")
      PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
    if (select == "L2")
      PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
    if (select == "max")
      PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
    if (select == "misc")
      PrintConvergenceTable(self.ErrorMisc,label+"_misc");
    if (select == "all")
    {
      PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
      PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
      PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
      PrintConvergenceTable(self.ErrorMisc,label+"_misc");
    }
  }),
       py::arg("label")="something",py::arg("select")="all"
       )
  ;

  m.def("CalcMaxDistance", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, int heapsize)
  {
    StatisticContainer dummy;
    LocalHeap lh (heapsize, "CalcDistance-Heap");
    if (lset_p1->GetMeshAccess()->GetDimension()==2)
      CalcDistances<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  dummy, lh, -1.0, false);
    else
      CalcDistances<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  dummy, lh, -1.0, false);
    return (double) dummy.ErrorMaxNorm[dummy.ErrorMaxNorm.Size()-1];
  } ),
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("heapsize")=1000000)
  ;

  m.def("CalcDistances", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
  {
    LocalHeap lh (heapsize, "CalcDistance-Heap");
    if (lset_p1.Get()->GetMeshAccess()->GetDimension()==2)
      CalcDistances<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  stats, lh, refine_threshold, absolute);
    else
      CalcDistances<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  stats, lh, refine_threshold, absolute);
  } ),
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("stats")=NULL,py::arg("heapsize")=1000000,py::arg("refine_threshold")=-1.0,py::arg("absolute")=false)
  ;

  m.def("CalcDeformationError", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, StatisticContainer & stats, double lower, double upper, int heapsize)
  {
    LocalHeap lh (heapsize, "CalcDeformationError-Heap");
    if (lset_p1.Get()->GetMeshAccess()->GetDimension()==2)
      CalcDeformationError<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), stats, lh, lower, upper);
    else
      CalcDeformationError<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), stats, lh, lower, upper);
  } ),
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL,py::arg("stats")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
  ;

  m.def("ProjectShift", FunctionPointer( [] (PyGF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, double lower, double upper, double threshold, int heapsize)
  {
    LocalHeap lh (heapsize, "ProjectShift-Heap");
    ProjectShift(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), lower, upper, threshold, lh);
  } ),
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("threshold")=1.0,py::arg("heapsize")=1000000)
  ;

// ProjectShift


  m.def("RefineAtLevelSet", FunctionPointer( [] (PyGF lset_p1, double lower, double upper, int heapsize)
  {
    LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
    RefineAtLevelSet(lset_p1.Get(), lower, upper, lh);
  } ),
        py::arg("lset_p1")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
  ;



  m.def("RefineAtLevelSet", FunctionPointer( [] (PyGF gf, double lower_lset_bound, double upper_lset_bound, int heapsize)
  {
    LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
    RefineAtLevelSet(gf.Get(),lower_lset_bound,upper_lset_bound,lh);
  } ),
        py::arg("gf"),py::arg("lower_lset_bound")=0.0,py::arg("upper_lset_bound")=0.0,py::arg("heapsize")=10000000)
  ;


  m.def("IntegrateX",
        FunctionPointer([](PyCF lset,
                           shared_ptr<MeshAccess> ma,
                           PyCF cf_neg,
                           PyCF cf_pos,
                           PyCF cf_interface,
                           int order, int subdivlvl, py::dict domains, int heapsize)
  {
    static Timer timer ("IntegrateX");
    static Timer timercutgeom ("IntegrateX::MakeCutGeom");
    static Timer timerevalcoef ("IntegrateX::EvalCoef");
    RegionTimer reg (timer);
    LocalHeap lh(heapsize, "lh-Integrate");

    Flags flags = py::extract<Flags> (domains)();

    Array<bool> tointon(3); tointon = false;
    Array<shared_ptr<CoefficientFunction>> cf(3); cf = nullptr;
    if (flags.GetDefineFlag("negdomain")) {
      tointon[int(NEG)] = true;
      if (cf_neg.Get() == nullptr)
        throw Exception("no coef for neg domain given");
      else
        cf[int(NEG)] = cf_neg.Get();
    }
    if (flags.GetDefineFlag("posdomain")) {
      tointon[int(POS)] = true;
      if (cf_pos.Get() == nullptr)
        throw Exception("no coef for pos domain given");
      else
        cf[int(POS)] = cf_pos.Get();
    }
    if (flags.GetDefineFlag("interface")) {
      tointon[int(IF)] = true;
      if (cf_interface.Get() == nullptr)
        throw Exception("no coef for interface domain given");
      else
        cf[int(IF)] = cf_interface.Get();
    }

    Vector<> domain_sum(3);                         // [val_neg,val_pos,val_interface]
    domain_sum = 0.0;
    int DIM = ma->GetDimension();
    ma->IterateElements
      (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
    {
      auto & trafo = ma->GetTrafo (el, lh);
      auto lset_eval
        = ScalarFieldEvaluator::Create(DIM,*(lset.Get()),trafo,lh);
      ELEMENT_TYPE eltype = el.GetType();
      timercutgeom.Start();
      shared_ptr<XLocalGeometryInformation> xgeom = nullptr;

      CompositeQuadratureRule<2> cquad2d;
      CompositeQuadratureRule<3> cquad3d;
      if (DIM == 2)
        xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT,
                                                  *lset_eval, cquad2d, lh,
                                                  order, 0, subdivlvl, 0);
      else
        xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT,
                                                  *lset_eval, cquad3d, lh,
                                                  order, 0, subdivlvl, 0);
      DOMAIN_TYPE element_domain = xgeom->MakeQuadRule();
      timercutgeom.Stop();
      if (element_domain == IF)
      {
        for (auto domtype : {NEG,POS})
        {
          if( ! tointon[int(domtype)] ) continue;
          double hsum = 0.0;
          // double value = 0.0;

          if (DIM == 2)
          {
            const QuadratureRule<2> & domain_quad = cquad2d.GetRule(domtype);
            IntegrationRule ir_domain (domain_quad.Size(),lh);
            for (int i = 0; i < ir_domain.Size(); ++i)
              ir_domain[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);

            MappedIntegrationRule<2,2> mir_domain(ir_domain, trafo,lh);
            FlatMatrix<> values(mir_domain.Size(), 1, lh);
            timerevalcoef.Start();
            cf[int(domtype)] -> Evaluate (mir_domain, values);
            timerevalcoef.Stop();

            for (int i = 0; i < domain_quad.Size(); ++i)
              hsum += values(i,0) * mir_domain[i].GetWeight();
          }
          else
          {
            const QuadratureRule<3> & domain_quad = cquad3d.GetRule(domtype);
            IntegrationRule ir_domain (domain_quad.Size(),lh);
            for (int i = 0; i < ir_domain.Size(); ++i)
              ir_domain[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);

            MappedIntegrationRule<3,3> mir_domain(ir_domain, trafo,lh);
            FlatMatrix<> values(mir_domain.Size(), 1, lh);
            timerevalcoef.Start();
            cf[int(domtype)] -> Evaluate (mir_domain, values);
            timerevalcoef.Stop();

            for (int i = 0; i < domain_quad.Size(); ++i)
              hsum += values(i,0) * mir_domain[i].GetWeight();
          }
          double & rsum = domain_sum(int(domtype));
          AsAtomic(rsum) += hsum;
        }

        if (tointon[int(IF)])
        {
          double hsum_if = 0.0;
          double value_if = 0.0;
          if (DIM == 2)
          {
            const QuadratureRuleCoDim1<2> & interface_quad(cquad2d.GetInterfaceRule());
            IntegrationRule ir_interface (interface_quad.Size(),lh);
            for (int i = 0; i < ir_interface.Size(); ++i)
              ir_interface[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i]);

            MappedIntegrationRule<2,2> mir_interface(ir_interface, trafo,lh);
            FlatMatrix<> values(mir_interface.Size(), 1, lh);
            timerevalcoef.Start();
            cf[int(IF)] -> Evaluate (mir_interface, values);
            timerevalcoef.Stop();

            for (int i = 0; i < interface_quad.Size(); ++i)
            {
              MappedIntegrationPoint<2,2> & mip(mir_interface[i]);

              Mat<2,2> Finv = mip.GetJacobianInverse();
              const double absdet = mip.GetMeasure();

              Vec<2> nref = interface_quad.normals[i];
              Vec<2> normal = absdet * Trans(Finv) * nref;
              double len = L2Norm(normal);
              const double weight = interface_quad.weights[i] * len;

              hsum_if += values(i,0) * weight;
            }
          }
          else
          {

            const QuadratureRuleCoDim1<3> & interface_quad(cquad3d.GetInterfaceRule());
            IntegrationRule ir_interface (interface_quad.Size(),lh);
            for (int i = 0; i < ir_interface.Size(); ++i)
              ir_interface[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i]);

            MappedIntegrationRule<3,3> mir_interface(ir_interface, trafo,lh);
            FlatMatrix<> values(mir_interface.Size(), 1, lh);
            timerevalcoef.Start();
            cf[int(IF)] -> Evaluate (mir_interface, values);
            timerevalcoef.Stop();

            for (int i = 0; i < interface_quad.Size(); ++i)
            {
              MappedIntegrationPoint<3,3> & mip(mir_interface[i]);

              Mat<3,3> Finv = mip.GetJacobianInverse();
              const double absdet = mip.GetMeasure();

              Vec<3> nref = interface_quad.normals[i];
              Vec<3> normal = absdet * Trans(Finv) * nref;
              double len = L2Norm(normal);
              const double weight = interface_quad.weights[i] * len;

              hsum_if += values(i,0) * weight;
            }
          }

          double & rsum_if = domain_sum(int(IF));
          AsAtomic(rsum_if) += hsum_if;
        }
      }
      else if( tointon[int(element_domain)] )
      {
        double hsum = 0.0;
        IntegrationRule ir(trafo.GetElementType(), order);
        BaseMappedIntegrationRule & mir = trafo(ir, lh);

        FlatMatrix<> values(ir.Size(), 1, lh);
        timerevalcoef.Start();
        cf[int(element_domain)] -> Evaluate (mir, values);
        timerevalcoef.Stop();
        for (int i = 0; i < values.Height(); i++)
          hsum += mir[i].GetWeight() * values(i,0);

        double & rsum = domain_sum(int(element_domain));
        AsAtomic(rsum) += hsum;
      }

    });

    py::dict resdict;
    if (tointon[int(NEG)])
      resdict["negdomain"] = py::cast(domain_sum(int(NEG)));
    if (tointon[int(POS)])
      resdict["posdomain"] = py::cast(domain_sum(int(POS)));
    if (tointon[int(IF)])
      resdict["interface"] = py::cast(domain_sum(int(IF)));

    return resdict;
    // bp::object result;
    // return  bp::list(bp::object(result_vec));
  }),
        py::arg("lset"), py::arg("mesh"),
        py::arg("cf_neg")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("cf_pos")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("cf_interface")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("order")=5, py::arg("subdivlvl")=0, py::arg("domains")=py::dict(), py::arg("heapsize")=1000000)
  ;


  typedef PyWrapper<BilinearFormIntegrator> PyBFI;
  typedef PyWrapper<LinearFormIntegrator> PyLFI;

  m.def("SymbolicCutBFI", FunctionPointer
          ([](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              py::object definedon)
          -> PyBFI
  {

    py::extract<Region> defon_region(definedon);
    if (defon_region.check())
      vb = VorB(defon_region());

    if (vb == BND)
      throw Exception("Symbolic cuts not yet (tested) for boundaries..");

    // check for DG terms
    bool has_other = false;
    cf.Get()->TraverseTree ([&has_other] (CoefficientFunction & cf)
    {
      if (dynamic_cast<ProxyFunction*> (&cf))
        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
          has_other = true;
    });
    if (has_other && !element_boundary && !skeleton)
      throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");
    if (element_boundary)
      throw Exception("No Facet BFI with Symbolic cuts..");

    shared_ptr<BilinearFormIntegrator> bfi;
    if (!has_other && !skeleton)
      bfi = make_shared<SymbolicCutBilinearFormIntegrator> (lset.Get(), cf.Get(), dt, order, subdivlvl);
    else
      bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (lset.Get(), cf.Get(), dt, order, subdivlvl);

    if (py::extract<py::list> (definedon).check())
      bfi -> SetDefinedOn (makeCArray<int> (definedon));

    if (defon_region.check())
    {
      cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
      bfi->SetDefinedOn(defon_region().Mask());
    }

    return PyBFI(bfi);
  }),
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("subdivlvl")=0,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=DummyArgument()
        );


  m.def("SymbolicCutLFI", FunctionPointer
          ([](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              py::object definedon)
          -> PyLFI
  {

    py::extract<Region> defon_region(definedon);
    if (defon_region.check())
      vb = VorB(defon_region());

    if (vb == BND)
      throw Exception("Symbolic cuts not yet (tested) for boundaries..");

    if (element_boundary || skeleton)
      throw Exception("No Facet LFI with Symbolic cuts..");

    shared_ptr<LinearFormIntegrator> lfi
      = make_shared<SymbolicCutLinearFormIntegrator> (lset.Get(), cf.Get(), dt, order, subdivlvl);

    if (py::extract<py::list> (definedon).check())
      lfi -> SetDefinedOn (makeCArray<int> (definedon));

    if (defon_region.check())
    {
      cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
      lfi->SetDefinedOn(defon_region().Mask());
    }

    return PyLFI(lfi);
  }),
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("subdivlvl")=0,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=py::bool_(false),
        py::arg("skeleton")=py::bool_(false),
        py::arg("definedon")=DummyArgument()
        );

  typedef PyWrapperDerived<ProxyFunction, CoefficientFunction> PyProxyFunction;
  m.def("dn", FunctionPointer
          ([] (const PyProxyFunction self, int order)
  {
    shared_ptr<DifferentialOperator> diffopdudnk;
    switch (order)
    {
    case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,1>>> (); break;
    case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,2>>> (); break;
    case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,3>>> (); break;
    case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    default : throw Exception("no order higher than 8 implemented yet");
    }

    auto adddiffop = make_shared<ProxyFunction> (self.Get()->IsTestFunction(), self.Get()->IsComplex(),
                                                 diffopdudnk, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self.Get()->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));
    return PyProxyFunction(adddiffop);
  }));

  m.def("dn", FunctionPointer
          ([](PyGF self, int order) -> PyCF
  {
    shared_ptr<DifferentialOperator> diffopdudnk;
    switch (order)
    {
    case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,1>>> (); break;
    case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,2>>> (); break;
    case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,3>>> (); break;
    case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
    case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,5>>> (); break;
    case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,6>>> (); break;
    case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,7>>> (); break;
    case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,8>>> (); break;
    default : throw Exception("no order higher than 8 implemented yet");
    }
    return PyCF(make_shared<GridFunctionCoefficientFunction> (self.Get(), diffopdudnk));
  }));

  // bp::def("GFCoeff", FunctionPointer( [] (shared_ptr<GridFunction> in) { return dynamic_pointer_cast<CoefficientFunction>(make_shared<GridFunctionCoefficientFunction>(in)); } ) );

  // bp::implicitly_convertible
  //   <shared_ptr<GridFunctionCoefficientFunction>,
  //   shared_ptr<CoefficientFunction> >();


  // void RefineAtLevelSet (PyGF gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh){


  // bp::docstring_options local_docstring_options(true, true, false);

  // std::string nested_name = "comp";
  // if( bp::scope() )
  //   nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");
  typedef PyWrapperDerived<SFESpace, FESpace> PySFES;
  // py::class_<PySFES, PyFES>
  //   (m, "SFESpace")
  //   .def("Something", FunctionPointer ([](PySFES self, PyCF cf)
  //                                      { throw Exception ("Something called"); }),
  //        "test something")
  //   .def("__init__", py::make_constructor
  //        (FunctionPointer ([](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
  //                          {
  //                            Flags flags = py::extract<Flags> (bpflags)();
  //                            shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset.Get(), order, flags);
  //                            LocalHeap lh (1000000, "SFESpace::Update-heap", true);
  //                            ret->Update(lh);
  //                            return ret;
  //                          })));


  m.def("SFESpace", FunctionPointer
          ([](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
          -> PyFES
  {
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset.Get(), order, flags);
    LocalHeap lh (1000000, "SFESpace::Update-heap", true);
    ret->Update(lh);
    ret->FinalizeUpdate(lh);
    return ret;
  }));
  // new implementation: only straight cuts - start with triangles only for a start!

  m.def("NewIntegrateX",
          FunctionPointer([](py::object lset,
                             shared_ptr<MeshAccess> ma, 
                             PyCF cf,
                             int order,
                             DOMAIN_TYPE dt,
                             int heapsize)
                          {
                            static Timer timer ("NewIntegrateX");
                            static Timer timercutgeom ("NewIntegrateX::MakeCutGeom");
                            static Timer timerevalintrule ("NewIntegrateX::EvaluateIntRule");
                            static Timer timeradding("NewIntegrateX::Adding");
                            static Timer timermapir("NewIntegrateX::MapingIntergrRule");
                            static Timer timergetdnums1("NewIntegrateX::GetDNums1");
                            static Timer timergetdnums2("NewIntegrateX::GetDNums2");

                            py::extract<PyGF> pygf(lset);
                            if (!pygf.check())
                              throw Exception("cast failed... need new candidates..");
                            shared_ptr<GridFunction> gf_lset = pygf().Get();

                            RegionTimer reg (timer);
                            LocalHeap lh(heapsize, "lh-New-Integrate");

                            double sum = 0.0;
                            int DIM = ma->GetDimension();

                            auto FESpace = gf_lset->GetFESpace();
                            cout << "I found a FESpace of order " << FESpace->GetOrder() << " in NewIntegrateX" << endl;

                            Array<int> dnums;
                            ma->IterateElements
                              (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
                               {
                                 //cout << endl << "NewIntegrateX() loop: Element Nr. " << el.Nr() << endl;
                                 auto & trafo = ma->GetTrafo (el, lh);

                                 timergetdnums1.Start();
                                 FESpace->GetDofNrs(el.Nr(),dnums);
                                 FlatVector<> elvec(dnums.Size(),lh);
                                 timergetdnums1.Stop();
                                 timergetdnums2.Start();
                                 gf_lset->GetVector().GetIndirect(dnums,elvec);
                                 timergetdnums2.Stop();

                                 timercutgeom.Start();
                                 const IntegrationRule * ir;
                                 ir = StraightCutIntegrationRule(elvec, trafo, dt, order, lh);

                                 timercutgeom.Stop();
                                 timerevalintrule.Start();
                                 if (ir != nullptr)
                                 {
                                   timermapir.Start();
                                   BaseMappedIntegrationRule & mir = trafo(*ir, lh);
                                   FlatMatrix<> val(mir.Size(), 1, lh);
                                   timermapir.Stop();


                                   cf -> Evaluate (mir, val);

                                   timeradding.Start();
                                   double lsum = 0.0;
                                   for (int i = 0; i < mir.Size(); i++)
                                     lsum += mir[i].GetWeight()*val(i,0);

                                   AsAtomic(sum) += lsum;
                                   timeradding.Stop();
                                 }
                                 timerevalintrule.Stop();
                               });

                            return sum;
                          }),
          py::arg("lset"),
           py::arg("mesh"),
           py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
           py::arg("order")=5,
           py::arg("domain_type")=IF,
           py::arg("heapsize")=1000000);
}

PYBIND11_PLUGIN(libngsxfem_py)
{
  cout << "importing ngs-xfem-" << NGSXFEM_VERSION << endl;
  py::module m("xfem", "pybind xfem");
  ExportNgsx(m);
  return m.ptr();
}
