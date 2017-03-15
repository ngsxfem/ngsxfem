//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../utils/bitarraycf.hpp"
#include "../xfem/cutinfo.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/sFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/ghostpenalty.hpp"
#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../cutint/straightcutrule.hpp"
#include "../cutint/xintegration.hpp"
// #include "../utils/error.hpp"

//using namespace ngcomp;

void ExportNgsx(py::module &m)
{




  typedef PyWrapper<FESpace> PyFES;
  typedef PyWrapper<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef PyWrapper<GF> PyGF;
  typedef PyWrapper<shared_ptr<CutInformation>> PyCI;
  typedef PyWrapper<shared_ptr<BitArray>> PyBA;

  py::enum_<DOMAIN_TYPE>(m, "DOMAIN_TYPE")
  .value("POS", POS)
  .value("NEG", NEG)
  .value("IF", IF)
  .export_values()
  ;

  // typedef PyWrapperDerived<CompoundFESpace, FESpace> PyCompFES;

  typedef PyWrapperDerived<XFESpace, FESpace> PyXFES;

  m.def("XToNegPos", FunctionPointer( [] (PyGF gfx, PyGF gfnegpos) {
    XFESpace::XToNegPos(gfx.Get(),gfnegpos.Get());
  } ) );

  py::class_<CutInformation, shared_ptr<CutInformation>>
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

  m.def("GetDofsOfElements",
        FunctionPointer( [] (PyFES fes,
                             PyBA a,
                             int heapsize)
  {
    LocalHeap lh (heapsize, "GetDofsOfElements-heap", true);
    return GetDofsOfElements(fes.Get(),a,lh);
  } ),
        py::arg("space"),
        py::arg("a"),
        py::arg("heapsize") = 1000000
        );

  m.def("CompoundBitArray",
        FunctionPointer( [] (py::list balist)
  {
    size_t cnt = 0;
    for( auto aba : balist )
    {
      shared_ptr<BitArray> ba = py::extract<PyBA>(aba)();
      cnt += ba->Size();
    }
    shared_ptr<BitArray> res = make_shared<BitArray>(cnt);
    res->Clear();
    size_t offset = 0;
    for( auto aba : balist )
    {
      shared_ptr<BitArray> ba = py::extract<PyBA>(aba)();
      for (size_t i = 0; i < ba->Size(); ++i)
      {
        if (ba->Test(i))
          res->Set(offset+i);
      }
      offset += ba->Size();
    }
    return res;
  } ),
        py::arg("balist")
        );



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
                                        py::object acutinfo,
                                        py::object alset,
                                        py::dict bpflags,
                                        int heapsize)
  {
    shared_ptr<CoefficientFunction> cf_lset = nullptr;
    shared_ptr<CutInformation> cutinfo = nullptr;
    if (py::extract<PyCI> (acutinfo).check())
      cutinfo = py::extract<PyCI>(acutinfo)();
    if (py::extract<PyCF> (acutinfo).check())
      cf_lset = py::extract<PyCF>(acutinfo)().Get();
    if (py::extract<PyCF> (alset).check())
      cf_lset = py::extract<PyCF>(alset)().Get();


    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<XFESpace> ret = nullptr;
    shared_ptr<MeshAccess> ma = basefes.Get()->GetMeshAccess();
    if (cutinfo)
    {
      if (ma->GetDimension()==2)
        ret = make_shared<T_XFESpace<2>> (ma, basefes.Get(), cutinfo, flags);
      else
        ret = make_shared<T_XFESpace<3>> (ma, basefes.Get(), cutinfo, flags);
    }
    else if (cf_lset)
    {
      if (ma->GetDimension()==2)
        ret = make_shared<T_XFESpace<2>> (ma, basefes.Get(), cf_lset, flags);
      else
        ret = make_shared<T_XFESpace<3>> (ma, basefes.Get(), cf_lset, flags);
    }
    else
      throw Exception("levelset and cutinfo are invalid");
    LocalHeap lh (heapsize, "XFESpace::Update-heap", true);
    ret->Update(lh);
    new (instance) PyFES(ret);
  }),
       py::arg("basefes"),
       py::arg("cutinfo") = DummyArgument(),
       py::arg("lset") = DummyArgument(),
       py::arg("flags") = py::dict(),
       py::arg("heapsize") = 1000000)
  .def("GetCutInfo", FunctionPointer ([](PyXFES self)
  {
    return self.Get()->GetCutInfo();
  }),
       "Get Information of cut geometry")
  .def("BaseDofOfXDof", FunctionPointer ([](PyXFES self, int i)
  {
    return self.Get()->GetBaseDofOfXDof(i);
  }),
       "get corresponding dof of base FESpace")
  .def("GetDomainOfDof", FunctionPointer ([](PyXFES self, int i)
  {
    return self.Get()->GetDomainOfDof(i);
  }),
       "get domain_type of degree of freedom")
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
          ([] (const PyProxyFunction self, int order, py::object comp, bool hdiv)
  {

    Array<int> comparr(0);
    if (py::extract<int> (comp).check())
    {
      int c = py::extract<int>(comp)();
      if (c != -1)
      {
        comparr.SetSize(1);
        comparr[0] = c;
      }
    }

    if (py::extract<py::list> (comp).check())
      comparr = makeCArray<int> (py::extract<py::list> (comp)());

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self.Get()->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopdudnk;
    if (! hdiv)
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
    else
      switch (order)
      {
      case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,1>>> (); break;
      case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,2>>> (); break;
      case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,3>>> (); break;
      case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,4>>> (); break;
      case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,5>>> (); break;
      case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,6>>> (); break;
      case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,7>>> (); break;
      case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,8>>> (); break;
      default : throw Exception("no order higher than 8 implemented yet");
      }

    for (int i = 0; i < comparr.Size(); ++i)
    {
      diffopdudnk = make_shared<CompoundDifferentialOperator> (diffopdudnk, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self.Get()->IsTestFunction(), self.Get()->IsComplex(),
                                                 diffopdudnk, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self.Get()->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
  }),
        py::arg("proxy"),
        py::arg("order"),
        py::arg("comp") = -1,
        py::arg("hdiv") = false
        );

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

  m.def("IntegrateX",
        FunctionPointer([](py::object lset,
                           shared_ptr<MeshAccess> ma,
                           PyCF cf,
                           int order,
                           DOMAIN_TYPE dt,
                           int subdivlvl,
                           int heapsize)
  {
    py::extract<PyCF> pycf(lset);
    if (!pycf.check())
      throw Exception("cast failed... need new candidates..");

    shared_ptr<GridFunction> gf_lset = nullptr;
    shared_ptr<CoefficientFunction> cf_lset = nullptr;
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(pycf().Get(),subdivlvl);

    LocalHeap lh(heapsize, "lh-IntegrateX");

    double sum = 0.0;
    int DIM = ma->GetDimension();

    Array<int> dnums;
    ma->IterateElements
      (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
    {
      auto & trafo = ma->GetTrafo (el, lh);

      const IntegrationRule * ir = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, order, lh, subdivlvl);

      if (ir != nullptr)
      {
        BaseMappedIntegrationRule & mir = trafo(*ir, lh);
        FlatMatrix<> val(mir.Size(), 1, lh);

        cf -> Evaluate (mir, val);

        double lsum = 0.0;
        for (int i = 0; i < mir.Size(); i++)
          lsum += mir[i].GetWeight()*val(i,0);

        AsAtomic(sum) += lsum;
      }
    });

    return sum;
  }),
        py::arg("lset"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("order")=5,
        py::arg("domain_type")=IF,
        py::arg("subdivlvl")=0,
        py::arg("heapsize")=1000000);
}

PYBIND11_PLUGIN(libngsxfem_py)
{
  cout << "importing ngs-xfem-" << NGSXFEM_VERSION << endl;
  py::module m("xfem", "pybind xfem");
  ExportNgsx(m);
  return m.ptr();
}
