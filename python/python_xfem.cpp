//#include "../ngstd/python_ngstd.hpp"
#include <regex>
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
#include "../utils/restrictedblf.hpp"
#include "../cutint/spacetimecutrule.hpp"
#include "../lsetcurving/shiftedevaluate.hpp"

#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"
#include "../spacetime/diffopDt.hpp"
#include "../spacetime/timecf.hpp"
// #include "../utils/error.hpp"

//using namespace ngcomp;

void ExportNgsx(py::module &m)
{




  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  typedef shared_ptr<CutInformation> PyCI;
  typedef shared_ptr<BitArray> PyBA;

  py::enum_<DOMAIN_TYPE>(m, "DOMAIN_TYPE")
  .value("POS", POS)
  .value("NEG", NEG)
  .value("IF", IF)
  .export_values()
  ;

  // typedef PyWrapperDerived<CompoundFESpace, FESpace> PyCompFES;

  typedef shared_ptr<XFESpace> PyXFES;

  m.def("XToNegPos",  [] (PyGF gfx, PyGF gfnegpos) {
    XFESpace::XToNegPos(gfx,gfnegpos);
  }  );

  py::class_<CutInformation, shared_ptr<CutInformation>>
    (m, "CutInfo")
  .def("__init__",  [] (CutInformation *instance,
                                        shared_ptr<MeshAccess> ma,
                                        py::object lset,
                                        int time_order,
                                        int heapsize)
  {
    new (instance) CutInformation (ma);
    if (py::extract<PyCF> (lset).check())
    {
      PyCF cflset = py::extract<PyCF>(lset)();
      LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
      instance->Update(cflset,time_order,lh);
    }
  },
       py::arg("mesh"),
       py::arg("levelset") = DummyArgument(),
       py::arg("time_order")=-1,
       py::arg("heapsize") = 1000000
       )
  .def("Update", [](CutInformation & self,
                                     PyCF lset,
                                     int time_order,
                                     int heapsize)
  {
    LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
    self.Update(lset,time_order,lh);
  },
       py::arg("levelset"),
       py::arg("time_order")=-1,
       py::arg("heapsize") = 1000000
       )
  .def("Mesh", [](CutInformation & self)
  {
    return self.GetMesh();
  }
       )
  .def("GetElementsOfType", [](CutInformation & self,
                                                DOMAIN_TYPE dt,
                                                VorB vb)
  {
    return self.GetElementsOfDomainType(dt,vb);
  },
       py::arg("domain_type") = IF,
       py::arg("VOL_or_BND") = VOL
       )
  .def("GetFacetsOfType", [](CutInformation & self,
                                              DOMAIN_TYPE dt)
  {
    return self.GetFacetsOfDomainType(dt);
  },
       py::arg("domain_type") = IF)

  .def("GetCutRatios", [](CutInformation & self,
                                           VorB vb)
  {
    return self.GetCutRatios(vb);
  },
       py::arg("VOL_or_BND") = VOL)
  ;


  m.def("GetFacetsWithNeighborTypes",
         [] (shared_ptr<MeshAccess> ma,
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
  } ,
        py::arg("mesh"),
        py::arg("a"),
        py::arg("bnd_val_a") = true,
        py::arg("bnd_val_b") = true,
        py::arg("use_and") = true,
        py::arg("b") = DummyArgument(),
        py::arg("heapsize") = 1000000
        );

  m.def("GetElementsWithNeighborFacets",
         [] (shared_ptr<MeshAccess> ma,
                             shared_ptr<BitArray> a,
                             int heapsize)
  {
    LocalHeap lh (heapsize, "GetElementsWithNeighborFacets-heap", true);
    return GetElementsWithNeighborFacets(ma,a,lh);
  } ,
        py::arg("mesh"),
        py::arg("a"),
        py::arg("heapsize") = 1000000
        );

  m.def("GetDofsOfElements",
         [] (PyFES fes,
                             PyBA a,
                             int heapsize)
  {
    LocalHeap lh (heapsize, "GetDofsOfElements-heap", true);
    return GetDofsOfElements(fes,a,lh);
  } ,
        py::arg("space"),
        py::arg("a"),
        py::arg("heapsize") = 1000000
        );

  m.def("CompoundBitArray",
         [] (py::list balist)
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
  } ,
        py::arg("balist")
        );



  //typedef shared_ptr<RestrictedBilinearForm> PyRBLF;
  //py::class_<RestrictedBilinearForm, PyRBLF, BilinearForm>
  //  (m, "CRestrictedBilinearForm");
  m.def("RestrictedBilinearForm",
         [](shared_ptr<FESpace> fes,
            const string & aname,
            py::object ael_restriction,
            py::object afac_restriction,
            bool check_unused,
            py::dict bpflags)
         {
           Flags flags = py::extract<Flags> (bpflags)();

           shared_ptr<BitArray> el_restriction = nullptr;
           shared_ptr<BitArray> fac_restriction = nullptr;
           if (py::extract<PyBA> (ael_restriction).check())
             el_restriction = py::extract<PyBA>(ael_restriction)();

           if (py::extract<PyBA> (afac_restriction).check())
             fac_restriction = py::extract<PyBA>(afac_restriction)();

           if (fes->IsComplex())
             throw Exception("RestrictedBilinearForm not implemented for complex fespace");
           
           shared_ptr<BilinearForm> biform = make_shared<RestrictedBilinearForm> (fes, aname, el_restriction, fac_restriction, flags);
           biform -> SetCheckUnused (check_unused);                             
           return biform;
         },
         py::arg("space"),
         py::arg("name") = "bfa",
         py::arg("element_restriction") = DummyArgument(),
         py::arg("facet_restriction") = DummyArgument(),
         py::arg("check_unused") = true,
         py::arg("flags") = py::dict()
      );

  
  typedef shared_ptr<BitArrayCoefficientFunction> PyBACF;
  py::class_<BitArrayCoefficientFunction, PyBACF, CoefficientFunction>
    (m, "BitArrayCF")
  .def("__init__",
       [](BitArrayCoefficientFunction *instance, shared_ptr<BitArray> ba)
  {
    new (instance) BitArrayCoefficientFunction (ba);
  },
       py::arg("bitarray")
       );

//   .def("__init__",  [] (XFESpace *instance,
  m.def("XFESpace", [] (
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
      cf_lset = py::extract<PyCF>(acutinfo)();
    if (py::extract<PyCF> (alset).check())
      cf_lset = py::extract<PyCF>(alset)();


    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<XFESpace> ret = nullptr;
    shared_ptr<MeshAccess> ma = basefes->GetMeshAccess();
    if (cutinfo)
    {
      if (ma->GetDimension()==2)
        ret = make_shared<T_XFESpace<2>> (ma, basefes, cutinfo, flags);
      else
        ret = make_shared<T_XFESpace<3>> (ma, basefes, cutinfo, flags);
    }
    else if (cf_lset)
    {
      if (ma->GetDimension()==2)
        ret = make_shared<T_XFESpace<2>> (ma, basefes, cf_lset, flags);
      else
        ret = make_shared<T_XFESpace<3>> (ma, basefes, cf_lset, flags);
    }
    else
      throw Exception("levelset and cutinfo are invalid");
    LocalHeap lh (heapsize, "XFESpace::Update-heap", true);
    ret->Update(lh);
    return ret;
  },
       py::arg("basefes"),
       py::arg("cutinfo") = DummyArgument(),
       py::arg("lset") = DummyArgument(),
       py::arg("flags") = py::dict(),
       py::arg("heapsize") = 1000000)
  ;
         

  py::class_<XFESpace, PyXFES, FESpace>
    (m, "CXFESpace")
  .def("GetCutInfo", [](PyXFES self)
  {
    return self->GetCutInfo();
  },
       "Get Information of cut geometry")
  .def("BaseDofOfXDof", [](PyXFES self, int i)
  {
    return self->GetBaseDofOfXDof(i);
  },
       "get corresponding dof of base FESpace")
  .def("GetDomainOfDof", [](PyXFES self, int i)
  {
    return self->GetDomainOfDof(i);
  },
       "get domain_type of degree of freedom")
  .def("GetDomainNrs",   [] (PyXFES self, int elnr) {
    Array<DOMAIN_TYPE> domnums;
    self->GetDomainNrs( elnr, domnums );
    return domnums;
  })
  ;

  m.def("InterpolateToP1",  [] (PyGF gf_ho, PyGF gf_p1, double eps_perturbation, int heapsize)
  {
    InterpolateP1 interpol(gf_ho, gf_p1);
    LocalHeap lh (heapsize, "InterpolateP1-Heap");
    interpol.Do(lh,eps_perturbation);
  } ,
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,
        py::arg("eps_perturbation")=1e-16,py::arg("heapsize")=1000000)
  ;

  m.def("InterpolateToP1",  [] (PyCF coef, PyGF gf_p1, double eps_perturbation, int heapsize)
  {
    InterpolateP1 interpol(coef, gf_p1);
    LocalHeap lh (heapsize, "InterpolateP1-Heap");
    interpol.Do(lh,eps_perturbation);
  } ,
        py::arg("coef"),py::arg("gf"),
        py::arg("eps_perturbation")=1e-16,py::arg("heapsize")=1000000)
  ;

  py::class_<StatisticContainer, shared_ptr<StatisticContainer>>(m, "StatisticContainer")
  .def(py::init<>())
  .def("Print", [](StatisticContainer & self, string label, string select)
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
  },
       py::arg("label")="something",py::arg("select")="all"
       )
  ;

  m.def("CalcMaxDistance",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, int heapsize)
  {
    StatisticContainer dummy;
    LocalHeap lh (heapsize, "CalcDistance-Heap");
    if (lset_p1->GetMeshAccess()->GetDimension()==2)
      CalcDistances<2>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
    else
      CalcDistances<3>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
    return (double) dummy.ErrorMaxNorm[dummy.ErrorMaxNorm.Size()-1];
  } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("heapsize")=1000000)
  ;

  m.def("CalcDistances",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
  {
    LocalHeap lh (heapsize, "CalcDistance-Heap");
    if (lset_p1->GetMeshAccess()->GetDimension()==2)
      CalcDistances<2>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
    else
      CalcDistances<3>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
  } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("stats")=NULL,py::arg("heapsize")=1000000,py::arg("refine_threshold")=-1.0,py::arg("absolute")=false)
  ;

  m.def("CalcDeformationError",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, StatisticContainer & stats, double lower, double upper, int heapsize)
  {
    LocalHeap lh (heapsize, "CalcDeformationError-Heap");
    if (lset_p1->GetMeshAccess()->GetDimension()==2)
      CalcDeformationError<2>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
    else
      CalcDeformationError<3>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
  } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL,py::arg("stats")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
  ;

  m.def("ProjectShift", [] (PyGF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, py::object pba, double lower, double upper, double threshold, int heapsize)
  {
    shared_ptr<BitArray> ba = nullptr;
    if (py::extract<PyBA> (pba).check())
      ba = py::extract<PyBA>(pba)();
    LocalHeap lh (heapsize, "ProjectShift-Heap");
    ProjectShift(lset_ho, lset_p1, deform, qn, ba, lower, upper, threshold, lh);
  } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL, py::arg("ba")=py::none(),py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("threshold")=1.0,py::arg("heapsize")=1000000)
  ;

// ProjectShift


  m.def("RefineAtLevelSet",  [] (PyGF lset_p1, double lower, double upper, int heapsize)
  {
    LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
    RefineAtLevelSet(lset_p1, lower, upper, lh);
  } ,
        py::arg("lset_p1")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
  ;



  m.def("RefineAtLevelSet",  [] (PyGF gf, double lower_lset_bound, double upper_lset_bound, int heapsize)
  {
    LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
    RefineAtLevelSet(gf,lower_lset_bound,upper_lset_bound,lh);
  } ,
        py::arg("gf"),py::arg("lower_lset_bound")=0.0,py::arg("upper_lset_bound")=0.0,py::arg("heapsize")=10000000)
  ;

  typedef shared_ptr<BilinearFormIntegrator> PyBFI;
  typedef shared_ptr<LinearFormIntegrator> PyLFI;

  m.def("SymbolicCutBFI", [](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int time_order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              py::object definedon,
              py::object definedonelem)
          -> PyBFI
  {

    py::extract<Region> defon_region(definedon);
    if (defon_region.check())
      vb = VorB(defon_region());

    if (vb == BND)
      throw Exception("Symbolic cuts not yet (tested) for boundaries..");
    // check for DG terms
    bool has_other = false;
    cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
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
    {
      auto bfime = make_shared<SymbolicCutBilinearFormIntegrator> (lset, cf, dt, order, subdivlvl,vb);
      bfime->SetTimeIntegrationOrder(time_order);
      bfi = bfime;
    }
    else
    {
      if (vb == BND)
        throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for boundaries..");
      
      bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (lset, cf, dt, order, subdivlvl);
    }
    if (py::extract<py::list> (definedon).check())
      bfi -> SetDefinedOn (makeCArray<int> (definedon));

    if (defon_region.check())
    {
      cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
      bfi->SetDefinedOn(defon_region().Mask());
    }

    if (! py::extract<DummyArgument> (definedonelem).check())
      bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

    return PyBFI(bfi);
  },
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("subdivlvl")=0,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=DummyArgument(),
        py::arg("definedonelements")=DummyArgument()
        );

  m.def("SymbolicFacetPatchBFI", [](PyCF cf,
                                    int order,
                                    int time_order,
                                    bool skeleton,
                                    py::object definedonelem)
        -> PyBFI
  {
    // check for DG terms
    bool has_other = false;
    cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
    {
      if (dynamic_cast<ProxyFunction*> (&cf))
        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
          has_other = true;
    });
    if (!has_other)
      cout << " no Other() used?!" << endl;

  
    shared_ptr<BilinearFormIntegrator> bfi;
    if (skeleton)
    {
      auto bfime = make_shared<SymbolicFacetBilinearFormIntegrator2> (cf, order);
      bfime->SetTimeIntegrationOrder(time_order);
      bfi = bfime;
    }
    else
    {
      // throw Exception("Patch facet blf not implemented yet: TODO(2)!");
      auto bfime = make_shared<SymbolicFacetPatchBilinearFormIntegrator> (cf, order);
      bfime->SetTimeIntegrationOrder(time_order);
      bfi = bfime;
    }

    if (! py::extract<DummyArgument> (definedonelem).check())
      bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

    return PyBFI(bfi);
  },
        py::arg("form"),
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("skeleton") = true,
        py::arg("definedonelements")=DummyArgument()
    );

  m.def("SymbolicCutLFI", [](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int time_order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              py::object definedon,
              py::object definedonelem)
          -> PyLFI
  {

    py::extract<Region> defon_region(definedon);
    if (defon_region.check())
      vb = VorB(defon_region());

    // if (vb == BND)
    //   throw Exception("Symbolic cuts not yet (tested) for boundaries..");
    if (vb == BND)
      throw Exception("Symbolic cuts not yet (tested) for boundaries..");

    if (element_boundary || skeleton)
      throw Exception("No Facet LFI with Symbolic cuts..");

    auto lfime  = make_shared<SymbolicCutLinearFormIntegrator> (lset, cf, dt, order, subdivlvl,vb);
    lfime->SetTimeIntegrationOrder(time_order);
    shared_ptr<LinearFormIntegrator> lfi = lfime;

    if (py::extract<py::list> (definedon).check())
      lfi -> SetDefinedOn (makeCArray<int> (definedon));

    if (defon_region.check())
    {
      cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
      lfi->SetDefinedOn(defon_region().Mask());
    }

    if (! py::extract<DummyArgument> (definedonelem).check())
      lfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

    return PyLFI(lfi);
  },
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("subdivlvl")=0,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=py::bool_(false),
        py::arg("skeleton")=py::bool_(false),
        py::arg("definedon")=DummyArgument(),
        py::arg("definedonelements")=DummyArgument()
        );

  typedef shared_ptr<ProxyFunction> PyProxyFunction;
  m.def("dn", [] (const PyProxyFunction self, int order, py::object comp, int dim_space, bool hdiv)
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

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopdudnk;
    if (! hdiv)
    {
      if (dim_space == 2)
      {
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
      }
      else
      {
      switch (order)
      {
      case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,1>>> (); break;
      case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,2>>> (); break;
      case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,3>>> (); break;
      case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,4>>> (); break;
      case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,5>>> (); break;
      case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,6>>> (); break;
      case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,7>>> (); break;
      case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,8>>> (); break;
      default : throw Exception("no order higher than 8 implemented yet");
      }
      }
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

    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopdudnk = make_shared<CompoundDifferentialOperator> (diffopdudnk, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->IsTestFunction(), self->IsComplex(),
                                                 diffopdudnk, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
  },
        py::arg("proxy"),
        py::arg("order"),
        py::arg("comp") = -1,
        py::arg("dim") = 2,
        py::arg("hdiv") = false
        );

  m.def("dn", [](PyGF self, int order) -> PyCF
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
    return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopdudnk));
  });

  // bp::def("GFCoeff", FunctionPointer( [] (shared_ptr<GridFunction> in) { return dynamic_pointer_cast<CoefficientFunction>(make_shared<GridFunctionCoefficientFunction>(in)); } ) );

  // bp::implicitly_convertible
  //   <shared_ptr<GridFunctionCoefficientFunction>,
  //   shared_ptr<CoefficientFunction> >();


  // void RefineAtLevelSet (PyGF gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh){


  // bp::docstring_options local_docstring_options(true, true, false);

  // std::string nested_name = "comp";
  // if( bp::scope() )
  //   nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");

  typedef shared_ptr<SFESpace> PySFES;

  m.def("SFESpace", [](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
          -> PyFES
  {
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset, order, flags);
    LocalHeap lh (1000000, "SFESpace::Update-heap", true);
    ret->Update(lh);
    ret->FinalizeUpdate(lh);
    return ret;
  });
  // new implementation: only straight cuts - start with triangles only for a start!

  m.def("DebugSpaceTimeCutIntegrationRule", [](){ DebugSpaceTimeCutIntegrationRule(); });
  m.def("IntegrateX",
        [](py::object lset,
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
    tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(pycf(),subdivlvl);

    LocalHeap lh(heapsize, "lh-IntegrateX");

    double sum = 0.0;
    int DIM = ma->GetDimension();

    Array<int> dnums;
    ma->IterateElements
      (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
    {
      auto & trafo = ma->GetTrafo (el, lh);

      const IntegrationRule * ir = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, order,-1, lh, subdivlvl);

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
  },
        py::arg("lset"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("order")=5,
        py::arg("domain_type")=IF,
        py::arg("subdivlvl")=0,
        py::arg("heapsize")=1000000);
  typedef shared_ptr<SpaceTimeFESpace> PySTFES;


  m.def("SpaceTimeFESpace", [] (
                                        PyFES basefes,
                                        shared_ptr<FiniteElement> fe,
                                        py::object dirichlet,
                                        py::dict bpflags,
                                        int heapsize)
  {


    shared_ptr<SpaceTimeFESpace> ret = nullptr;
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<MeshAccess> ma = basefes->GetMeshAccess();

    if (py::isinstance<py::list>(dirichlet)) {
        flags.SetFlag("dirichlet", makeCArray<double>(py::list(dirichlet)));

    }

    if (py::isinstance<py::str>(dirichlet))
    {
          std::regex pattern(dirichlet.cast<string>());
          Array<double> dirlist;
          for (int i = 0; i < ma->GetNBoundaries(); i++)
             if (std::regex_match (ma->GetMaterial(BND, i), pattern))
               dirlist.Append (i+1);
               flags.SetFlag("dirichlet", dirlist);
    }

    auto tfe = dynamic_pointer_cast<ScalarFiniteElement<1>>(fe);
    cout << tfe << endl;
    if(tfe == nullptr)
      cout << "Warning!" << endl;

    ret = make_shared<SpaceTimeFESpace> (ma, basefes,tfe, flags);



    LocalHeap lh (heapsize, "SpaceTimeFESpace::Update-heap", true);
    ret->Update(lh);
    ret->FinalizeUpdate(lh);
    return ret;
  },
       py::arg("spacefes"),
       py::arg("timefe"),
       py::arg("dirichlet") = DummyArgument(),
       py::arg("flags") = py::dict(),
        py::arg("heapsize") = 1000000);

  py::class_<SpaceTimeFESpace, PySTFES, FESpace>
    (m, "CSpaceTimeFESpace")
  .def("SetTime", [](PySTFES self, double t)
  {
    self->SetTime(t);
  },
       "Set the time variable")
  .def("SetOverrideTime", [](PySTFES self, bool override)
  {
    self->SetOverrideTime(override);
  },
       "Set flag to or not to override the time variable")
  .def("k_t", [](PySTFES self)
  {
     return self->order_time();
  },
     "Return order of the time FE")
  .def("TimeFE_nodes", [](PySTFES self)
  {
      Vector<double> intp_pts(self->order_time() + 1);
      self->TimeFE_nodes(intp_pts);
      return intp_pts;
   },
     "Return nodes of the time FE")
  ;

  m.def("ScalarTimeFE", []( int order)
  {
    BaseScalarFiniteElement * fe = nullptr;

    fe = new NodalTimeFE(order);


    return shared_ptr<BaseScalarFiniteElement>(fe);
  },
        "creates nodal FE in time based on Gauss-Lobatto integration points"
        )
   ;


  // DiffOpDt

  m.def("dt", [] (const PyProxyFunction self,py::object comp)
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

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopdt;
    diffopdt = make_shared<T_DifferentialOperator<DiffOpDt>> ();

    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopdt = make_shared<CompoundDifferentialOperator> (diffopdt, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->IsTestFunction(), self->IsComplex(),
                                                 diffopdt, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
    },
          py::arg("proxy"),
          py::arg("comp") = -1
          );

  m.def("dt", [](PyGF self) -> PyCF
  {
    shared_ptr<DifferentialOperator> diffopdt;
    diffopdt = make_shared<T_DifferentialOperator<DiffOpDt>> ();

    return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopdt));
  });

   m.def("ReferenceTimeVariable", []() -> PyCF
   {
     return PyCF(make_shared<TimeVariableCoefficientFunction> ());
   });


   // DiffOpFixt

  m.def("fix_t", [] (const PyProxyFunction self, double time, py::object comp, bool use_FixAnyTime )
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

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopfixt;

    if(!use_FixAnyTime && (time == 0.0 || time == 1.0))
    {
      switch (int(time))
      {
        case 0 : diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<0>>> (); break;
        case 1 : diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<1>>> (); break;
        default : throw Exception("Requested time not implemented yet.");
      }
    }
    else {
      cout << "Calling DiffOpFixAnyTime" << endl;
      diffopfixt = make_shared<DiffOpFixAnyTime> (time);
    }


    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopfixt = make_shared<CompoundDifferentialOperator> (diffopfixt, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->IsTestFunction(), self->IsComplex(),
                                                 diffopfixt, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
    },
          py::arg("proxy"),
          py::arg("time"),
          py::arg("comp") = -1,
          py::arg("use_FixAnyTime") = false
          );

   m.def("fix_t", [](PyGF self, double time, bool use_FixAnyTime) -> PyCF
   {
     shared_ptr<DifferentialOperator> diffopfixt;

     if(!use_FixAnyTime && (time == 0.0 || time == 1.0))
     {
       switch (int(time))
       {
         case 0 : diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<0>>> (); break;
         case 1 : diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<1>>> (); break;
         default : throw Exception("Requested time not implemented yet.");
       }
     }
     else {
       cout << "Calling DiffOpFixAnyTime" << endl;
       diffopfixt = make_shared<DiffOpFixAnyTime> (time);
     }


     return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopfixt));
   });

   m.def("shifted_eval", [](PyGF self, PyGF back,PyGF forth ) -> PyCF
   {

     auto diffop = make_shared<DiffOpShiftedEval> (back,forth);


     return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffop));
   });


   m.def("ReferenceTimeVariable", []() -> PyCF
   {
     return PyCF(make_shared<TimeVariableCoefficientFunction> ());
   });

   m.def("RestrictToTime", [](PyGF st_GF,double time) -> PyGF
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FE = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     return st_FE->CreateRestrictedGF(st_GF,time);
   });
}

PYBIND11_PLUGIN(ngsxfem_py)
{
  cout << "importing ngs-xfem-" << NGSXFEM_VERSION << endl;
  py::module m("xfem", "pybind xfem");
  ExportNgsx(m);
  return m.ptr();
}
