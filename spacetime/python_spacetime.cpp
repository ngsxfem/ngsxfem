//#include "../ngstd/python_ngstd.hpp"
#include <regex>
#include <python_ngstd.hpp>
// #include "../utils/bitarraycf.hpp"
// #include "../xfem/cutinfo.hpp"
// #include "../xfem/xFESpace.hpp"
// #include "../xfem/sFESpace.hpp"
// #include "../xfem/symboliccutbfi.hpp"
// #include "../xfem/symboliccutlfi.hpp"
// #include "../xfem/ghostpenalty.hpp"
// #include "../lsetcurving/p1interpol.hpp"
// #include "../lsetcurving/calcgeomerrors.hpp"
// #include "../lsetcurving/lsetrefine.hpp"
// #include "../lsetcurving/projshift.hpp"
// #include "../cutint/straightcutrule.hpp"
// #include "../cutint/xintegration.hpp"
// #include "../utils/restrictedblf.hpp"
#include "../cutint/spacetimecutrule.hpp"
// #include "../lsetcurving/shiftedevaluate.hpp"
// #include "../utils/error.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"
#include "../spacetime/diffopDt.hpp"
#include "../spacetime/timecf.hpp"
#include "../spacetime/spacetime_vtk.hpp"

using namespace ngcomp;
using namespace xintegration;

void ExportNgsx_spacetime(py::module &m)
{




  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  typedef shared_ptr<BitArray> PyBA;
  typedef shared_ptr<SpaceTimeFESpace> PySTFES;
  typedef shared_ptr<ProxyFunction> PyProxyFunction;

  m.def("SpaceTimeFESpace", [] (
                                        PyFES basefes,
                                        shared_ptr<FiniteElement> fe,
                                        py::object dirichlet,
                                        int heapsize, 
                                        py::kwargs kwargs)
  {


    shared_ptr<SpaceTimeFESpace> ret = nullptr;
    //Flags flags = py::extract<Flags> (bpflags)();
    auto flags = CreateFlagsFromKwArgs(kwargs);


    shared_ptr<MeshAccess> ma = basefes->GetMeshAccess();

    if (py::isinstance<py::list>(dirichlet)) {
        flags.SetFlag("dirichlet", makeCArray<double>(py::list(dirichlet)));

    }

    if (py::isinstance<py::str>(dirichlet))
    {
        Array<double> dirlist;
        Region dir(ma, BND, dirichlet.cast<string>());
        for (int i = 0; i < ma->GetNBoundaries(); i++)
            if (dir.Mask()[i])
              dirlist.Append (i+1);
        flags.SetFlag("dirichlet", dirlist);
    }

    auto tfe = dynamic_pointer_cast<ScalarFiniteElement<1>>(fe);
    //cout << tfe << endl;
    if(tfe == nullptr)
      cout << IM(1) << "Warning! tfe == nullptr" << endl;

    ret = make_shared<SpaceTimeFESpace> (ma, basefes,tfe, flags);



    LocalHeap lh (heapsize, "SpaceTimeFESpace::Update-heap", true);
    ret->Update();
    ret->FinalizeUpdate();
    return ret;
  },
       py::arg("spacefes"),
       py::arg("timefe"),
       py::arg("dirichlet")=py::none(),
       py::arg("heapsize") = 1000000,
       docu_string(R"raw_string(
This function creates a SpaceTimeFiniteElementSpace based on a spacial FE space and a time Finite element
Roughly, this is the tensor product between those two arguments. Further arguments specify several details.

Parameters

spacefes : ngsolve.FESpace
  This is the spacial finite element used for the space-time discretisation.
  Both scalar and vector valued spaces might be used. An example would be
  spacefes = H1(mesh, order=order) for given mesh and order.

timefe : ngsolve.FiniteElement
  This is the time finite element for the space-time discretisation. That is
  essentially a simple finite element on the unit interval. There is a class
  ScalarTimeFE to create something fitting here. For example, one could call
  timefe = ScalarTimeFE(order) to create a time finite element of order order.

dirichlet : list or string
  The boundary of the space domain which should have Dirichlet boundary values.
  Specification policy is the same as with the usual space finite element spaces.

heapsize : int
  Size of the local heap of this class. Increase this if you observe errors which look
  like a heap overflow.

dgjumps : bool  
  )raw_string")
   );

  py::class_<SpaceTimeFESpace, PySTFES, FESpace>
    (m, "CSpaceTimeFESpace")
  .def("SetTime", [](PySTFES self, double t)
  {
    self->SetTime(t);
  },
       "Set the time variable\n Also sets override time")
  .def_property_readonly("spaceFES",
                [](PySTFES self) { return self->GetSpaceFESpace(); },
               "get space FESpace")
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
      Array<double> & nodesr = self->TimeFE_nodes();
      py::list nodes (nodesr.Size());
      for (int i = 0; i < nodesr.Size(); i++)
        nodes[i] = nodesr[i];
      return nodes;
   },
     "Return nodes of time FE")
  .def("IsTimeNodeActive", [](PySTFES self, int i)
  {
      return self->IsTimeNodeActive(i);
   },
     "Return bool whether node is active")
  ;


  py::class_<NodalTimeFE, shared_ptr<NodalTimeFE>,FiniteElement>(m, "ScalarTimeFE")
  .def(py::init([] ( int order, bool skip_first_nodes, bool only_first_nodes, bool skip_first_node, bool only_first_node) -> shared_ptr<NodalTimeFE>
  {
    if ((skip_first_nodes||skip_first_node) && (only_first_nodes||only_first_node))
      throw Exception("can't skip and keep first node at the same time.");
    return make_shared<NodalTimeFE>(order, skip_first_nodes || skip_first_node, only_first_nodes || only_first_node);
  }),
  py::arg("order") = 0,
  py::arg("skip_first_nodes") = false,
  py::arg("only_first_nodes") = false,
  py::arg("skip_first_node") = false,
  py::arg("only_first_node") = false,
  docu_string(R"raw_string(
Creates a nodal Finite element in time on the interval [0,1].
Internally, Gauss-Lobatto integration points are exploited for that.

Parameters

order : int
The polynomial order of the discretisation. That controlls the number of
points in the time interval. See Gauss-Lobatto points for further details.
Orders up to 5 are given by explicit closed formulas, beyond that an
iterative construction is applied.

skip_first_nodes : bool
This will create the time finite element without the first node at t=0.
That feature comes in handy for several CG like implementations in time.
Also see only_first_node.

only_first_nodes : bool
This will create the time finite element with only the first node at t=0.
That feature comes in handy for several CG like implementations in time.
Also see skip_first_node.
  )raw_string")
   )
  // .def("__rmul__", [](shared_ptr<NodalTimeFE> self, shared_ptr<FESpace> fes)
  // {
  //   Flags flags(fes->GetFlags());
  //   auto ret = make_shared<SpaceTimeFESpace> (fes->GetMeshAccess(), fes, self, fes->GetFlags());

  //   LocalHeap lh (1000000, "SpaceTimeFESpace::Update-heap", true);
  //   ret->Update();
  //   ret->FinalizeUpdate();
  //   return ret;
  // })
  .def("__mul__", [](shared_ptr<NodalTimeFE> self, shared_ptr<FESpace> fes)
  {
    Flags flags(fes->GetFlags());
    auto ret = make_shared<SpaceTimeFESpace> (fes->GetMeshAccess(), fes, self, fes->GetFlags());

    LocalHeap lh (1000000, "SpaceTimeFESpace::Update-heap", true);
    ret->Update();
    ret->FinalizeUpdate();
    return ret;
  })
  ;

  py::class_<GCC3FE, shared_ptr<GCC3FE>,NodalTimeFE>(m, "GCC3FE")
  .def(py::init([] ( bool skip_first_nodes, bool only_first_nodes) -> shared_ptr<GCC3FE>
  {
    if ((skip_first_nodes) && (only_first_nodes))
      throw Exception("can't skip and keep first node at the same time.");
    return make_shared<GCC3FE>(skip_first_nodes, only_first_nodes);
  }),
  py::arg("skip_first_nodes") = false,
  py::arg("only_first_nodes") = false,
  docu_string(R"raw_string(
docu missing
  )raw_string")
  );
  


  typedef shared_ptr<TimeVariableCoefficientFunction> PyTimeVariableCF;

  py::class_<TimeVariableCoefficientFunction, PyTimeVariableCF, CoefficientFunction>(m, "TimeVariableCoefficientFunction")
          .def("__init__", [] () -> PyTimeVariableCF { return make_shared<TimeVariableCoefficientFunction>(); })
          .def("FixTime", &TimeVariableCoefficientFunction::FixTime)
          .def("UnfixTime", &TimeVariableCoefficientFunction::UnfixTime)
          .def("IsFixed", [] (PyTimeVariableCF self) -> bool {
            try {
                self->Evaluate(BaseMappedIntegrationPoint());
            } catch (...) {
                 return false;
            }
            return true;
          });

   m.def("ReferenceTimeVariable", []() -> PyTimeVariableCF
   {
     return make_shared<TimeVariableCoefficientFunction> ();
   }, docu_string(R"raw_string(
This is the time variable. Call tref = ReferenceTimeVariable() to have a symbolic variable
for the time like x,y,z for space. That can be used e.g. in lset functions for unfitted methods.
Note that one would typically use tref in [0,1] as one time slab, leading to a call like
t = told + delta_t * tref, when tref is our ReferenceTimeVariable.
ngsxfem.__init__ defines tref.
)raw_string")
);


   // DiffOpFixt

  m.def("fix_tref_proxy", [] (const PyProxyFunction self, double time, py::object comp, bool use_FixAnyTime )
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

    auto scd= dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator());
    if (comparr.Size()== 0 && scd)
    {
      comparr.SetSize(1);
      comparr[0] = scd->Component();
      //throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopfixt;
    const int SpaceD = self->GetFESpace()->GetSpatialDimension();
    if(!use_FixAnyTime && (time == 0.0 || time == 1.0))
    {
      //if (SpaceD < 2)
        //throw Exception("SpaceD < 2 not implemented yet.");
      Switch<2> (int(time), [&] (auto TT) {
        if (SpaceD == 2)
        // Switch<2> (SpaceD-2, [&] (auto SD) {
          diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, TT>>> ();
        else if(SpaceD == 1)
            diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<1, TT>>> ();
        else
          diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, TT>>> ();
        // });
      });
    } 
    else {
      cout << IM(4) << "Calling DiffOpFixAnyTime" << endl;
      Switch<2> (SpaceD-2, [&] (auto SD) {
        diffopfixt = make_shared<DiffOpFixAnyTime<SD+2>> (time);
      });
    }


    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopfixt = make_shared<CompoundDifferentialOperator> (diffopfixt, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->GetFESpace(), self->IsTestFunction(), self->IsComplex(), diffopfixt, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
    },
          py::arg("proxy"),
          py::arg("time"),
          py::arg("comp") = -1,
          py::arg("use_FixAnyTime") = false
          );

  m.def("fix_tref_coef", [](PyCF self, py::object time) -> PyCF
  {
    shared_ptr<ParameterCoefficientFunction<double>> t = nullptr;
    auto t1 = py::extract<shared_ptr<ParameterCoefficientFunction<double>>> (time);
    if (t1.check())
    {
      t = t1();
    }
    else
    { 
      auto t2 = py::extract<double> (time);
      if (t2.check())
        t = make_shared<ParameterCoefficientFunction<double>>(t2());
      else  
        throw Exception("time object not valid");
    }
    return PyCF(make_shared<FixTimeCoefficientFunction> (self, t));
  },
  docu_string(R"raw_string(
fix_t fixes the evaluated time to a fixed value.

Parameters

self: ngsolve.CoefficientFunction
  CoefficientFunction in which the time should be fixed
  
time: Parameter or double
  Value the time should become (if Parameter, the value can be adjusted later on)

)raw_string")
     );
     
   m.def("fix_tref_gf", [](PyGF self, double time) -> PyCF
   {
     shared_ptr<DifferentialOperator> diffopfixt;
     const int SpaceD = self->GetFESpace()->GetSpatialDimension();
     if(time == 0.0 || time == 1.0)
     {
      if (SpaceD < 2) throw Exception("Requested time not implemented yet.");
      Switch<2> (int(time), [&] (auto TT) {
        if (SpaceD == 2)
        // Switch<2> (SpaceD-2, [&] (auto SD) {
          diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, TT>>> ();
        else
          diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, TT>>> ();
        // });
      });
     }
     else {
       cout << IM(4) << "Calling DiffOpFixAnyTime" << endl;
       Switch<2> (SpaceD-2, [&] (auto SD) {
         diffopfixt = make_shared<DiffOpFixAnyTime<SD+2>> (time);
       });
     }
     return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopfixt));
   },
   py::arg("gf"),
   py::arg("time") = 0.0,
    docu_string(R"raw_string(
fix_t fixes the time (ReferenceTimeVariable) of a given expression.
This is the variant for a gridfunction.

Parameters

self: ngsolve.GridFunction
  Gridfunction in which the time should be fixed
  
time: double
  Value the time should become
  
)raw_string")
);

   m.def("CreateTimeRestrictedGF", [](PyGF st_GF,double time) -> PyGF
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     return st_FES->CreateRestrictedGF(st_GF,time);
   },
   py::arg("gf"),
   py::arg("reference_time") = 0.0,
   "Create spatial-only Gridfunction corresponding to a fixed time.");

   m.def("RestrictGFInTime", [](PyGF st_GF,double time,PyGF s_GF)
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     switch (st_FES->GetDimension())
     {
       case 1:
         st_FES->RestrictGFInTime<double>(st_GF,time,s_GF);
         break;
       case 2:
         st_FES->RestrictGFInTime<Vec<2>>(st_GF,time,s_GF);
         break;
       case 3:
         st_FES->RestrictGFInTime<Vec<3>>(st_GF,time,s_GF);
         break;
       default:
         throw Exception("cannot handle GridFunction type (dimension too large?).");
         break;
     }
   }, 
   py::arg("spacetime_gf"),
   py::arg("reference_time") = 0.0,
   py::arg("space_gf"),
   "Extract Gridfunction corresponding to a fixed time from a space-time GridFunction.");

   m.def("SpaceTimeInterpolateToP1", [](PyCF st_CF, PyCF tref, PyGF st_GF)
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     if (!st_FES) throw Exception("not a spacetime gridfunction");
     st_FES->InterpolateToP1(st_CF,tref,st_GF);
   }, 
   py::arg("spacetime_cf"),
   py::arg("time"),
   py::arg("spacetime_gf"),
   "Interpolate nodal in time (possible high order) and nodal in space (P1).");


   py::class_<SpaceTimeVTKOutput, shared_ptr<SpaceTimeVTKOutput>>(m, "SpaceTimeVTKOutput")
    .def(py::init([] (shared_ptr<MeshAccess> ma, py::list coefs_list,
                      py::list names_list, string filename,
                      int subdivisionx, int subdivisiont, int only_element)
         -> shared_ptr<SpaceTimeVTKOutput>
         {
           Array<shared_ptr<CoefficientFunction> > coefs
             = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
           Array<string > names
             = makeCArray<string> (names_list);
           shared_ptr<SpaceTimeVTKOutput> ret;
           if (ma->GetDimension() == 2)
             ret = make_shared<SpaceTimeVTKOutput> (ma, coefs, names, filename, subdivisionx, subdivisiont, only_element);
           else
             throw Exception("Space time VTK only for 2D(+time)");
           return ret;
         }),
         py::arg("ma"),
         py::arg("coefs")= py::list(),
         py::arg("names") = py::list(),
         py::arg("filename") = "vtkout",
         py::arg("subdivision_x") = 0,
         py::arg("subdivision_t") = 0,
         py::arg("only_element") = -1
         )
     .def("Do", [](shared_ptr<SpaceTimeVTKOutput> self, VorB vb, double t_start, double t_end)
          { 
            size_t global_heapsize = 10000000;
            LocalHeap glh(global_heapsize, "spacetime-lh", true);
            self->Do(glh, vb, nullptr, t_start, t_end);
          },
          py::arg("vb")=VOL,
          py::arg("t_start") = 0,
          py::arg("t_end") = 1,
          py::call_guard<py::gil_scoped_release>())
     .def("Do", [](shared_ptr<SpaceTimeVTKOutput> self, VorB vb, double t_start, double t_end, const BitArray * drawelems)
          { 
            size_t global_heapsize = 10000000;
            LocalHeap glh(global_heapsize, "spacetime-lh", true);
            self->Do(glh, vb, drawelems, t_start, t_end);
          },
          py::arg("vb")=VOL,
          py::arg("t_start") = 0,
          py::arg("t_end") = 1,
          py::arg("drawelems"),
          py::call_guard<py::gil_scoped_release>())
     ;

   
}
