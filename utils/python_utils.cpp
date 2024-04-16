#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../utils/bitarraycf.hpp"
#include "../utils/restrictedblf.hpp"
#include "../utils/p1interpol.hpp"
#include "../utils/xprolongation.hpp"
#include "../utils/restrictedfespace.hpp"
#include "../utils/ngsxstd.hpp"

GlobalNgsxfemVariables globxvar;

using namespace ngcomp;
typedef shared_ptr<BitArray> PyBA;

 auto rblf_string_T = docu_string(R"raw_string(
A restricted bilinear form is a bilinear form with a reduced MatrixGraph
compared to the usual BilinearForm. BitArray(s) define on which elements/facets entries will be
created.

Use cases:

 * ghost penalty type stabilization:
    Facet-stabilization that are introduced only act on a few facets in the mesh. By providing the
    information on the corresponding facets, these additional couplings will only be introduced
    where necessary.

 * fictitious domain methods:
    When PDE problems are only solved on a part of a domain while a finite element space is used
    that is still defined on the whole domain, a BitArray can be used to mark the 'active' part of
    the mesh.

Parameters

space (trialspace) : ngsolve.FESpace
  finite element space on which the bilinear form is defined 
  (trial space and (if no test space is defined) test space).

testspace : ngsolve.FESpace
  finite element space on which the bilinear form is defined
  (test space).

name : string
  name of the bilinear form

element_restriction : ngsolve.BitArray
  BitArray defining the 'active mesh' element-wise

facet_restriction : ngsolve.BitArray
  BitArray defining the 'active facets'. This is only relevant if FESpace has DG-terms (dgjumps=True)

kwargs : keyword args 
  additional arguments that are passed to bilinear form (in form of flags)
)raw_string");

template <class TM,class TV>
void declare_RestrictedBilinearForm(py::module &m, std::string const &typestr) {
using Rbfi_TT = RestrictedBilinearForm<TM,TV>;
std::string pyclass_name_T = std::string("RestrictedBilinearForm") + typestr;

py::class_<Rbfi_TT, shared_ptr<Rbfi_TT>, BilinearForm> rblf_T(m, pyclass_name_T.c_str(),docu_string(R"raw_string(BilinearForm restricted on a set of elements and facets.
)raw_string") , py::dynamic_attr());

rblf_T.def(py::init([](shared_ptr<FESpace> fes,
      const string & aname,
      py::object ael_restriction,
      py::object afac_restriction,
      py::kwargs kwargs)
  {
    auto flags = CreateFlagsFromKwArgs(kwargs);

    shared_ptr<BitArray> el_restriction = nullptr;
    shared_ptr<BitArray> fac_restriction = nullptr;

    if ((!ael_restriction.is_none()) && py::extract<PyBA> (ael_restriction).check())
      el_restriction = py::extract<PyBA>(ael_restriction)();

    if ((!afac_restriction.is_none()) && py::extract<PyBA> (afac_restriction).check())
      fac_restriction = py::extract<PyBA>(afac_restriction)();

    auto biform = make_shared<Rbfi_TT> (fes, aname, el_restriction, fac_restriction, flags);
    return biform;
  }),
  py::arg("space"),
  py::arg("name") = "bfa",
  py::arg("element_restriction")=py::none(),
  py::arg("facet_restriction")=py::none(),
  rblf_string_T)
.def(py::init([](shared_ptr<FESpace> fes1,
   shared_ptr<FESpace> fes2,
   const string & aname,
   py::object ael_restriction,
   py::object afac_restriction,
   py::kwargs kwargs)
{
   auto flags = CreateFlagsFromKwArgs(kwargs);
   shared_ptr<BitArray> el_restriction = nullptr;
   shared_ptr<BitArray> fac_restriction = nullptr;
   if ((!ael_restriction.is_none()) && py::extract<PyBA> (ael_restriction).check())
     el_restriction = py::extract<PyBA>(ael_restriction)();

   if ((!afac_restriction.is_none()) && py::extract<PyBA> (afac_restriction).check())
     fac_restriction = py::extract<PyBA>(afac_restriction)();

   auto biform = make_shared<Rbfi_TT> (fes1, fes2, aname, el_restriction, fac_restriction, flags);
   return biform;
}),
py::arg("trialspace"),
py::arg("testspace"),
py::arg("name") = "bfa",
py::arg("element_restriction")=py::none(),
py::arg("facet_restriction")=py::none(),
rblf_string_T)        
.def_property("element_restriction", 
	  &Rbfi_TT::GetElementRestriction,
	  &Rbfi_TT::SetElementRestriction, "element restriction")
.def_property("facet_restriction", 
	  &Rbfi_TT::GetFacetRestriction,
	  &Rbfi_TT::SetFacetRestriction, "facet restriction");
}



void ExportNgsx_utils(py::module &m)
{
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;

  m.def("InterpolateToP1",  [] (PyGF gf_ho, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(gf_ho, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh, eps_perturbation);
        } ,
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,
        py::arg("eps_perturbation")=globxvar.EPS_INTERPOLATE_TO_P1,
        py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Takes the vertex values of a GridFunction (also possible with a CoefficentFunction) and puts them
into a piecewise (multi-) linear function.

Parameters

gf_ho : ngsolve.GridFunction
  Function to interpolate

gf_p1 : ngsolve.GridFunction
  Function to interpolate to (should be P1)

eps_perturbation : float
  If the absolute value if the function is smaller than eps_perturbation, it will be set to
  eps_perturbation. Thereby, exact and close-to zeros at vertices are avoided (Useful to reduce cut
  configurations for level set based methods).

heapsize : int
  heapsize of local computations.
)raw_string")
    )
    ;

  m.def("InterpolateToP1",  [] (PyCF coef, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(coef, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh, eps_perturbation);
        } ,
        py::arg("coef"),py::arg("gf"),
        py::arg("eps_perturbation")=globxvar.EPS_INTERPOLATE_TO_P1, py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Takes the vertex values of a CoefficentFunction) and puts them into a piecewise (multi-) linear
function.

Parameters

coef : ngsolve.CoefficientFunction
  Function to interpolate

gf_p1 : ngsolve.GridFunction
  Function to interpolate to (should be P1)

eps_perturbation : float
  If the absolute value if the function is smaller than eps_perturbation, it will be set to
  eps_perturbation. Thereby, exact and close-to zeros at vertices are avoided (Useful to reduce cut
  configurations for level set based methods).

heapsize : int
  heapsize of local computations.
)raw_string")
    )
    ;

  // Export RestrictedBilinearForm for different data types 
  declare_RestrictedBilinearForm<double,double>(m,"Double");
  declare_RestrictedBilinearForm<Complex,Complex>(m,"Complex");

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
                res->SetBit(offset+i);
            }
            offset += ba->Size();
          }
          return res;
        } ,
        py::arg("balist"),
        docu_string(R"raw_string(
Takes a list of BitArrays and merges them to one larger BitArray. Can be useful for
CompoundFESpaces.
)raw_string")
    );

  py::class_<GlobalNgsxfemVariables>(m, "GlobalNgsxfemVariables",
              docu_string(R"raw_string(
The class GlobalNgsxfemVariables provides Python-access to several internal
parameters and options used by different subprocedures of ngsxfem. For "mainstream"
application cases, it should not be required to change parameters here. Most cases
where this class is practically relevant will be debugging or special applications,
like investigations in a regime of total error below ~1e-8.

Properties:

eps_spacetime_lset_perturbation : double
    When handling cut topologies, it is sometimes cumbersome to include the case
    of a lset value of exactly 0. Hence, the value will be set to eps_spacetime_lset_perturbation
    in the routine for generating space-time quadrature rules in case its absolute value is smaller.
    Default: 1e-14

eps_spacetime_cutrule_bisection : double
    For high temporal orders, the space-time quadrature rule will apply a bisection
    method to find those time points with topology changes. This parameters controls
    how small 2 times the value must be in order to be counted as a root.
    Default: 1e-15

eps_P1_perturbation : double
    Similar to eps_spacetime_lset_perturbation, but for the P1 interpolation routine.
    Default: 1e-14

eps_spacetime_fes_node : double
    When a Gridfunction is restricted, the given time point is compared to the nodes
    of the finite element, such that those node values can be extracted directly in
    a matching case. This parameters controlls how far a deviation will still be counted
    as coincidence.
    Default: 1e-9


 )raw_string"))
          .def_readwrite("eps_spacetime_lset_perturbation", &GlobalNgsxfemVariables::EPS_STCR_LSET_PERTUBATION)
          .def_readwrite("eps_spacetime_cutrule_bisection", &GlobalNgsxfemVariables::EPS_STCR_ROOT_SEARCH_BISECTION)
          .def_readwrite("eps_P1_perturbation", &GlobalNgsxfemVariables::EPS_INTERPOLATE_TO_P1)
          .def_readwrite("eps_spacetime_fes_node", &GlobalNgsxfemVariables::EPS_STFES_RESTRICT_GF)
          .def_readwrite("eps_shifted_eval", &GlobalNgsxfemVariables::EPS_SHIFTED_EVAL)
          .def_readwrite("eps_facetpatch_ips", &GlobalNgsxfemVariables::EPS_FACET_PATCH_INTEGRATOR)
          .def_readwrite("newton_maxiter", &GlobalNgsxfemVariables::NEWTON_ITER_TRESHOLD)
          .def_readwrite("max_dist_newton", &GlobalNgsxfemVariables::MAX_DIST_NEWTON)
          .def_readwrite("fixed_point_maxiter_shifted_eval", &GlobalNgsxfemVariables::FIXED_POINT_ITER_TRESHOLD)
          .def_readwrite("do_naive_timeint", &GlobalNgsxfemVariables::DO_NAIVE_TIMEINT)
          .def_readwrite("naive_timeint_order", &GlobalNgsxfemVariables::NAIVE_TIMEINT_ORDER)
          .def_readwrite("naive_timeint_subdivs", &GlobalNgsxfemVariables::NAIVE_TIMEINT_SUBDIVS)
          .def_readwrite("non_conv_warn_msg_lvl", &GlobalNgsxfemVariables::NON_CONV_WARN_MSG_LVL)
          .def_readwrite("simd_eval", &GlobalNgsxfemVariables::SIMD_EVAL)
          .def("MultiplyAllEps", &GlobalNgsxfemVariables::MultiplyAllEps)
          .def("Output", &GlobalNgsxfemVariables::Output)
          .def("SetDefaults", &GlobalNgsxfemVariables::SetDefaults)
          .def("SwitchSIMD", &GlobalNgsxfemVariables::SwitchSIMD);

  
  m.attr("ngsxfemglobals") = py::cast(&globxvar);

  typedef shared_ptr<BitArrayCoefficientFunction> PyBACF;
  py::class_<BitArrayCoefficientFunction, PyBACF, CoefficientFunction>
    (m, "BitArrayCF",
        docu_string(R"raw_string(
CoefficientFunction that evaluates a BitArray. On elements with an index i where the BitArray
evaluates to true the CoefficientFunction will evaluate as 1, otherwise as 0.

Similar functionality (also for facets) can be obtained with IndicatorCF.
)raw_string"))
    .def("__init__",
         [](BitArrayCoefficientFunction *instance, shared_ptr<BitArray> ba)
         {
           new (instance) BitArrayCoefficientFunction (ba);
         },
         py::arg("bitarray")
      );

  py::class_<RestrictedFESpace, shared_ptr<RestrictedFESpace>, CompressedFESpace>(m, "Restrict",
	docu_string(R"delimiter(Wrapper Finite Element Spaces.
The restricted fespace is a wrapper around a standard fespace which removes dofs from marked elements.

Parameters:

fespace : ngsolve.comp.FESpace
    finite element space

active_els : BitArray or None
    Only use dofs from these elements
)delimiter"))
    .def(py::init([] (shared_ptr<FESpace> & fes,
                      py::object active_els)
                  {
                    // shared_ptr<CompoundFESpace> compspace = dynamic_pointer_cast<CompoundFESpace> (fes);
                    // if (compspace)
                      // cout << "yes, we can also compress a CompoundFESpace" << endl;
                    // throw py::type_error("cannot apply compression on CompoundFESpace - Use CompressCompound(..)");
                    auto ret = make_shared<RestrictedFESpace> (fes);
                    shared_ptr<BitArray> actdofs = nullptr;
                    if (!active_els.is_none())
                      dynamic_pointer_cast<RestrictedFESpace>(ret)->SetActiveElements(py::extract<shared_ptr<BitArray>>(active_els)());
                    ret->Update();
                    ret->FinalizeUpdate();
                    return ret;                    
                  }), py::arg("fespace"), py::arg("active_elements")=py::none())
    .def("GetBaseSpace", [](RestrictedFESpace & self)
         {
           return self.GetBaseSpace();
         })
    .def(py::pickle([](const RestrictedFESpace* restr_fes)
                    {
                      return py::make_tuple(restr_fes->GetBaseSpace(),restr_fes->GetActiveElements());
                    },
                    [] (py::tuple state) -> shared_ptr<RestrictedFESpace>
                    {
                      auto fes = make_shared<RestrictedFESpace>(state[0].cast<shared_ptr<FESpace>>());
                      if (state[1].cast<shared_ptr<BitArray>>())
                        fes->SetActiveElements(state[1].cast<shared_ptr<BitArray>>());
                      fes->Update();
                      fes->FinalizeUpdate();
                      return fes;
                    }))
    .def_property("active_elements", 
                    &RestrictedFESpace::GetActiveElements,
                    &RestrictedFESpace::SetActiveElements, "active elements")

    ;









  typedef shared_ptr<P1Prolongation> PyP1P;
  py::class_<P1Prolongation, PyP1P, Prolongation>
    (m, "P1Prolongation",
        docu_string(R"raw_string(
Prolongation for P1-type spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P1Prolongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P1Prolongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P1Prolongation> p1p, shared_ptr<FESpace> fes)
         {
           p1p -> Update(*fes);
         },
         py::arg("space")
      );

      typedef shared_ptr<P2Prolongation> PyP2P;
  py::class_<P2Prolongation, PyP2P, Prolongation>
    (m, "P2Prolongation",
        docu_string(R"raw_string(
Prolongation for P2 spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P2Prolongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P2Prolongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P2Prolongation> p2p, shared_ptr<FESpace> fes)
         {
           p2p -> Update(*fes);
         },
         py::arg("space")
      );


      typedef shared_ptr<P2CutProlongation> PyP2CutP;
  py::class_<P2CutProlongation, PyP2CutP, Prolongation>
    (m, "P2CutProlongation",
        docu_string(R"raw_string(
Prolongation for P2 spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P2CutProlongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P2CutProlongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P2CutProlongation> p2p, shared_ptr<FESpace> fes)
         {
           p2p -> Update(*fes);
         },
         py::arg("space")
      );

    typedef shared_ptr<CompoundProlongation> PyCProl;
    py::class_< CompoundProlongation, PyCProl, Prolongation>
    ( m, "CompoundProlongation", 
     docu_string(R"raw_string(prolongation for compound spaces)raw_string"))
     .def("__init__",
          [](CompoundProlongation *instance, const FESpace *fes)
          {
            new (instance) CompoundProlongation( dynamic_cast<const CompoundFESpace*>(fes) );
          },
          py::arg("compoundFESpace"))
      .def("Update", &CompoundProlongation::Update, py::arg("fespace"))
      .def ("Prolongate", &CompoundProlongation::ProlongateInline, py::arg("finelevel"), py::arg("vec"))
      .def ("Restrict", &CompoundProlongation::RestrictInline, py::arg("finelevel"), py::arg("vec"))
      .def ("AddProlongation", [](shared_ptr<CompoundProlongation> cprol, shared_ptr<Prolongation> prol )
            { cprol -> AddProlongation( prol ); }, py::arg("p1prol")
          );
      //.def ("AddProlongation" &CompoundProlongation::AddProlongation, py::arg("prol"));


}
