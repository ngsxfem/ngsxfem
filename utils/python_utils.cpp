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

check_unused : boolean
  Check if some degrees of freedoms are not considered during assembly

flags : ngsolve.Flags
  additional bilinear form flags
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
      bool check_unused,
      py::kwargs kwargs)
  {
    auto flags = CreateFlagsFromKwArgs(kwargs);

    shared_ptr<BitArray> el_restriction = nullptr;
    shared_ptr<BitArray> fac_restriction = nullptr;
    if (py::extract<PyBA> (ael_restriction).check())
      el_restriction = py::extract<PyBA>(ael_restriction)();

    if (py::extract<PyBA> (afac_restriction).check())
      fac_restriction = py::extract<PyBA>(afac_restriction)();

    auto biform = make_shared<Rbfi_TT> (fes, aname, el_restriction, fac_restriction, flags);
    biform -> SetCheckUnused (check_unused);
    return biform;
  }),
  py::arg("space"),
  py::arg("name") = "bfa",
  py::arg("element_restriction") = DummyArgument(),
  py::arg("facet_restriction") = DummyArgument(),
  py::arg("check_unused") = true,
  rblf_string_T)
.def(py::init([](shared_ptr<FESpace> fes1,
   shared_ptr<FESpace> fes2,
   const string & aname,
   py::object ael_restriction,
   py::object afac_restriction,
   bool check_unused,
   py::kwargs kwargs)
{
   auto flags = CreateFlagsFromKwArgs(kwargs);

   shared_ptr<BitArray> el_restriction = nullptr;
   shared_ptr<BitArray> fac_restriction = nullptr;
   if (py::extract<PyBA> (ael_restriction).check())
     el_restriction = py::extract<PyBA>(ael_restriction)();

   if (py::extract<PyBA> (afac_restriction).check())
     fac_restriction = py::extract<PyBA>(afac_restriction)();

   auto biform = make_shared<Rbfi_TT> (fes1, fes2, aname, el_restriction, fac_restriction, flags);
   biform -> SetCheckUnused (check_unused);
   return biform;
}),
py::arg("trialspace"),
py::arg("testspace"),
py::arg("name") = "bfa",
py::arg("element_restriction") = DummyArgument(),
py::arg("facet_restriction") = DummyArgument(),
py::arg("check_unused") = true,
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
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,
        py::arg("eps_perturbation")=1e-14,py::arg("heapsize")=1000000,
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
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("coef"),py::arg("gf"),
        py::arg("eps_perturbation")=1e-14,py::arg("heapsize")=1000000,
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
                    if (! py::extract<DummyArgument> (active_els).check())
                      dynamic_pointer_cast<RestrictedFESpace>(ret)->SetActiveElements(py::extract<shared_ptr<BitArray>>(active_els)());
                    ret->Update();
                    ret->FinalizeUpdate();
                    return ret;                    
                  }), py::arg("fespace"), py::arg("active_elements")=DummyArgument())
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

