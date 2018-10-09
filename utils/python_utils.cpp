#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../utils/bitarraycf.hpp"
#include "../utils/restrictedblf.hpp"
#include "../utils/p1interpol.hpp"

using namespace ngcomp;

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

  
  typedef shared_ptr<BitArray> PyBA;

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
        py::arg("flags") = py::dict(),
        docu_string(R"raw_string(
A restricted bilinear form is a (so far real-valued) bilinear form with a reduced MatrixGraph
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

space : ngsolve.FESpace
  finite element space on which the bilinear form is defined.

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
)raw_string"));

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

}

PYBIND11_MODULE(ngsxfem_utils_py,m)
{
  cout << "importing ngsxfem-utils lib" << endl;
  ExportNgsx_utils(m);
}
