#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../utils/bitarraycf.hpp"
#include "../utils/restrictedblf.hpp"

using namespace ngcomp;

void ExportNgsx_utils(py::module &m)
{
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
        py::arg("flags") = py::dict()
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
  
}

PYBIND11_PLUGIN(ngsxfem_utils_py)
{
  cout << "importing ngsxfem-utils lib" << endl;
  py::module m("ngsxfem-utils", "pybind ngsxfem-utils");
  ExportNgsx_utils(m);
  return m.ptr();
}
