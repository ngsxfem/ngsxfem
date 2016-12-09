#include <python_ngstd.hpp>
PYBIND11_PLUGIN(libngsxfem_tracefem) 
{
  py::module m("tracefem", "pybind tracefem");
  return m.ptr();
}


