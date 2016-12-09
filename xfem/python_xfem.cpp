#include <python_ngstd.hpp>
PYBIND11_PLUGIN(libngsxfem_xfem) 
{
  py::module m("xfem_xfem", "pybind xfem_xfem");
  return m.ptr();
}
