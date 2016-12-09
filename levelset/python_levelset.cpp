#include <python_ngstd.hpp>
PYBIND11_PLUGIN(libngsxfem_levelset) 
{
  py::module m("xfem_levelset", "pybind xfem_levelset");
  return m.ptr();
}
