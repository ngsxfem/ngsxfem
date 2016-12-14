#include <python_ngstd.hpp>
PYBIND11_PLUGIN(libngsxfem_xstokes) 
{
  py::module m("xfem_xstokes", "pybind xfem_xstokes");
  return m.ptr();
}
