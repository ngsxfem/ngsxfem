#include <python_ngstd.hpp>
PYBIND11_PLUGIN(libngsxfem_osmosis) 
{
  py::module m("xfem_osmosis", "pybind xfem_osmosis");
  return m.ptr();
}
