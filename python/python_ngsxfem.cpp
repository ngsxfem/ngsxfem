#include <core/utils.hpp>
#include <python_ngstd.hpp>

#include "../cutint/python_cutint.cpp"
#include "../utils/python_utils.cpp"
#include "../xfem/python_xfem.cpp"
#include "../spacetime/python_spacetime.cpp"
#include "../lsetcurving/python_lsetcurving.cpp"

#include "../utils/ngsxstd.hpp"

void ExportNgsx(py::module &m)
{

  py::enum_<DOMAIN_TYPE>(m, "DOMAIN_TYPE")
    .value("POS", POS)
    .value("NEG", NEG)
    .value("IF", IF)
    .export_values()
    ;

  py::enum_<COMBINED_DOMAIN_TYPE>(m, "COMBINED_DOMAIN_TYPE")
    .value("NO", CDOM_NO)
    .value("CDOM_NEG", CDOM_NEG)
    .value("CDOM_POS", CDOM_POS)
    .value("UNCUT", CDOM_UNCUT)
    .value("CDOM_IF", CDOM_IF)
    .value("HASNEG", CDOM_HASNEG)
    .value("HASPOS", CDOM_HASPOS)
    .value("ANY", CDOM_ANY)
    .export_values()
    ;
  
  py::enum_<SWAP_DIMENSIONS_POLICY>(m, "QUAD_DIRECTION_POLICY")
    .value("FIRST", FIRST_ALLOWED)
    .value("OPTIMAL", FIND_OPTIMAL)
    .value("FALLBACK", ALWAYS_NONE)
    .export_values()
    ;

  py::enum_<TIME_DOMAIN_TYPE>(m, "TIME_DOMAIN_TYPE")
    .value("BOTTOM", BOTTOM)
    .value("TOP", TOP)
    .value("INTERVAL", INTERVAL)
    .export_values()
    ;


}

PYBIND11_MODULE(xfem, m)
{
  cout << "importing ngsxfem-" << NGSXFEM_VERSION << endl;
  ExportNgsx(m);
  ExportNgsx_cutint(m);
  ExportNgsx_utils(m);
  ExportNgsx_xfem(m);
  ExportNgsx_spacetime(m);
  ExportNgsx_lsetcurving(m);
}
