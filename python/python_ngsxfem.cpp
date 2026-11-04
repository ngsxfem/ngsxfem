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

  py::enum_<CUT_QUADRATURE>(m, "CUT_QUADRATURE", docu_string(R"raw_string(
Quadrature method on cut elements (cf. ngsxfemglobals.cut_quadrature and the
argument cut_quadrature of dCut). Wherever a CUT_QUADRATURE is expected, the
strings "standard", "algoim" and "auto" can be used as well.

STANDARD:
    ngsxfem's standard rules: straight cuts of piecewise linear (multilinear
    on quads/hexes) level set approximations, subdivision for level set
    CoefficientFunctions, space-time and multiple level set rules.
ALGOIM:
    rules based on the Algoim library (R. I. Saye, https://algoim.github.io/),
    which represent the level set function as element-wise polynomial of
    higher degree (level set GridFunctions exactly, other CoefficientFunctions
    by interpolation, cf. ngsxfemglobals.algoim_lset_order). Used wherever
    available (single level sets on segments, triangles, quadrilaterals,
    tetrahedra and hexahedra; on boundary elements only for NEG/POS), the
    standard rules are used otherwise (multiple level sets, space-time).
AUTO:
    as ALGOIM, but only for higher order level sets (GridFunctions of order
    > 1, general CoefficientFunctions) and only on elements without mesh
    deformation (with a deformation, the isoparametric approach with a P1
    level set is assumed), the standard rules otherwise.
)raw_string"))
    .value("STANDARD", CUTQUAD_STANDARD)
    .value("ALGOIM", CUTQUAD_ALGOIM)
    .value("AUTO", CUTQUAD_AUTO)
    .def(py::init([](const string & name) { return CutQuadratureFromString(name); }))
    ;
  py::implicitly_convertible<py::str, CUT_QUADRATURE>();


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
