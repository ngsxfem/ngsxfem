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
    on quads/hexes) level set GridFunctions, also in the isoparametric
    approach; space-time rules. Other level sets are only accepted with
    subdivlvl > 0 (single level sets, piecewise linear interpolation on
    subdivided simplices), otherwise an exception is raised.
ALGOIM:
    rules based on the Algoim library (R. I. Saye, https://algoim.github.io/),
    which represent the level set function(s) as element-wise polynomial(s)
    of higher degree (level set GridFunctions exactly, other
    CoefficientFunctions by interpolation, cf. ngsxfemglobals.algoim_lset_order).
    Available for element integrals (segments, triangles, quadrilaterals,
    tetrahedra, hexahedra, prisms; also on boundary elements, vb=BND), for
    facet and element boundary integrals (skeleton=True, element_boundary=True),
    for multiple level sets (DomainTypeArray, at most one interface condition
    per domain tuple), for space-time integrals (time_order >= 0, level sets on
    SpaceTimeFESpaces or CoefficientFunctions of tref; elements and skeleton
    facets) and at fixed reference times (tref). The standard rules are used
    otherwise (e.g. co-dimension 2 parts of multiple level set domains).
AUTO (default):
    the standard rules for piecewise linear (in space) level set GridFunctions,
    ALGOIM for all other level sets (higher order GridFunctions, general
    CoefficientFunctions) and for cases that the standard rules do not support
    (prisms; multiple level sets on quadrilaterals, hexahedra, prisms).
The use of Algoim is reported once per run in the terminal (with a link and
the reference to cite), which also mentions an automatic choice (AUTO, if not
chosen explicitly).
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
