#pragma once

/*
  Cut integration rules based on the bundled Algoim library
  (R. I. Saye, https://algoim.github.io/, see external_dependencies/algoim).

  The level set function is represented on each element by a polynomial in
  Bernstein form on the unit hypercube [0,1]^D, which contains the reference
  element (only affine transformations are involved):
    * quadrilaterals / hexahedra (and segments): the reference element is the
      unit hypercube itself;
    * triangles / tetrahedra: the reference simplex {x_0+...+x_{D-1} < 1} is
      the corner of the unit hypercube; it is described by a second (linear)
      polynomial in Algoim's multi-polynomial quadrature, and the level set
      function is evaluated (polynomially extended) in the whole hypercube.
      The hypercube is subdivided once (the sub-boxes that intersect the
      simplex are used) to keep the boxes Algoim works on comparable in size
      to the simplex.
      Since only affine mappings are involved, the quadrature rules are
      exact for polynomial integrands on domains cut by planes. (Collapsed or
      multilinear maps of the hypercube onto (parts of) the simplex render the
      integrands non-polynomial and spoil the h-convergence of the rules.)
  The polynomial (of degree k in each variable for a polynomial of total
  degree k on simplices) is obtained by interpolation in Algoim's
  interpolation points. This is exact for level set GridFunctions, whose
  element-wise degree is used. Level set CoefficientFunctions that are not
  GridFunctions are interpolated with a prescribed degree.

  The returned integration rules follow the conventions of
  StraightCutIntegrationRule: points and weights refer to the reference
  element; scaling a weight with the measure of the element transformation
  gives the physical weight. For interface rules (dt == IF) this means that the
  weight of a reference surface quadrature point x with reference normal n is
  multiplied by |F^{-T}(x) n|, F the Jacobian of the element transformation.
  On boundary elements, the interface has co-dimension 2 in physical space;
  the weights are such that the scaling with the measure of the (boundary)
  element transformation gives the weight w.r.t. the physical point / length
  measure.

  Space-time: the level set function (a GridFunction on a SpaceTimeFESpace or
  a CoefficientFunction of the reference time variable) is represented on the
  space-time cell element x [0,1] (reference time) by a polynomial in D+1
  variables (for triangles / tetrahedra: on [0,1]^(D+1), the prism is the
  product of the simplex with the time interval). Algoim's quadrature is
  applied in space-time (up to 3D+t, i.e. dimension 4); the time is preferred
  as the variable of the outermost integral (see
  algoim_wrapper::MultiCutQuadrature). The rules follow the conventions of the
  standard space-time rules (SpaceTimeCutIntegrationRule).
*/

#include "xintegration.hpp"

namespace xintegration
{
  /// Algoim-based cut rules are available for (volume and boundary) elements of
  /// type SEGM, TRIG, QUAD, TET, HEX.
  bool AlgoimCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt);

  /// Cut rule for a domain described by multiple level sets (the union of the
  /// parts given by the domain type tuples of lsetintdom; at most one interface
  /// condition per tuple) computed with Algoim (Algoim's multi-polynomial
  /// quadrature w.r.t. all level sets that cut the element). Conventions as for
  /// AlgoimCutIntegrationRule.
  tuple<const IntegrationRule *, Array<double>>
  AlgoimMultiLevelsetCutIntegrationRule (const LevelsetIntegrationDomain & lsetintdom,
                                         const ElementTransformation & trafo,
                                         int cf_degree, LocalHeap & lh);

  /// Algoim-based cut rules on facets are available for facets of type POINT,
  /// SEGM, TRIG, QUAD.
  bool AlgoimCutFacetRuleAvailable (ELEMENT_TYPE etfacet);

  /// Cut rule on the facet facetnr of the element (trafo, transform: the map from
  /// the facet to the element) computed with Algoim w.r.t. the level set function
  /// restricted to the facet. Points are given in reference coordinates of the facet,
  /// with the weight conventions of the cut facet integrators (cf.
  /// StraightCutIntegrationRuleUntransformed):
  ///  * NEG/POS: weights w.r.t. the reference facet (to be scaled with the measure of
  ///    the facet),
  ///  * IF: weights w.r.t. the physical measure of the interface on the facet
  ///    (points on 1D facets: weight 1, curves on 2D facets: physical length element).
  /// Returns nullptr if the part is empty.
  const IntegrationRule *
  AlgoimCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                 const ElementTransformation & trafo,
                                 Facet2ElementTrafo & transform, int facetnr,
                                 DOMAIN_TYPE dt, int intorder, int cf_degree, int subdivlvl,
                                 LocalHeap & lh);

  /// Space-time version of AlgoimCutFacetIntegrationRule: cut rule on the facet facetnr
  /// times the reference time interval [0,1] w.r.t. the space-time level set function (a
  /// GridFunction on a SpaceTimeFESpace or a CoefficientFunction that depends on the time
  /// variable), with the conventions of SpaceTimeCutIntegrationRuleUntransformed: points in
  /// reference coordinates of the facet, the reference time in the weight of the (space-time)
  /// integration points, the weights in the array (NEG/POS: w.r.t. the reference facet; IF:
  /// w.r.t. the iterated physical measure int_0^1 int_{Gamma(t) on the facet} dt).
  tuple<const IntegrationRule *, Array<double>>
  AlgoimSpaceTimeCutFacetIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                          const ElementTransformation & trafo,
                                          Facet2ElementTrafo & transform, int facetnr,
                                          DOMAIN_TYPE dt, int intorder, int time_intorder,
                                          int cf_degree, int subdivlvl, LocalHeap & lh);

  /// Polynomial degree that is used to represent the level set function on
  /// the element: the element-wise order of a (scalar) level set GridFunction,
  /// otherwise cf_degree. For GridFunctions on a SpaceTimeFESpace: the order in
  /// space, or (spacetime == true) the maximum of the orders in space and time.
  int AlgoimLevelsetDegree (const shared_ptr<CoefficientFunction> & cf_lset,
                            const ElementTransformation & trafo,
                            int cf_degree,
                            LocalHeap & lh,
                            bool spacetime = false);

  /// Algoim-based space-time cut rules are available for the element types of
  /// AlgoimCutRuleAvailable (space-time interfaces only on volume elements).
  bool AlgoimSpaceTimeCutRuleAvailable (const ElementTransformation & trafo, DOMAIN_TYPE dt);

  /// Cut integration rule (and weights) for the part of the element given by
  /// dt (NEG: {phi<0}, POS: {phi>0}, IF: {phi=0}) computed with Algoim.
  /// intorder: integration order (for polynomial integrands and planar cuts
  ///   the rule is exact up to this order),
  /// cf_degree: interpolation degree for level set CoefficientFunctions that
  ///   are not GridFunctions,
  /// subdivlvl: number of uniform subdivisions of the element (of its
  ///   bounding box [0,1]^D) before Algoim is applied, see
  ///   algoim_wrapper::CutQuadrature (useful on meshes that are coarse
  ///   compared to the curvature of the interface).
  /// tref: fixed reference time for space-time level set functions (then, the
  ///   reference time is stored in the weights of the (space-time) integration
  ///   points, the weights are given in the array only).
  /// Returns (nullptr, empty array) if the part is empty. Uncut elements get the
  /// standard integration rule of the element.
  tuple<const IntegrationRule *, Array<double>>
  AlgoimCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                            const ElementTransformation & trafo,
                            DOMAIN_TYPE dt,
                            int intorder,
                            int cf_degree,
                            int subdivlvl,
                            LocalHeap & lh,
                            optional<double> tref = nullopt);

  /// Space-time cut integration rule on the element times the reference time
  /// interval [0,1] w.r.t. the space-time level set function (a GridFunction on a
  /// SpaceTimeFESpace or a CoefficientFunction that depends on the time variable),
  /// represented by a polynomial on [0,1]^(D+1) (of degree max(order in space,
  /// order in time) per variable). Algoim's quadrature is applied in space-time,
  /// the conventions are those of SpaceTimeCutIntegrationRule: the integration
  /// points are spatial (reference) points, with the reference time in the weight
  /// of the (space-time) integration point, the weights are given in the array.
  /// For interfaces (dt == IF), the weights refer to the iterated measure
  /// int_0^1 int_{Gamma(t)} ds dt (as for the standard space-time rules), i.e. the
  /// weights of the space-time surface rule are multiplied with
  /// |grad_x phi| / |grad_(x,t) phi| (coarea formula).
  /// intorder / time_intorder: integration orders in space / time.
  tuple<const IntegrationRule *, Array<double>>
  AlgoimSpaceTimeCutIntegrationRule (const shared_ptr<CoefficientFunction> & cf_lset,
                                     const ElementTransformation & trafo,
                                     DOMAIN_TYPE dt,
                                     int intorder,
                                     int time_intorder,
                                     int cf_degree,
                                     int subdivlvl,
                                     LocalHeap & lh);
}
