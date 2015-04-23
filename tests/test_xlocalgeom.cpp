#define BOOST_TEST_MAIN
#if !defined( WIN32 )
    #define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>

#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

/// from ngxfem
#include "../cutint/xintegration.hpp"

// --- Cut situation to be tested: ---
//
//   (0,1) 
//     \
//     |\
//     | \
//     |  \
//     |   \
//     +    \
//     |\    \
//     | \    \
//     |  \    \
//     |NEG\ POS\
//     +----+----+
//   (0,0)     (1,0)
//
// meas_2(POS) = 3/8, meas_2(NEG) = 1/8, meas_1(IF) = 1/2

BOOST_AUTO_TEST_CASE( xlocalgeom ) {

  ngstd::LocalHeap lh(10000, "test_xlocal");
  ngfem::ELEMENT_TYPE eltype = ET_TRIG;
  ngfem::ELEMENT_TYPE et_time = ET_POINT;

  ngfem::H1HighOrderFE<ET_TRIG> fe(1);
  FlatVector<> fv(3,lh);

  fv(0) = -1; fv(1) = 1; fv(2) = 1;

  const int order_space = 1;
  const int order_time = 0;
  const int ref_lvl_space = 0;
  const int ref_lvl_time = 0;

  xintegration::CompositeQuadratureRule<2> cquad;

  //linear function which is -1 at (0,0) and 1 at the other vertices...
  ScalarFieldEvaluator* lset_eval_p = ScalarFieldEvaluator::Create(2, fe, fv, lh); //allocated on lh
  

  auto xgeom
    = xintegration::XLocalGeometryInformation::Create
      (eltype, et_time, *lset_eval_p,
       cquad, lh,
       2*order_space, 2*order_time,
       ref_lvl_space, ref_lvl_time); //memory has to be taken care of...

  xintegration::DOMAIN_TYPE dt = xgeom->MakeQuadRule();

  BOOST_CHECK_EQUAL(dt, xintegration::IF);


  double int_neg = 0.0;
  const xintegration::QuadratureRule<2> & quad_neg(cquad.GetRule(xintegration::NEG));
  for (int i = 0; i < quad_neg.Size(); ++i)
      int_neg += quad_neg.weights[i];

  BOOST_CHECK_CLOSE(int_neg, 0.125, 1e-12);

  double int_pos = 0.0;
  const xintegration::QuadratureRule<2> & quad_pos(cquad.GetRule(xintegration::POS));
  for (int i = 0; i < quad_pos.Size(); ++i)
      int_pos += quad_pos.weights[i];

  BOOST_CHECK_CLOSE(int_pos, 0.375, 1e-12);

  double int_if = 0.0;
  const xintegration::QuadratureRuleCoDim1<2> & quad_if(cquad.GetInterfaceRule());
  for (int i = 0; i < quad_if.Size(); ++i)
      int_if += quad_if.weights[i];

  BOOST_CHECK_CLOSE(int_if, 0.5, 1e-12);
}
