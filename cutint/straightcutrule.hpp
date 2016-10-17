#pragma once
#include "xintegration.hpp"

using namespace ngfem;
namespace xintegration
{
  const IntegrationRule * StraightCutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                                     const ElementTransformation & trafo,
                                                     DOMAIN_TYPE dt,
                                                     int intorder,
                                                     LocalHeap & lh, vector<int>& sign_of_lset_at_vertex,
                                                     vector<int>& my_vert_idx);
}
