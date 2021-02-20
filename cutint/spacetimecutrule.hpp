#pragma once
#include "xintegration.hpp"
#include "straightcutrule.hpp"
#include <algorithm>
#include <memory>
#include <numeric>

using namespace ngfem;


namespace xintegration
{
  tuple<const IntegrationRule *, Array<double>> SpaceTimeCutIntegrationRule(FlatVector<> cf_lset_at_element,
                                                     const ElementTransformation & trafo, //To be added
                                                     ScalarFiniteElement<1>* fe_time,
                                                     DOMAIN_TYPE dt,
                                                     int order_time,
                                                     int order_space,
                                                     SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                     LocalHeap & lh);
}
