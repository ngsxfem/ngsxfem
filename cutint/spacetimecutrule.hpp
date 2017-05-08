#pragma once
#include "xintegration.hpp"
#include "straightcutrule.hpp"
#include <algorithm>
#include <memory>
#include <numeric>

using namespace ngfem;


namespace xintegration
{
  void DebugSpaceTimeCutIntegrationRule();

  const IntegrationRule * SpaceTimeCutIntegrationRule(FlatVector<> cf_lset_at_element,
                                                     //const ElementTransformation & trafo, //To be added
                                                     ELEMENT_TYPE et_space, DOMAIN_TYPE dt,
                                                     int order_time,
                                                     int order_space,
                                                     LocalHeap & lh);
}
