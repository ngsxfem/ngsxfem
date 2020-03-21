#pragma once

#include "xintegration.hpp"

using namespace ngfem;
using ngfem::ELEMENT_TYPE;

namespace xintegration
{

  /// Mlset-integrationRule Factory
  tuple<const IntegrationRule *, Array<double> > CreateCutIntegrationRule(Array<shared_ptr<GridFunction>> & gflsets,
                                                                          const ElementTransformation & trafo,
                                                                          Array<DOMAIN_TYPE> & dts,
                                                                          int intorder,
                                                                          int time_intorder,
                                                                          LocalHeap & lh,
                                                                          SWAP_DIMENSIONS_POLICY quad_dir_policy = FIND_OPTIMAL);

  tuple<const IntegrationRule *, Array<double> > CreateMultiLevelsetCutIntegrationRule(const LevelsetIntegrationDomain & lsetintdom,
                                                                                       const ElementTransformation & trafo,
                                                                                       LocalHeap & lh);
}

