#pragma once

#include <fem.hpp>

const int SPACETIME_SANITY_CHECK_NR = -9;

INLINE void MarkAsSpaceTimeIntegrationPoint(ngfem::IntegrationPoint & ip)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  ip.SetNr(SPACETIME_SANITY_CHECK_NR);
#endif  
}

INLINE void MarkAsSpaceTimeIntegrationRule(ngfem::IntegrationRule & ir)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  for (auto & ip : ir)
    ip.SetNr(SPACETIME_SANITY_CHECK_NR);
#endif  
}

INLINE bool IsSpaceTimeIntegrationPoint(const ngfem::IntegrationPoint & ip)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  return (ip.Nr() == SPACETIME_SANITY_CHECK_NR);
#else
  return true;
#endif
}
