#pragma once

// #undef LAPACK
/// from ngsolve
#include <comp.hpp>

using namespace ngbla;

enum SWAP_DIMENSIONS_POLICY {FIRST_ALLOWED, FIND_OPTIMAL, ALWAYS_NONE};

// /// domain types: two domains: POS/NEG and the diving interface IF
// enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};

/// domain types: two domains, POS/NEG and the diving interface IF
enum DOMAIN_TYPE {  NEG = 0,
                    POS = 1,
                    IF = 2,
};

/// domain types: two domains, POS/NEG and the diving interface IF and all combinations of these:
enum COMBINED_DOMAIN_TYPE {  CDOM_NO = 0,
                             CDOM_NEG = 1,
                             CDOM_POS = 2,
                             CDOM_UNCUT = 3,
                             CDOM_IF = 4,
                             CDOM_HASNEG = 5,
                             CDOM_HASPOS = 6,
                             CDOM_ANY = 7
};

const int N_COMBINED_DOMAIN_TYPES = 8;
const ArrayMem<COMBINED_DOMAIN_TYPE,
            N_COMBINED_DOMAIN_TYPES> all_cdts = {CDOM_NO,CDOM_NEG,CDOM_POS,CDOM_UNCUT,
                                                 CDOM_IF,CDOM_HASNEG,CDOM_HASPOS,CDOM_ANY};



INLINE COMBINED_DOMAIN_TYPE TO_CDT(DOMAIN_TYPE a) {
  if (a == NEG)
    return CDOM_NEG;
  else if (a == POS)
    return CDOM_POS;
  else
    return CDOM_IF;
}

INLINE COMBINED_DOMAIN_TYPE operator|(DOMAIN_TYPE a,DOMAIN_TYPE b)
{
  auto p = char(TO_CDT(a));
  auto q = char(TO_CDT(b));
  return COMBINED_DOMAIN_TYPE(p|q);
}

INLINE COMBINED_DOMAIN_TYPE operator&(DOMAIN_TYPE a,DOMAIN_TYPE b)
{
  auto p = char(TO_CDT(a));
  auto q = char(TO_CDT(b));
  return COMBINED_DOMAIN_TYPE(p&q);
}

INLINE COMBINED_DOMAIN_TYPE operator~(COMBINED_DOMAIN_TYPE a)
{
  return COMBINED_DOMAIN_TYPE(~char(a));
}

INLINE DOMAIN_TYPE INVERT( DOMAIN_TYPE dt)
{
  if (dt == IF)
  {
    return IF;
  }
  else
  {
    if (dt == NEG)
      return POS;
    else
      return NEG;
  }
}

/// time domain types: Interval or one of the two of its ends (bottom / top)
enum TIME_DOMAIN_TYPE {  BOTTOM = 0,
                         TOP = 1,
                         INTERVAL = 2
};


void IterateRange (int ne, LocalHeap & clh, const function<void(int,LocalHeap&)> & func);

const int SPACETIME_SANITY_CHECK_NR = -9;

INLINE void MarkAsSpaceTimeIntegrationPoint(ngcomp::IntegrationPoint & ip)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  ip.SetNr(SPACETIME_SANITY_CHECK_NR);
#endif  
}

INLINE void MarkAsSpaceTimeIntegrationRule(ngcomp::IntegrationRule & ir)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  for (auto & ip : ir)
    ip.SetNr(SPACETIME_SANITY_CHECK_NR);
#endif  
}

INLINE bool IsSpaceTimeIntegrationPoint(const ngcomp::IntegrationPoint & ip)
{
#ifdef SPACETIME_SANITY_CHECKS                    
  return (ip.Nr() == SPACETIME_SANITY_CHECK_NR);
#else
  return true;
#endif
}

ostream & operator<< (ostream & ost, DOMAIN_TYPE dt);
ostream & operator<< (ostream & ost, COMBINED_DOMAIN_TYPE cdt);

class GlobalNgsxfemVariables {
public:
    double EPS_STCR_LSET_PERTUBATION;
    double EPS_STCR_ROOT_SEARCH_BISECTION;
    double EPS_INTERPOLATE_TO_P1;
    double EPS_STFES_RESTRICT_GF;
    double EPS_SHIFTED_EVAL;
    double EPS_FACET_PATCH_INTEGRATOR;
    double MAX_DIST_NEWTON;

    int NEWTON_ITER_TRESHOLD;
    int FIXED_POINT_ITER_TRESHOLD;

    bool DO_NAIVE_TIMEINT;
    int NAIVE_TIMEINT_SUBDIVS;
    int NAIVE_TIMEINT_ORDER;

    int NON_CONV_WARN_MSG_LVL;

    bool SIMD_EVAL;

    GlobalNgsxfemVariables();
    
    void SetDefaults();
    void MultiplyAllEps(double fac);
    void SwitchSIMD(bool simd);
    void Output();
};

extern GlobalNgsxfemVariables globxvar;
