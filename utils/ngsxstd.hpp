#pragma once

/// from ngsolve
#include <comp.hpp>

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


void IterateRange (int ne, LocalHeap & clh, const function<void(int,LocalHeap&)> & func);
