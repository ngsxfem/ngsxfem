#pragma once

/// from ngsolve
#include <comp.hpp>


/// domain types: two domains: POS/NEG and the diving interface IF
enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};

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
