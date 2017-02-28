#pragma once


/// domain types: two domains: POS/NEG and the diving interface IF
enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};


/// from ngsolve
#include <comp.hpp>

void IterateRange (int ne, LocalHeap & clh, const function<void(int,LocalHeap&)> & func);
