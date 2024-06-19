#include "../utils/ngsxstd.hpp"

void IterateRange (int ne, LocalHeap & clh,
                   const function<void(int,LocalHeap&)> & func)
{
#ifndef WIN32
  if (task_manager)
  {
    SharedLoop2 sl(ne);
    task_manager -> CreateJob
      ( [&] (const TaskInfo & ti)
    {
      LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
      for (int elnr : sl)
      {
        HeapReset hr(lh);
        func (elnr,lh);
      }

    } );
  }
  else
#endif // WIN32
  {
    for (int elnr = 0; elnr < ne; elnr++)
    {
      HeapReset hr(clh);
      func (elnr,clh);
    }
  }
}

ostream & operator<< (ostream & ost, DOMAIN_TYPE dt)
{
  switch (dt)
    {
    case NEG: ost << "NEG"; break;
    case POS:  ost << "POS"; break;
    case IF:  ost << "IF"; break;
    };
  return ost;
}

ostream & operator<< (ostream & ost, COMBINED_DOMAIN_TYPE cdt)
{
  switch (cdt)
    {
    case CDOM_NO: ost << "NO"; break;
    case CDOM_NEG: ost << "NEG"; break;
    case CDOM_POS:  ost << "POS"; break;
    case CDOM_UNCUT:  ost << "UNCUT"; break;
    case CDOM_IF:  ost << "IF"; break;
    case CDOM_HASNEG:  ost << "HASNEG"; break;
    case CDOM_HASPOS:  ost << "HASPOS"; break;
    case CDOM_ANY:  ost << "ANY"; break;
    };
  return ost;
}

//GlobalNgsxfemVariables globxvar;

GlobalNgsxfemVariables::GlobalNgsxfemVariables() {
    SetDefaults();
}

void GlobalNgsxfemVariables::SetDefaults()  {
    EPS_STCR_LSET_PERTUBATION = 1e-14;
    EPS_STCR_ROOT_SEARCH_BISECTION = 1e-15;
    EPS_INTERPOLATE_TO_P1 = 1e-14;
    EPS_STFES_RESTRICT_GF = 1e-9;
    EPS_SHIFTED_EVAL = 1e-8;
    EPS_FACET_PATCH_INTEGRATOR = 1e-12;
    NEWTON_ITER_TRESHOLD = 15;
    MAX_DIST_NEWTON = 10;
    FIXED_POINT_ITER_TRESHOLD = 100;
    DO_NAIVE_TIMEINT = false;
    NAIVE_TIMEINT_ORDER = 3;
    NAIVE_TIMEINT_SUBDIVS = 10;
    NON_CONV_WARN_MSG_LVL = 3;
    SIMD_EVAL = true;
    cout << IM(3) << "All NGSXFEM eps values have been set to their default values" << endl;
}

void GlobalNgsxfemVariables::MultiplyAllEps(double fac) {
    EPS_STCR_LSET_PERTUBATION *= fac;
    EPS_STCR_ROOT_SEARCH_BISECTION *= fac;
    EPS_INTERPOLATE_TO_P1 *= fac;
    EPS_STFES_RESTRICT_GF *= fac;
    EPS_SHIFTED_EVAL *= fac;
    EPS_FACET_PATCH_INTEGRATOR *= fac;
}

void GlobalNgsxfemVariables::Output(){
    cout << "Report of GlobalNgsxfemVariables: " << endl;
    cout << "EPS_STCR_LSET_PERTUBATION = " << EPS_STCR_LSET_PERTUBATION << endl;
    cout << "EPS_STCR_ROOT_SEARCH_BISECTION = " << EPS_STCR_ROOT_SEARCH_BISECTION << endl;
    cout << "EPS_INTERPOLATE_TO_P1 = " << EPS_INTERPOLATE_TO_P1 << endl;
    cout << "EPS_STFES_RESTRICT_GF = " << EPS_STFES_RESTRICT_GF << endl;
    cout << "EPS_SHIFTED_EVAL = " << EPS_SHIFTED_EVAL << endl;
    cout << "EPS_FACET_PATCH_INTEGRATOR = " << EPS_FACET_PATCH_INTEGRATOR << endl;
    cout << "NEWTON_ITER_TRESHOLD = " << NEWTON_ITER_TRESHOLD << endl;
    cout << "MAX_DIST_NEWTON = " << MAX_DIST_NEWTON << endl;
    cout << "FIXED_POINT_ITER_TRESHOLD = " << FIXED_POINT_ITER_TRESHOLD << endl;
    cout << "DO_NAIVE_TIMEINT = " << DO_NAIVE_TIMEINT << endl;
    cout << "NAIVE_TIMEINT_ORDER = " << NAIVE_TIMEINT_ORDER << endl;
    cout << "NAIVE_TIMEINT_SUBDIVS = " << NAIVE_TIMEINT_SUBDIVS << endl;
    cout << "NON_CONV_WARN_MSG_LVL = " << NON_CONV_WARN_MSG_LVL << endl;
}
void GlobalNgsxfemVariables::SwitchSIMD(bool simd){
  SIMD_EVAL=simd;
}
