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

eps_collection_class eps_collection;

eps_collection_class::eps_collection_class() {
    cout << IM(4) << "Hello from the eps_collection_class::eps_collection_class constructor" << endl;
    SetAll(1e-14);
}

void eps_collection_class::SetAll(double val)  {
        EPS_STCR_LSET_PERTUBATION = val;
        EPS_STCR_ROOT_SEARCH_BISECTION = val;
        EPS_INTERPOLATE_TO_P1 = val;
        EPS_STFES_RESTRICT_GF = val;
        EPS_SHIFTED_EVAL = val;
        EPS_FACET_PATCH_INTEGRATOR = val;
        cout << IM(3) << "All NGSXFEM eps values have been set to " << val << endl;
}
