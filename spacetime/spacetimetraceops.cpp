
#include "spacetimefespace.hpp"
// #include "spacetimefe.hpp"

namespace ngcomp
{

class TraceEvaluator
{
protected:
  const GridFunction & u_st;
  GridFunction & u_tr;
public:
  TraceEvaluator(shared_ptr<GridFunction> a_u_st, shared_ptr<GridFunction> a_u_trace)
    : u_st(a_u_st), u_tr(a_u_trace)
  {
    ;
  }

  void EvaluateTrace(double time, LocalHeap & lh)
  {
    // apply trace operation
    BaseVector & v_st = u_st->GetVector();
    BaseVector & v_tr = u_tr->GetVector();
    
    const int ndof_st = v_st.Size();
    const int ndof_tr = v_tr.Size();

    shared_ptr<SpaceTimeFESpace> fes_st = dynamic_pointer_cast<SpaceTimeFESpace> (u_st->GetFESpace());
    shared_ptr<FESpace> fes_tr = u_tr->GetFESpace();
    const int time_order = fes_st->OrderTime();

    if (ndof_st != ndof_tr * (time_order+1))
      throw Exception ("SpaceTime and Trace Gridfunction are not / no longer compatible in size!");

    DGFiniteElement<1> & fel_time = dynamic_cast<const DGFiniteElement<1> &> (fes_st->GetTimeFE(lh));

    // shape_time
    // -> fel_time.CalcShape(time, ..) 

    for (int i = 0; i < ndof_tr; ++i)
    {
      for (int j = 0; j < ndof_time; ++j)
      {
        // ....
      }
      
    }
    
  }


};



}
