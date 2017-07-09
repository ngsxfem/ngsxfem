#pragma once
#include <comp.hpp>

namespace ngcomp
{

  class RestrictedBilinearForm : public T_BilinearForm<double,double>
  {
    shared_ptr<BitArray> el_restriction = nullptr;
    shared_ptr<BitArray> fac_restriction = nullptr;
  public:
    /// generate a bilinear-form
    // RestrictedBilinearForm () ;
    // /// generate a bilinear-form
    RestrictedBilinearForm (shared_ptr<FESpace> afespace,
                            const string & aname,
                            shared_ptr<BitArray> el_restriction,
                            shared_ptr<BitArray> fac_restriction,
                            const Flags & flags);

    /// generate a bilinear-form
    // BilinearForm (shared_ptr<FESpace> afespace, 
    //     	  shared_ptr<FESpace> afespace2, 
    //               shared_ptr<BitArray> el_restriction,
    //               shared_ptr<BitArray> fac_restriction,
    //     	  const string & aname,
    //     	  const Flags & flags);

    virtual MatrixGraph * GetGraph (int level, bool symmetric);
  };

}
