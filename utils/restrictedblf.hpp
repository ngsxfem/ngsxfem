#pragma once
#include <comp.hpp>

namespace ngcomp
{

  template <class TM, class TV> 
  class  RestrictedBilinearForm : public T_BilinearForm<TM,TV>
  {
    shared_ptr<BitArray> el_restriction = nullptr;
    shared_ptr<BitArray> fac_restriction = nullptr;
  public:
    /// generate a bilinear-form
    // RestrictedBilinearForm () ;
    // /// generate a bilinear-form
    RestrictedBilinearForm (shared_ptr<FESpace> fespace,
                            const string & name,
                            shared_ptr<BitArray> ael_restriction,
                            shared_ptr<BitArray> afac_restriction,
                            const Flags & flags)
      : T_BilinearForm<TM,TV>(fespace, name, flags),
      el_restriction(ael_restriction),
      fac_restriction(afac_restriction)
    {
      ;
    }


    RestrictedBilinearForm (shared_ptr<FESpace> fespace,
                            shared_ptr<FESpace> fespace2,
                            const string & name,
                            shared_ptr<BitArray> ael_restriction,
                            shared_ptr<BitArray> afac_restriction,
                            const Flags & flags)
    : T_BilinearForm<TM,TV>(fespace, fespace2, name, flags),
      el_restriction(ael_restriction),
      fac_restriction(afac_restriction)
    {
     ;
    }

    void SetElementRestriction(shared_ptr<BitArray> _el_restriction){ el_restriction = _el_restriction;}
    void SetFacetRestriction(shared_ptr<BitArray> _fac_restriction){ fac_restriction = _fac_restriction;}
    shared_ptr<BitArray> GetElementRestriction(){ return el_restriction; }
    shared_ptr<BitArray> GetFacetRestriction(){ return fac_restriction; }

    /// generate a bilinear-form
    // BilinearForm (shared_ptr<FESpace> afespace, 
    //     	  shared_ptr<FESpace> afespace2, 
    //               shared_ptr<BitArray> el_restriction,
    //               shared_ptr<BitArray> fac_restriction,
    //     	  const string & aname,
    //     	  const Flags & flags);

    virtual MatrixGraph GetGraph (int level, bool symmetric); 

  };

}
