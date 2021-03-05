#pragma once
#include "../utils/ngsxstd.hpp"
#include <python_ngstd.hpp>

using namespace ngfem;
using namespace ngcomp;

namespace xintegration
{


  class LevelsetIntegrationDomain {
  private:
    Array<shared_ptr<GridFunction>> gfs_lset;
    Array<shared_ptr<CoefficientFunction>> cfs_lset;
    Array<Array<DOMAIN_TYPE>> dts;
    int intorder = -1;
    int force_intorder = -1;
    int time_intorder = -1;
    int subdivlvl = 0;
    SWAP_DIMENSIONS_POLICY quad_dir_policy = FIND_OPTIMAL;
    optional<double> tref;
  public:
    LevelsetIntegrationDomain( const Array<shared_ptr<CoefficientFunction>> & cfs_lset_in,
                               const Array<shared_ptr<GridFunction>> & gfs_lset_in,
                               const Array<Array<DOMAIN_TYPE>> & dts_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);
    LevelsetIntegrationDomain( const Array<shared_ptr<GridFunction>> & gfs_lset_in,
                               const Array<Array<DOMAIN_TYPE>> & dts_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);
    LevelsetIntegrationDomain( const Array<shared_ptr<CoefficientFunction>> & cfs_lset_in,
                               const Array<Array<DOMAIN_TYPE>> & dts_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);
    LevelsetIntegrationDomain( const shared_ptr<CoefficientFunction> & cf_lset_in,
                               const shared_ptr<GridFunction> & gf_lset_in,
                               DOMAIN_TYPE dt_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);
    LevelsetIntegrationDomain( const shared_ptr<GridFunction> & gf_lset_in,
                               DOMAIN_TYPE dt_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);
    LevelsetIntegrationDomain( const shared_ptr<CoefficientFunction> & cf_lset_in,
                               DOMAIN_TYPE dt_in,
                               int intorder_in = -1,
                               int time_intorder_in = -1,
                               int subdivlvl_in = 0,
                               SWAP_DIMENSIONS_POLICY quad_dir_policy_in = FIND_OPTIMAL,
                               optional<double> tref = nullopt);

    bool IsMultiLevelsetDomain () const
    {
      return ((gfs_lset.Size() > 1) || (dts.Size() > 1) || (dts[0].Size() > 1));
    }

    optional<double> OptionalReferenceTime() const {return tref;}
    bool HasReferenceTime() const { if (tref) return true; else return false; }
    double ReferenceTime() const 
    { 
      if (tref)
      {
        return *tref;
      }
      else
        throw Exception("no reference time stored.");
    }

    shared_ptr<CoefficientFunction> GetLevelsetCF () const
    {
      if (IsMultiLevelsetDomain ())
        throw Exception ("LevelsetIntegrationDomain is a MultiLevelsetDomain. ");
      else
        if (cfs_lset.Size() > 0)
          return cfs_lset[0];
        else
          return nullptr;
    }

    shared_ptr<GridFunction> GetLevelsetGF () const
    {
      if (IsMultiLevelsetDomain ())
        throw Exception ("LevelsetIntegrationDomain is a MultiLevelsetDomain. ");
      else
        if (gfs_lset.Size() > 0)
          return gfs_lset[0];
        else
          return nullptr;
    }

    DOMAIN_TYPE GetDomainType () const
    {
      if (IsMultiLevelsetDomain ())
        throw Exception ("LevelsetIntegrationDomain is a MultiLevelsetDomain. ");
      else
        if (dts.Size() > 0)
          return dts[0][0];
        else
          throw Exception ("dts empty.");
          
    }

    const Array<shared_ptr<CoefficientFunction>> & GetLevelsetCFs () const
    {
      return cfs_lset;
    }

    const Array<shared_ptr<GridFunction>> & GetLevelsetGFs () const
    {
      return gfs_lset;
    }

    const Array<Array<DOMAIN_TYPE>> & GetDomainTypes () const
    {
      return dts;
    }

    int GetIntegrationOrder () const
    {
      return intorder;
    }

    void SetIntegrationOrder (int order) 
    {
      intorder = order;
    }
    
    int GetTimeIntegrationOrder () const
    {
      return time_intorder;
    }

    void SetTimeIntegrationOrder (int tiorder) 
    {
      time_intorder = tiorder;
    }

    int GetNSubdivisionLevels () const
    {
      return subdivlvl;
    }

    SWAP_DIMENSIONS_POLICY GetSwapDimensionPolicy () const
    {
      return quad_dir_policy;
    }
    
  private:
  };

  ostream & operator<< (ostream & ost, const LevelsetIntegrationDomain & cdt);

  std::tuple<shared_ptr<CoefficientFunction>,shared_ptr<GridFunction>> CF2GFForStraightCutRule(shared_ptr<CoefficientFunction> cflset, int subdivlvl = 0);
  
  shared_ptr<LevelsetIntegrationDomain> PyDict2LevelsetIntegrationDomain(py::dict dictionary);

  
}

