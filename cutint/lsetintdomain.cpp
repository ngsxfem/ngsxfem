#include "lsetintdomain.hpp"

namespace xintegration
{
  // using ngfem::INT;

  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const Array<shared_ptr<CoefficientFunction>> & cfs_lset_in,
                                                        const Array<shared_ptr<GridFunction>> & gfs_lset_in,
                                                        const Array<Array<DOMAIN_TYPE>> & dts_in,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(gfs_lset_in),
    cfs_lset(cfs_lset_in),
    dts(dts_in),
    intorder(intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    ;
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const Array<shared_ptr<GridFunction>> & gfs_lset_in,
                                                        const Array<Array<DOMAIN_TYPE>> & dts_in,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(gfs_lset_in),
    cfs_lset(0),
    dts(dts_in),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    ;
  }

  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const Array<shared_ptr<CoefficientFunction>> & cfs_lset_in,
                                                        const Array<Array<DOMAIN_TYPE>> & dts_in,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(0),
    cfs_lset(cfs_lset_in),
    dts(dts_in),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    ;
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const shared_ptr<CoefficientFunction> & cf_lset_in,
                                                        const shared_ptr<GridFunction> & gf_lset_in,
                                                        DOMAIN_TYPE dt,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(1),
    cfs_lset(1),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    gfs_lset[0] = gf_lset_in;    
    cfs_lset[0] = cf_lset_in;    
    dts[0].SetSize(1); dts[0][0] = dt;    
  }
  
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const shared_ptr<CoefficientFunction> & cf_lset_in,
                                                        DOMAIN_TYPE dt,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(0),
    cfs_lset(1),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    cfs_lset[0] = cf_lset_in;    
    dts[0].SetSize(1); dts[0][0] = dt;    
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const shared_ptr<GridFunction> & gf_lset_in,
                                                        DOMAIN_TYPE dt,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in)
  : gfs_lset(1),
    cfs_lset(0),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
  {
    gfs_lset[0] = gf_lset_in;    
    dts[0].SetSize(1); dts[0][0] = dt;    
  }
  
} // end of namespace
