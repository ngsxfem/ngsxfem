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
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(gfs_lset_in),
    cfs_lset(cfs_lset_in),
    dts(dts_in),
    intorder(intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)
  {
    ;
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const Array<shared_ptr<GridFunction>> & gfs_lset_in,
                                                        const Array<Array<DOMAIN_TYPE>> & dts_in,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(gfs_lset_in),
    cfs_lset(0),
    dts(dts_in),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)
  {
    ;
  }

  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const Array<shared_ptr<CoefficientFunction>> & cfs_lset_in,
                                                        const Array<Array<DOMAIN_TYPE>> & dts_in,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(0),
    cfs_lset(cfs_lset_in),
    dts(dts_in),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)
  {
    ;
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const shared_ptr<CoefficientFunction> & cf_lset_in,
                                                        const shared_ptr<GridFunction> & gf_lset_in,
                                                        DOMAIN_TYPE dt,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(1),
    cfs_lset(1),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)
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
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(1),
    cfs_lset(1),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)
  {
    tie(cfs_lset[0],gfs_lset[0]) = CF2GFForStraightCutRule(cf_lset_in,subdivlvl);
    if (cfs_lset[0] == nullptr)
      cfs_lset.SetSize(0);
    if (gfs_lset[0] == nullptr)
      gfs_lset.SetSize(0);
    dts[0].SetSize(1); dts[0][0] = dt;    
  }
    
  LevelsetIntegrationDomain::LevelsetIntegrationDomain( const shared_ptr<GridFunction> & gf_lset_in,
                                                        DOMAIN_TYPE dt,
                                                        int intorder_in,
                                                        int time_intorder_in,
                                                        int subdivlvl_in,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy_in,
                                                        optional<double> tref_in)
  : gfs_lset(1),
    cfs_lset(0),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in),
    tref(tref_in)    
  {
    gfs_lset[0] = gf_lset_in;    
    dts[0].SetSize(1); dts[0][0] = dt;    
  }

  ostream & operator<< (ostream & ost, const LevelsetIntegrationDomain & lsetintdom)
  {
    if (lsetintdom.IsMultiLevelsetDomain())
    {
      ost << "MultiLevelsetDomain" << endl;
      ost << "GridFunctions: \n " << lsetintdom.GetLevelsetGFs() << endl;
      ost << "CoefficientFunctions: \n " << lsetintdom.GetLevelsetCFs() << endl;
      ost << "DomainTypes: \n " << lsetintdom.GetDomainTypes() << endl;
    }
    else
    {
      ost << "SingleLevelsetDomain" << endl;
      ost << "GridFunction: \n " << lsetintdom.GetLevelsetGF() << endl;
      ost << "CoefficientFunction: \n " << lsetintdom.GetLevelsetCF() << endl;
      ost << "DomainType: \n " << lsetintdom.GetDomainType() << endl;
    }
    ost << "IntegrationOrder: \n " << lsetintdom.GetIntegrationOrder() << endl;
    ost << "Time IntegrationOrder: \n " << lsetintdom.GetTimeIntegrationOrder() << endl;
    ost << "Number of subdivision levels: \n " << lsetintdom.GetNSubdivisionLevels() << endl;
    ost << "Policy on Quads/Hexes: \n " << lsetintdom.GetSwapDimensionPolicy() << endl;
    if (lsetintdom.HasReferenceTime())
      ost << "Fixed reference time : \n " << lsetintdom.ReferenceTime() << endl;
    return ost;
  }

  std::tuple<shared_ptr<CoefficientFunction>,shared_ptr<GridFunction>> CF2GFForStraightCutRule(shared_ptr<CoefficientFunction> cflset, int subdivlvl)
  {
    if (subdivlvl != 0)
      return make_tuple(cflset, nullptr);
    else
    {
      shared_ptr<GridFunction> ret = dynamic_pointer_cast<GridFunction>(cflset);
      if ((ret != nullptr) && (ret->GetFESpace()->GetOrder() <= 1) 
          && (    (ret->GetFESpace()->GetClassName() == "H1HighOrderFESpace")
               || (ret->GetFESpace()->GetClassName() == "DiscontinuousH1HighOrderFESpace")
               || (ret->GetFESpace()->GetClassName() == "SpaceTimeFESpace(H1HighOrderFESpace)")
             )
         )
        return make_tuple(nullptr, ret);
      else
        return make_tuple(cflset, nullptr);
    }
  }
  

} // end of namespace
