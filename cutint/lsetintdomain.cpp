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
  : gfs_lset(1),
    cfs_lset(1),
    dts(1),
    intorder(intorder_in),
    time_intorder(time_intorder_in),
    subdivlvl(subdivlvl_in),
    quad_dir_policy(quad_dir_policy_in)
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


  std::tuple<shared_ptr<CoefficientFunction>,shared_ptr<GridFunction>> CF2GFForStraightCutRule(shared_ptr<CoefficientFunction> cflset, int subdivlvl)
  {
    if (subdivlvl != 0)
      return make_tuple(cflset, nullptr);
    else
    {
      shared_ptr<GridFunction> ret = dynamic_pointer_cast<GridFunction>(cflset);
      if ((ret != nullptr) && (ret->GetFESpace()->GetOrder() <= 1) && ( (ret->GetFESpace()->GetClassName() == "H1HighOrderFESpace") || (ret->GetFESpace()->GetClassName() == "SpaceTimeFESpace")))
        return make_tuple(nullptr, ret);
      else
        return make_tuple(cflset, nullptr);
    }
  }
  

  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;

  shared_ptr<LevelsetIntegrationDomain> PyDict2LevelsetIntegrationDomain(py::dict dictionary)
  {
    if (!dictionary.contains("levelset"))
      throw Exception("You need to provide (a) levelset(s).");
    
    if (!dictionary.contains("domain_type"))
      throw Exception("You need to provide (a) domain type(s).");
      
    py::object lset = dictionary["levelset"];
    py::object dt_in = dictionary["domain_type"];

    int subdivlvl = 0;
    if (dictionary.contains("subdivlvl"))
    {
      auto subdivlvl_ = py::extract<int>(dictionary["subdivlvl"]);
      if (!subdivlvl_.check())
        throw Exception("data type for subdivlvl not admissible.");
      subdivlvl = subdivlvl_();
    }
    
    int order = -1;
    if (dictionary.contains("order"))
    {
      auto order_ = py::extract<int>(dictionary["order"]);
      if (!order_.check())
        throw Exception("data type for order not admissible.");
      order = order_();
    }
    
    int time_order = -1;
    if (dictionary.contains("time_order"))
    {
      auto time_order_ = py::extract<int>(dictionary["time_order"]);
      if (!time_order_.check())
        throw Exception("data type for time_order not admissible.");
      time_order= time_order_();
    }

    SWAP_DIMENSIONS_POLICY quad_dir_policy = FIND_OPTIMAL;
    if (dictionary.contains("quad_dir_policy"))
    {
      auto quad_dir_policy_ = py::extract<SWAP_DIMENSIONS_POLICY>(dictionary["quad_dir_policy"]);
      if (!quad_dir_policy_.check())
        throw Exception("data type for quad_dir_policy not admissible.");
      quad_dir_policy = quad_dir_policy_();
    }
    
    if (py::extract<DOMAIN_TYPE> (dt_in).check())
    {
      py::extract<PyCF> pycf(lset);
      py::extract<int> dt(dt_in);
      if (!dt.check())
        throw Exception("dt is not a domain type");
      shared_ptr<GridFunction> gf_lset = nullptr;
      shared_ptr<CoefficientFunction> cf_lset = nullptr;
      tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(pycf(),subdivlvl);
      return make_shared<LevelsetIntegrationDomain>(cf_lset,gf_lset,DOMAIN_TYPE(dt()),order,time_order,subdivlvl,quad_dir_policy);
    }
    else
    {
      if (subdivlvl > 0)
        throw Exception("multi level sets don't work with subdivlvl > 0");
      if (time_order > -1)
        throw Exception("multi level sets don't work with time_order > 0");
      
      py::extract<py::list> dts_list_(dt_in);
      if (!dts_list_.check())
        throw Exception("domain_type is neither a DOMAIN_TYPE nor a list ... need new candidates..");
      auto dts_list(dts_list_());

      int common_length = -1; //not a list
      for (int i = 0; i < py::len(dts_list); i++)
      {
        auto dts_list_entry(dts_list[i]);
        py::extract<py::list> dta(dts_list_entry);
        if (!dta.check())
        {
          if (common_length != -1)
            throw Exception("domain_type arrays are incompatible");
        }
        else
        {
          if ((i>0) && (common_length != py::len(dta())))
            throw Exception("domain_type arrays have different length");
          else
            common_length = py::len(dta());
        }
      }

      py::extract<py::list> lset_list(lset);
      if (!lset_list.check())
        throw Exception("lset is neither a level set nor a list ... need new candidates..");
      Array<shared_ptr<GridFunction>> gf_lsets;
      gf_lsets = makeCArray<shared_ptr<GridFunction>> (lset_list());
      
      if (common_length == -1) // not a list of lists
      {
        Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (dts_list);
        //TODO: check if entries are GF or only CF
        Array<Array<DOMAIN_TYPE>> dts(1);
        dts[0] = dts_;
        return make_shared<LevelsetIntegrationDomain>(gf_lsets,dts,order,time_order,subdivlvl,quad_dir_policy);
      }
      else
      {
        Array<Array<DOMAIN_TYPE>> dtas(py::len(dts_list));
        for (int i = 0; i < py::len(dts_list); i++)
        {
          py::extract<py::list> dta(dts_list[i]);
          dtas[i] = makeCArray<DOMAIN_TYPE> (dta());
        }
        return make_shared<LevelsetIntegrationDomain>(gf_lsets,dtas,order,time_order,subdivlvl,quad_dir_policy);
        
      }
    }
  }
  
} // end of namespace
