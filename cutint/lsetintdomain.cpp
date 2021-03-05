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

    optional<double> tref = nullopt;

    if (dictionary.contains("tref"))
      tref = py::cast<double>(dictionary["tref"]);

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
      return make_shared<LevelsetIntegrationDomain>(cf_lset,gf_lset,DOMAIN_TYPE(dt()),order,time_order,subdivlvl,quad_dir_policy,tref);
    }
    else
    {
      if (subdivlvl > 0)
        throw Exception("multi level sets don't work with subdivlvl > 0");
      if (time_order > -1)
        throw Exception("multi level sets don't work with time_order > 0");

      py::extract<py::tuple> lset_tuple_(lset);
      if (!lset_tuple_.check())
        throw Exception("lset is neither a level set nor a tuple ... need new candidates..");
      py::tuple lset_tuple(lset_tuple_());
      for (int i = 0; i < py::len(lset_tuple); i++)
        if (!(py::extract<shared_ptr<GridFunction>>(lset_tuple[i]).check()))
          throw Exception("lsets need to be GridFunctions!");
      
      Array<shared_ptr<GridFunction>> gf_lsets;
      gf_lsets = makeCArray<shared_ptr<GridFunction>> (lset_tuple);

      
      if (py::isinstance<py::tuple>(dt_in))
      {
        Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (py::extract<py::tuple>(dt_in)());
        Array<Array<DOMAIN_TYPE>> dts(1);
        dts[0] = dts_;
        return make_shared<LevelsetIntegrationDomain>(gf_lsets,dts,order,time_order,subdivlvl,quad_dir_policy,tref);
      }
      
      py::list dts_list;
      if (py::hasattr(dt_in, "as_list") && py::isinstance<py::list>(dt_in.attr("as_list")))
        dts_list = dt_in.attr("as_list");      
      else if (py::isinstance<py::list>(dt_in))
        dts_list = dt_in;
      else
        throw Exception("domain_type is neither a tuple nor a list nor a DomainTypeArray.");

      Array<Array<DOMAIN_TYPE>> dtas(py::len(dts_list));
      int common_length = -1; //not a list
      for (int i = 0; i < py::len(dts_list); i++)
      {
        auto dta(dts_list[i]);

        // Check for valid input
        if (!py::isinstance<py::tuple>(dta))
        {
          throw Exception("domain_type arrays are incompatible. Maybe you used a list instead of a tuple?");
        }
        else
        {
          if ((i>0) && (common_length != py::len(dta)))
            throw Exception("domain_type arrays have different length");
          else
            common_length = py::len(dta);
        }

        // Input valid. Pass input on 
        dtas[i] = makeCArray<DOMAIN_TYPE> (dta);
      }

      return make_shared<LevelsetIntegrationDomain>(gf_lsets,dtas,order,time_order,subdivlvl,quad_dir_policy,tref);
    }
  }
  
} // end of namespace
