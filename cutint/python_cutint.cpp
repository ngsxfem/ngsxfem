#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../cutint/straightcutrule.hpp"
#include "../cutint/xintegration.hpp"
#include "../cutint/mlsetintegration.hpp"

using namespace xintegration;

typedef shared_ptr<CoefficientFunction> PyCF;
typedef GridFunction GF;
typedef shared_ptr<GF> PyGF;
  

shared_ptr<LevelsetIntegrationDomain> PyDict2LevelsetIntegrationDomain(py::dict dictionary)
{
  py::object lset = dictionary["levelset"];
  py::object dt_in = dictionary["domain_type"];
  auto subdivlvl_ = py::extract<int>(dictionary["subdivlvl"]);
  if (!subdivlvl_.check())
    throw Exception("data type for subdivlvl not admissible.");
  int subdivlvl = subdivlvl_();
  auto order_ = py::extract<int>(dictionary["order"]);
  if (!order_.check())
    throw Exception("data type for order not admissible.");
  int order = order_();
  auto time_order_ = py::extract<int>(dictionary["time_order"]);
  if (!time_order_.check())
    throw Exception("data type for time_order not admissible.");
  int time_order = time_order_();
  auto quad_dir_policy_ = py::extract<SWAP_DIMENSIONS_POLICY>(dictionary["quad_dir_policy"]);
  if (!quad_dir_policy_.check())
    throw Exception("data type for quad_dir_policy not admissible.");
  SWAP_DIMENSIONS_POLICY quad_dir_policy = quad_dir_policy_();

  if (py::extract<DOMAIN_TYPE> (dt_in).check())
  {
    cout << "not a list" << endl;
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
    cout << " a list" << endl;
    py::extract<py::list> dts_list(dt_in);
    if (!dts_list.check())
      throw Exception("domain_type is neither a DOMAIN_TYPE nor a list ... need new candidates..");
    //TODO: deal with Array<Array<DOMAIN_TYPE>> instead of Array<DOMAIN_TYPE> only
    Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (py::extract<py::list> (dt_in)());
    py::extract<py::list> lset_list(lset);
    if (!lset_list.check())
      throw Exception("lset is neither a level set nor a list ... need new candidates..");
    Array<shared_ptr<GridFunction>> gf_lsets;
    gf_lsets = makeCArray<shared_ptr<GridFunction>> (lset_list());
    //TODO: check if entries are GF or only CF
    Array<Array<DOMAIN_TYPE>> dts(1);
    dts[0] = dts_;
    if (subdivlvl > 0)
      throw Exception("multi level sets only work with grid functions and subdivlvl == 0.");
    return make_shared<LevelsetIntegrationDomain>(gf_lsets,dts,order,time_order,0,quad_dir_policy);
  }
}


void ExportNgsx_cutint(py::module &m)
{

//   m.def("IntegrateX_old",
//         [](py::object lset,
//            shared_ptr<MeshAccess> ma,
//            PyCF cf,
//            int order,
//            py::object dt_in,
//            int subdivlvl,
//            int time_order,
//            SWAP_DIMENSIONS_POLICY quad_dir_policy,
//            int heapsize)
//         {
//           static Timer t ("IntegrateX"); RegionTimer reg(t);
//           shared_ptr<LevelsetIntegrationDomain> lsetintdom = nullptr;
          
//           if (py::extract<DOMAIN_TYPE> (dt_in).check())
//           {
//             cout << "not a list" << endl;
//             py::extract<PyCF> pycf(lset);
//             py::extract<int> dt(dt_in);
//             if (!dt.check())
//               throw Exception("dt is not a domain type");
//             shared_ptr<GridFunction> gf_lset = nullptr;
//             shared_ptr<CoefficientFunction> cf_lset = nullptr;
//             tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(pycf(),subdivlvl);
//             lsetintdom = make_shared<LevelsetIntegrationDomain>(cf_lset,gf_lset,DOMAIN_TYPE(dt()),order,time_order,subdivlvl,quad_dir_policy);
//           }
//           else
//           {
//             cout << " a list" << endl;
//             py::extract<py::list> dts_list(dt_in);
//             if (!dts_list.check())
//               throw Exception("domain_type is neither a DOMAIN_TYPE nor a list ... need new candidates..");
//             //TODO: deal with Array<Array<DOMAIN_TYPE>> instead of Array<DOMAIN_TYPE> only
//             Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (py::extract<py::list> (dt_in)());
//             py::extract<py::list> lset_list(lset);
//             if (!lset_list.check())
//               throw Exception("lset is neither a level set nor a list ... need new candidates..");
//             Array<shared_ptr<GridFunction>> gf_lsets;
//             gf_lsets = makeCArray<shared_ptr<GridFunction>> (lset_list());
//             //TODO: check if entries are GF or only CF
//             Array<Array<DOMAIN_TYPE>> dts(1);
//             dts[0] = dts_;
//             if (subdivlvl > 0)
//               throw Exception("multi level sets only work with grid functions and subdivlvl == 0.");
//             lsetintdom = make_shared<LevelsetIntegrationDomain>(gf_lsets,dts,order,time_order,0,quad_dir_policy);
//           }
          
//           LocalHeap lh(heapsize, "lh-IntegrateX");

//           double sum = 0.0;
//           int DIM = ma->GetDimension();

//           Array<int> dnums;
//           ma->IterateElements
//             (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
//              {
//                auto & trafo = ma->GetTrafo (el, lh);

//                const IntegrationRule * ir;
//                Array<double> wei_arr;
//                tie (ir, wei_arr) = CreateCutIntegrationRule(*lsetintdom,trafo,lh);

//                if (ir != nullptr)
//                {
//                  BaseMappedIntegrationRule & mir = trafo(*ir, lh);
//                  FlatMatrix<> val(mir.Size(), 1, lh);

//                  cf -> Evaluate (mir, val);

//                  double lsum = 0.0;
//                  for (int i = 0; i < mir.Size(); i++)
//                      lsum += mir[i].GetMeasure()*wei_arr[i]*val(i,0);
//                  AtomicAdd(sum,lsum);
//                }
//              });

//           sum = ma->GetCommunicator().AllReduce(sum, MPI_SUM);
          
//           return sum;
//         },
//         py::arg("lset"),
//         py::arg("mesh"),
//         py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
//         py::arg("order")=5,
//         py::arg("domain_type")=IF,
//         py::arg("subdivlvl")=0,
//         py::arg("time_order")=-1,
//         py::arg("quad_dir_policy")=FIND_OPTIMAL,
//         py::arg("heapsize")=1000000,
//         docu_string(R"raw_string(
// Integrate on a level set domains. The accuracy of the integration is 'order' w.r.t. a (multi-)linear
// approximation of the level set function. At first, this implies that the accuracy will, in general,
// only be second order. However, if the isoparametric approach is used (cf. lsetcurving functionality)
// this will be improved.

// Parameters

// lset : ngsolve.CoefficientFunction or a list thereof
//   CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
//   FESpace with scalar continuous piecewise (multi-) linear basis functions.

// mesh : 
//   Mesh to integrate on (on some part) 

// cf : ngsolve.CoefficientFunction
//   the integrand

// order : int
//   integration order.

// domain_type : {NEG,POS,IF} (ENUM) or a list (of lists) thereof
//   Integration on the domain where either:
//   * the level set function is negative (NEG)
//   * the level set function is positive (POS)
//   * the level set function is zero     (IF )

// subdivlvl : int
//   On simplex meshes a subtriangulation is created on which the level set function lset is
//   interpolated piecewise linearly. Based on this approximation, the integration rule is
//   constructed. Note: this argument only works on simplices.

// time_order : int
//   integration order in time for space-time integration

// heapsize : int
//   heapsize for local computations.

// quad_dir_policy : int
//   policy for the selection of the order of integration directions
// )raw_string"));


    m.def("IntegrateX",
        [](py::dict lsetdom,
           shared_ptr<MeshAccess> ma,
           PyCF cf,
           int heapsize)
        {
          static Timer t ("IntegrateX"); RegionTimer reg(t);
          shared_ptr<LevelsetIntegrationDomain> lsetintdom = PyDict2LevelsetIntegrationDomain(lsetdom);
          
          LocalHeap lh(heapsize, "lh-IntegrateX");

          double sum = 0.0;
          int DIM = ma->GetDimension();

          Array<int> dnums;
          ma->IterateElements
            (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
             {
               auto & trafo = ma->GetTrafo (el, lh);

               const IntegrationRule * ir;
               Array<double> wei_arr;
               tie (ir, wei_arr) = CreateCutIntegrationRule(*lsetintdom,trafo,lh);

               if (ir != nullptr)
               {
                 BaseMappedIntegrationRule & mir = trafo(*ir, lh);
                 FlatMatrix<> val(mir.Size(), 1, lh);

                 cf -> Evaluate (mir, val);

                 double lsum = 0.0;
                 for (int i = 0; i < mir.Size(); i++)
                     lsum += mir[i].GetMeasure()*wei_arr[i]*val(i,0);
                 AtomicAdd(sum,lsum);
               }
             });

          sum = ma->GetCommunicator().AllReduce(sum, MPI_SUM);
          
          return sum;
        },
        py::arg("levelset_domain"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Integrate on a level set domains. The accuracy of the integration is 'order' w.r.t. a (multi-)linear
approximation of the level set function. At first, this implies that the accuracy will, in general,
only be second order. However, if the isoparametric approach is used (cf. lsetcurving functionality)
this will be improved.

Parameters

levelset_domain : dictionary which provides levelsets, domain_types and integration specifica. 
  [TODO: more detailed documentation (see old documentation)]

mesh : 
  Mesh to integrate on (on some part) 

cf : ngsolve.CoefficientFunction
  the integrand

order : int
  integration order.

heapsize : int
  heapsize for local computations.
)raw_string"));




  
}
