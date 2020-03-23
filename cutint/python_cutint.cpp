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
  

void ExportNgsx_cutint(py::module &m)
{
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

levelset_domain : dictionary which provides levelsets, domain_types and integration specifica:
  important keys are "levelset", "domain_type", "order", the remainder are additional:

    "levelset" : ngsolve.CoefficientFunction or a list thereof
      CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
      FESpace with scalar continuous piecewise (multi-) linear basis functions.


    "order" : int
      integration order.

    "domain_type" : {NEG,POS,IF} (ENUM) or a list (of lists) thereof
      Integration on the domain where either:
      * the level set function is negative (NEG)
      * the level set function is positive (POS)
      * the level set function is zero     (IF )

    "subdivlvl" : int
      On simplex meshes a subtriangulation is created on which the level set function lset is
      interpolated piecewise linearly. Based on this approximation, the integration rule is
      constructed. Note: this argument only works on simplices without space-time and without 
      multiple levelsets.

    "time_order" : int
      integration order in time for space-time integration

    "quad_dir_policy" : int
      policy for the selection of the order of integration directions

mesh : 
  Mesh to integrate on (on some part) 

cf : ngsolve.CoefficientFunction
  the integrand

heapsize : int
  heapsize for local computations.
)raw_string"));




  
}
