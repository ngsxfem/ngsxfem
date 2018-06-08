#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../cutint/straightcutrule.hpp"
#include "../cutint/xintegration.hpp"

using namespace xintegration;

void ExportNgsx_cutint(py::module &m)
{

  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  
  m.def("IntegrateX",
        [](py::object lset,
           shared_ptr<MeshAccess> ma,
           PyCF cf,
           int order,
           DOMAIN_TYPE dt,
           int subdivlvl,
           int time_order,
           SWAP_DIMENSIONS_POLICY quad_dir_policy,
           int heapsize)
        {
          py::extract<PyCF> pycf(lset);
          if (!pycf.check())
            throw Exception("cast failed... need new candidates..");

          shared_ptr<GridFunction> gf_lset = nullptr;
          shared_ptr<CoefficientFunction> cf_lset = nullptr;
          tie(cf_lset,gf_lset) = CF2GFForStraightCutRule(pycf(),subdivlvl);

          LocalHeap lh(heapsize, "lh-IntegrateX");

          double sum = 0.0;
          int DIM = ma->GetDimension();

          Array<int> dnums;
          ma->IterateElements
            (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
             {
               auto & trafo = ma->GetTrafo (el, lh);

               const IntegrationRule * ir = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, order, time_order, lh, subdivlvl, quad_dir_policy);

               if (ir != nullptr)
               {
                 BaseMappedIntegrationRule & mir = trafo(*ir, lh);
                 FlatMatrix<> val(mir.Size(), 1, lh);

                 cf -> Evaluate (mir, val);

                 double lsum = 0.0;
                 for (int i = 0; i < mir.Size(); i++)
                   lsum += mir[i].GetWeight()*val(i,0);

                 AsAtomic(sum) += lsum;
               }
             });

          return sum;
        },
        py::arg("lset"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("order")=5,
        py::arg("domain_type")=IF,
        py::arg("subdivlvl")=0,
        py::arg("time_order")=-1,
        py::arg("quad_dir_policy")=FIND_OPTIMAL,
        py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Integrate on a level set domains. The accuracy of the integration is 'order' w.r.t. a (multi-)linear
approximation of the level set function. At first, this implies that the accuracy will, in general,
only be second order. However, if the isoparametric approach is used (cf. lsetcurving functionality)
this will be improved.

Parameters

lset : ngsolve.CoefficientFunction
  CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
  FESpace with scalar continuous piecewise (multi-) linear basis functions.

mesh : 
  Mesh to integrate on (on some part) 

cf : ngsolve.CoefficientFunction
  the integrand

order : int
  integration order.

domain_type : {NEG,POS,IF} (ENUM)
  Integration on the domain where either:
  * the level set function is negative (NEG)
  * the level set function is positive (POS)
  * the level set function is zero     (IF )

subdivlvl : int
  On simplex meshes a subtriangulation is created on which the level set function lset is
  interpolated piecewise linearly. Based on this approximation, the integration rule is
  constructed. Note: this argument only works on simplices.

time_order : int
  integration order in time for space-time integration

heapsize : int
  heapsize for local computations.

quad_dir_policy : int
  policy for the selection of the order of integration directions
)raw_string"));

}

PYBIND11_MODULE(ngsxfem_cutint_py,m)
{
  cout << "importing ngsxfem-cutint lib" << endl;
  ExportNgsx_cutint(m);
}
