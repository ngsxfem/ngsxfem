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

               const IntegrationRule * ir = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, order, lh, subdivlvl);

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
        py::arg("heapsize")=1000000);

}

PYBIND11_PLUGIN(ngsxfem_cutint_py)
{
  cout << "importing ngsxfem-cutint lib" << endl;
  py::module m("ngsxfem-cutint", "pybind ngsxfem-cutint");
  ExportNgsx_cutint(m);
  return m.ptr();
}
