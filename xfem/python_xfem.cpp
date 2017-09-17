#include <python_ngstd.hpp>
#include "../xfem/sFESpace.hpp"

//using namespace ngcomp;

void ExportNgsx_xfem(py::module &m)
{
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<SFESpace> PySFES;

  m.def("SFESpace", [](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
          -> PyFES
  {
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset, order, flags);
    LocalHeap lh (1000000, "SFESpace::Update-heap", true);
    ret->Update(lh);
    ret->FinalizeUpdate(lh);
    return ret;
  });
}

PYBIND11_PLUGIN(ngsxfem_xfem_py)
{
  cout << "importing ngsxfem-xfem lib" << endl;
  py::module m("ngsxfem-xfem", "pybind ngsxfem-xfem");
  ExportNgsx_xfem(m);
  return m.ptr();
}
