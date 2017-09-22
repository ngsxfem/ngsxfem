#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../lsetcurving/shiftedevaluate.hpp"

using namespace ngcomp;

void ExportNgsx_lsetcurving(py::module &m)
{
  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  typedef shared_ptr<BitArray> PyBA;


  m.def("InterpolateToP1",  [] (PyGF gf_ho, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(gf_ho, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,
        py::arg("eps_perturbation")=1e-16,py::arg("heapsize")=1000000)
    ;

  m.def("InterpolateToP1",  [] (PyCF coef, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(coef, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("coef"),py::arg("gf"),
        py::arg("eps_perturbation")=1e-16,py::arg("heapsize")=1000000)
    ;

  py::class_<StatisticContainer, shared_ptr<StatisticContainer>>(m, "StatisticContainer")
    .def(py::init<>())
    .def("Print", [](StatisticContainer & self, string label, string select)
         {
           if (select == "L1")
             PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
           if (select == "L2")
             PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
           if (select == "max")
             PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
           if (select == "misc")
             PrintConvergenceTable(self.ErrorMisc,label+"_misc");
           if (select == "all")
           {
             PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
             PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
             PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
             PrintConvergenceTable(self.ErrorMisc,label+"_misc");
           }
         },
         py::arg("label")="something",py::arg("select")="all"
      )
    ;

  m.def("CalcMaxDistance",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, int heapsize)
        {
          StatisticContainer dummy;
          LocalHeap lh (heapsize, "CalcDistance-Heap");
          if (lset_p1->GetMeshAccess()->GetDimension()==2)
            CalcDistances<2>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
          else
            CalcDistances<3>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
          return (double) dummy.ErrorMaxNorm[dummy.ErrorMaxNorm.Size()-1];
        } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("heapsize")=1000000)
    ;

  m.def("CalcDistances",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
        {
          LocalHeap lh (heapsize, "CalcDistance-Heap");
          if (lset_p1->GetMeshAccess()->GetDimension()==2)
            CalcDistances<2>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
          else
            CalcDistances<3>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
        } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("stats")=NULL,py::arg("heapsize")=1000000,py::arg("refine_threshold")=-1.0,py::arg("absolute")=false)
    ;

  m.def("CalcDeformationError",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, StatisticContainer & stats, double lower, double upper, int heapsize)
        {
          LocalHeap lh (heapsize, "CalcDeformationError-Heap");
          if (lset_p1->GetMeshAccess()->GetDimension()==2)
            CalcDeformationError<2>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
          else
            CalcDeformationError<3>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
        } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL,py::arg("stats")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
    ;

  m.def("ProjectShift",  [] (PyGF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, PyCF blending,
                             double lower, double upper, double threshold, int heapsize)
        {
          LocalHeap lh (heapsize, "ProjectShift-Heap");
          ProjectShift(lset_ho, lset_p1, deform, qn, blending, lower, upper, threshold, lh);
        } ,
        py::arg("lset_ho")=NULL,
        py::arg("lset_p1")=NULL,
        py::arg("deform")=NULL,
        py::arg("qn")=NULL,
        py::arg("blending")=NULL,
        py::arg("lower")=0.0,
        py::arg("upper")=0.0,
        py::arg("threshold")=1.0,
        py::arg("heapsize")=1000000)
    ;

// ProjectShift


  m.def("RefineAtLevelSet",  [] (PyGF lset_p1, double lower, double upper, int heapsize)
        {
          LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
          RefineAtLevelSet(lset_p1, lower, upper, lh);
        } ,
        py::arg("lset_p1")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
    ;



  m.def("RefineAtLevelSet",  [] (PyGF gf, double lower_lset_bound, double upper_lset_bound, int heapsize)
        {
          LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
          RefineAtLevelSet(gf,lower_lset_bound,upper_lset_bound,lh);
        } ,
        py::arg("gf"),py::arg("lower_lset_bound")=0.0,py::arg("upper_lset_bound")=0.0,py::arg("heapsize")=10000000)
    ;

  m.def("shifted_eval", [](PyGF self, PyGF back,PyGF forth ) -> PyCF
        {

          auto diffop = make_shared<DiffOpShiftedEval> (back,forth);


          return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffop));
        });


  
}

PYBIND11_PLUGIN(ngsxfem_lsetcurving_py)
{
  cout << "importing ngsxfem-lsetcurving lib" << endl;
  py::module m("ngsxfem-lsetcurving", "pybind ngsxfem-lsetcurving");
  ExportNgsx_lsetcurving(m);
  return m.ptr();
}
