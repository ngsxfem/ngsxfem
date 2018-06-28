#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

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



  // py::class_<StatisticContainer, shared_ptr<StatisticContainer>>(m, "StatisticContainer")
  //   .def(py::init<>())
  //   .def("Print", [](StatisticContainer & self, string label, string select)
  //        {
  //          if (select == "L1")
  //            PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
  //          if (select == "L2")
  //            PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
  //          if (select == "max")
  //            PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
  //          if (select == "misc")
  //            PrintConvergenceTable(self.ErrorMisc,label+"_misc");
  //          if (select == "all")
  //          {
  //            PrintConvergenceTable(self.ErrorL1Norm,label+"_L1");
  //            PrintConvergenceTable(self.ErrorL2Norm,label+"_L2");
  //            PrintConvergenceTable(self.ErrorMaxNorm,label+"_max");
  //            PrintConvergenceTable(self.ErrorMisc,label+"_misc");
  //          }
  //        },
  //        py::arg("label")="something",py::arg("select")="all"
  //     )
  //   ;

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
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Compute approximated distance between of the isoparametrically obtained geometry.

  G_h = { phi_lin o Psi^{-1} }

and 

  G   = { phi = 0 }

as 

  max_{x in G_h} | phi(x) |

where 

  phi = lset_ho
  phi_lin = lset_p1
  Psi = Id + deform

The approximation is obtained as the maximum that is only computed on the integration points.

Parameters

lset_ho : ngsolve.CoefficientFunction
  level set (high order) function

lset_p1 : ngsolve.GridFunction
  P1 approximation of level set function

deform : ngsolve.GridFunction
  Deformation field describing a transformation

heapsize : int
  heapsize of local computations.
)raw_string")
    )
    ;

  
  m.def("CalcDistances",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
        {
          LocalHeap lh (heapsize, "CalcDistance-Heap");
          if (lset_p1->GetMeshAccess()->GetDimension()==2)
            CalcDistances<2>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
          else
            CalcDistances<3>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
        } ,
        py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("stats")=NULL,py::arg("heapsize")=1000000,py::arg("refine_threshold")=-1.0,py::arg("absolute")=false,
        docu_string(R"raw_string(
This is an internal function (and should be removed after some refactoring at some point)!
)raw_string")
    )
    ;

  // m.def("CalcDeformationError",  [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, StatisticContainer & stats, double lower, double upper, int heapsize)
  //       {
  //         LocalHeap lh (heapsize, "CalcDeformationError-Heap");
  //         if (lset_p1->GetMeshAccess()->GetDimension()==2)
  //           CalcDeformationError<2>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
  //         else
  //           CalcDeformationError<3>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
  //       } ,
  //       py::arg("lset_ho")=NULL,py::arg("lset_p1")=NULL,py::arg("deform")=NULL,py::arg("qn")=NULL,py::arg("stats")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000)
  //   ;

  m.def("ProjectShift",  [] (PyGF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn,
                             py::object active_elems_in,
                             PyCF blending,
                             double lower, double upper, double threshold, int heapsize)
        {
          shared_ptr<BitArray> active_elems = nullptr;
          if (py::extract<PyBA> (active_elems_in).check())
            active_elems = py::extract<PyBA>(active_elems_in)();
          
          LocalHeap lh (heapsize, "ProjectShift-Heap");
          ProjectShift(lset_ho, lset_p1, deform, qn, active_elems, blending, lower, upper, threshold, lh);
        } ,
        py::arg("lset_ho")=NULL,
        py::arg("lset_p1")=NULL,
        py::arg("deform")=NULL,
        py::arg("qn")=NULL,
        py::arg("active_elements")=DummyArgument(),
        py::arg("blending")=NULL,
        py::arg("lower")=0.0,
        py::arg("upper")=0.0,
        py::arg("threshold")=1.0,
        py::arg("heapsize")=1000000),
        docu_string(R"raw_string(
Computes the shift between points that are on the (P1 ) approximated level set function and its
higher order accurate version. This is only applied on elements where a level value inside
(lower,upper) exists. The result is put into deform (D) which is computed pointwise as

1)phi_lin( Psi(x) ) = phi_h(x)

  with Psi(x) = x + d(x) qn(x) =: x + D(x)

for all x on 'cut' elements

with

  phi_h : lset_ho
    the higher order level set function

  phi_lin : lset_p1
    the P1 level set function

  Psi : Id + deform
    the resulting deformation

  qn : normal direction field

Parameters

lset_ho : ngsolve.CoefficientFunction
  Scalar (higher order approximation) level set fct.

lset_p1 : ngsolve.GridFunction
  Scalar piecewise (multi-)linear Gridfunction

deform : ngsolve.GridFunction
  vector valued GridFunction to store the resulting deformation

active_elements : ngsolve.BitArray / None
  explicit marking of elements on which the transformation should be applied. If this is not None
  lower and upper will be ignored.

blending : ngsolve.CoefficientFunction
  Option to apply the mesh deformation more localized on cut elements. Setting blending function to
  0 (CoefficientFunction(0.0)) corresponds to applying the mapping on all points on cut elements
  completely. Using a blending function as a CoefficientFunction allows for a transition between the
  full application of the mapping (value 0) and no application of the mapping (value 1).

  This argument can be left away. Otherwise the mapping 1) is changed to

2)phi_lin(Psi(x))=phi_h(x)+b(x)·(phi_lin-phi_h)(x) 

  with a blending function b(x). Note that b(x) should be 0 where phi_lin(x) = 0

lower: float
  smallest relevant level set value to define the 'cut' elements where the mapping should be applied

upper: float
  highest relevant level set value to define the 'cut' elements where the mapping should be applied

threshold: float
  maximum (pointwise) value for d(x)/h in the mapping
    Psi(x) = x + d(x) qn(x)
  This might be necessary if the geometry is only coarsely approximated to avoid irregular meshes
  after a corresponding mesh deformation.

heapsize : int
  heapsize of local computations.
)raw_string")
    ;

// ProjectShift


  m.def("RefineAtLevelSet",  [] (PyGF lset_p1, double lower, double upper, int heapsize)
        {
          LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
          RefineAtLevelSet(lset_p1, lower, upper, lh);
        } ,
        py::arg("gf")=NULL,py::arg("lower")=0.0,py::arg("upper")=0.0,py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Mark mesh for refinement on all elements where the piecewise linear level set function lset_p1 has
values in the interval [lower,upper] (default [0,0]).

Parameters

gf : ngsolve.GridFunction
  Scalar piecewise (multi-)linear Gridfunction

lower : float
  smallest level set value of interest

upper : float
  largest level set value of interest

heapsize : int
  heapsize of local computations.
)raw_string"));

  m.def("shifted_eval", [](PyGF self,
                           py::object back_in,
                           py::object forth_in)
        -> PyCF
        {
          PyGF back = nullptr;
          if (py::extract<PyGF> (back_in).check())
            back = py::extract<PyGF>(back_in)();
          PyGF forth = nullptr;
          if (py::extract<PyGF> (forth_in).check())
            forth = py::extract<PyGF>(forth_in)();

          shared_ptr<DifferentialOperator> diffop  = nullptr;
          
          if (self->GetFESpace()->GetDimension() == 1)
            diffop = make_shared<DiffOpShiftedEval<1>> (back,forth);
          else if (self->GetFESpace()->GetDimension() == 2)
            diffop = make_shared<DiffOpShiftedEval<2>> (back,forth);
          else
            throw Exception("shifted_eval only for dim <= 2 so far");

          return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffop));
        },
        py::arg("gf"),
        py::arg("back") = DummyArgument(),
        py::arg("forth") = DummyArgument(),
        docu_string(R"raw_string(
Returns a CoefficientFunction that evaluates Gridfunction gf at a shifted location, s.t. the
original function to gf, gf: x -> f(x) is changed to cf: x -> f(s(x)) where z = s(x) is the shifted
location that is computed ( pointwise ) from:

     Psi_back(z) = Psi_forth(x),
< = >            z = Inv(Psi_back)( Psi_forth(x) )
< = >            s = Inv(Psi_back) o Psi_forth(x)

To compute z = s(x) a fixed point iteration is used.

ATTENTION: 
==========

If s(x) leaves the the element that the integration point x is defined on, it will *NOT* change the
element but result in an integration point that lies outside of the physical element.

Parameters

back : ngsolve.GridFunction
  transformation describing Psi_back as I + d_back where d_back is the deformation (can be None).

forth : ngsolve.GridFunction
  transformation describing Psi_forth as I + d_forth where d_forth is the deformation (can be None).

ASSUMPTIONS: 
============
- 2D mesh
- Gridfunction of dim=1 or dim=2 (ScalarFE behind it)
)raw_string"));
  
}

PYBIND11_MODULE(ngsxfem_lsetcurving_py,m)
{
  cout << "importing ngsxfem-lsetcurving lib" << endl;
  ExportNgsx_lsetcurving(m);
}
