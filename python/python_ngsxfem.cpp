
#include <regex>
#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../cutint/python_cutint.cpp"
#include "../utils/python_utils.cpp"
#include "../xfem/python_xfem.cpp"

//using namespace ngcomp;
#include "../utils/ngsxstd.hpp"
#include "../cutint/straightcutrule.hpp" // ?
#include "../cutint/xintegration.hpp"    // ?
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../lsetcurving/shiftedevaluate.hpp"

#include "../cutint/spacetimecutrule.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../spacetime/SpaceTimeFESpace.hpp"
#include "../spacetime/diffopDt.hpp"
#include "../spacetime/timecf.hpp"

#include "../utils/bitarraycf.hpp"
#include "../utils/restrictedblf.hpp"
#include "../utils/p1interpol.hpp"
#include "../utils/xprolongation.hpp"

#include "../xfem/sFESpace.hpp"
#include "../xfem/cutinfo.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/ghostpenalty.hpp"


using namespace xintegration;
using namespace ngcomp;

void ExportNgsx(py::module &m)
{

  py::enum_<DOMAIN_TYPE>(m, "DOMAIN_TYPE")
    .value("POS", POS)
    .value("NEG", NEG)
    .value("IF", IF)
    .export_values()
    ;

  py::enum_<COMBINED_DOMAIN_TYPE>(m, "COMBINED_DOMAIN_TYPE")
    .value("NO", CDOM_NO)
    .value("CDOM_NEG", CDOM_NEG)
    .value("CDOM_POS", CDOM_POS)
    .value("UNCUT", CDOM_UNCUT)
    .value("CDOM_IF", CDOM_IF)
    .value("HASNEG", CDOM_HASNEG)
    .value("HASPOS", CDOM_HASPOS)
    .value("ANY", CDOM_ANY)
    .export_values()
    ;
  
  py::enum_<SWAP_DIMENSIONS_POLICY>(m, "QUAD_DIRECTION_POLICY")
    .value("FIRST", FIRST_ALLOWED)
    .value("OPTIMAL", FIND_OPTIMAL)
    .value("FALLBACK", ALWAYS_NONE)
    .export_values()
    ;


  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;



  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<BitArray> PyBA;



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
          
          if (self->GetFESpace()->GetDimension() == 1) {
              if (self->GetFESpace()->GetSpatialDimension() == 1)
                diffop = make_shared<DiffOpShiftedEval<1,1>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 2)
                diffop = make_shared<DiffOpShiftedEval<1,2>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 3)
                diffop = make_shared<DiffOpShiftedEval<1,3>> (back,forth);
              else throw Exception("shifted_eval only for space dim = 1,2,3 so far");
          }
          else if (self->GetFESpace()->GetDimension() == 2) {
              if (self->GetFESpace()->GetSpatialDimension() == 1)
                diffop = make_shared<DiffOpShiftedEval<2,1>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 2)
                diffop = make_shared<DiffOpShiftedEval<2,2>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 3)
                diffop = make_shared<DiffOpShiftedEval<2,3>> (back,forth);
              else throw Exception("shifted_eval only for space dim = 1,2,3 so far");
          }
          else if (self->GetFESpace()->GetDimension() == 3) {
              if (self->GetFESpace()->GetSpatialDimension() == 1)
                diffop = make_shared<DiffOpShiftedEval<3,1>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 2)
                diffop = make_shared<DiffOpShiftedEval<3,2>> (back,forth);
              else if (self->GetFESpace()->GetSpatialDimension() == 3)
                diffop = make_shared<DiffOpShiftedEval<3,3>> (back,forth);
              else throw Exception("shifted_eval only for space dim = 1,2,3 so far");
          }
          else
            throw Exception("shifted_eval only for dim <= 3 so far");

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

  typedef shared_ptr<SpaceTimeFESpace> PySTFES;
  typedef shared_ptr<ProxyFunction> PyProxyFunction;

  m.def("SpaceTimeFESpace", [] (
                                        PyFES basefes,
                                        shared_ptr<FiniteElement> fe,
                                        py::object dirichlet,
                                        py::dict bpflags,
                                        int heapsize)
  {


    shared_ptr<SpaceTimeFESpace> ret = nullptr;
    Flags flags = py::extract<Flags> (bpflags)();
    shared_ptr<MeshAccess> ma = basefes->GetMeshAccess();

    if (py::isinstance<py::list>(dirichlet)) {
        flags.SetFlag("dirichlet", makeCArray<double>(py::list(dirichlet)));

    }

    if (py::isinstance<py::str>(dirichlet))
    {
          std::regex pattern(dirichlet.cast<string>());
          Array<double> dirlist;
          for (int i = 0; i < ma->GetNBoundaries(); i++)
             if (std::regex_match (ma->GetMaterial(BND, i), pattern))
               dirlist.Append (i+1);
               flags.SetFlag("dirichlet", dirlist);
    }

    auto tfe = dynamic_pointer_cast<ScalarFiniteElement<1>>(fe);
    //cout << tfe << endl;
    if(tfe == nullptr)
      cout << "Warning! tfe == nullptr" << endl;

    ret = make_shared<SpaceTimeFESpace> (ma, basefes,tfe, flags);



    LocalHeap lh (heapsize, "SpaceTimeFESpace::Update-heap", true);
    ret->Update(lh);
    ret->FinalizeUpdate(lh);
    return ret;
  },
       py::arg("spacefes"),
       py::arg("timefe"),
       py::arg("dirichlet") = DummyArgument(),
       py::arg("flags") = py::dict(),
       py::arg("heapsize") = 1000000,
       docu_string(R"raw_string(
This function creates a SpaceTimeFiniteElementSpace based on a spacial FE space and a time Finite element
Roughly, this is the tensor product between those two arguments. Further arguments specify several details.

Parameters

spacefes : ngsolve.FESpace
  This is the spacial finite element used for the space-time discretisation.
  Both scalar and vector valued spaces might be used. An example would be
  spacefes = H1(mesh, order=order) for given mesh and order.

timefe : ngsolve.FiniteElement
  This is the time finite element for the space-time discretisation. That is
  essentially a simple finite element on the unit interval. There is a class
  ScalarTimeFE to create something fitting here. For example, one could call
  timefe = ScalarTimeFE(order) to create a time finite element of order order.

dirichlet : list or string
  The boundary of the space domain which should have Dirichlet boundary values.
  Specification policy is the same as with the usual space finite element spaces.

flags : dict
  Dictionary of additional flags for the finite element space. In the constructor
  of the SpaceTimeFESpace, those will be forwarded to the contructor of the general
  base class FESpace. An example would be flags = {"dgjumps": True}.

heapsize : int
  Size of the local heap of this class. Increase this if you observe errors which look
  like a heap overflow.

       )raw_string")
               );

  py::class_<SpaceTimeFESpace, PySTFES, FESpace>
    (m, "CSpaceTimeFESpace")
  .def("SetTime", [](PySTFES self, double t)
  {
    self->SetTime(t);
  },
       "Set the time variable\n Also sets override time")
  .def("SetOverrideTime", [](PySTFES self, bool override)
  {
    self->SetOverrideTime(override);
  },
       "Set flag to or not to override the time variable")
  .def("k_t", [](PySTFES self)
  {
     return self->order_time();
  },
     "Return order of the time FE")
  .def("TimeFE_nodes", [](PySTFES self)
  {
      Array<double> & nodesr = self->TimeFE_nodes();
      py::list nodes (nodesr.Size());
      for (int i = 0; i < nodesr.Size(); i++)
        nodes[i] = nodesr[i];
      return nodes;
   },
     "Return nodes of time FE")
  .def("IsTimeNodeActive", [](PySTFES self, int i)
  {
      return self->IsTimeNodeActive(i);
   },
     "Return bool whether node is active")
  ;

  m.def("ScalarTimeFE", []( int order, bool skip_first_node, bool only_first_node)
  {
    BaseScalarFiniteElement * fe = nullptr;

    if (skip_first_node && only_first_node)
      throw Exception("can't skip and keep first node at the same time.");
    fe = new NodalTimeFE(order, skip_first_node, only_first_node);
    return shared_ptr<BaseScalarFiniteElement>(fe);
  },
  py::arg("order") = 0,
  py::arg("skip_first_node") = false,
  py::arg("only_first_node") = false,
  docu_string(R"raw_string(
Creates a nodal Finite element in time on the interval [0,1].
Internally, Gauss-Lobatto integration points are exploited for that.

Parameters

order : int
The polynomial order of the discretisation. That controlls the number of
points in the time interval. See Gauss-Lobatto points for further details.
Currently, orders up to 5 are available.

skip_first_node : bool
This will create the time finite element without the first node at t=0.
That feature comes in handy for several CG like implementations in time.
Also see only_first_node.

only_first_node : bool
This will create the time finite element with only the first node at t=0.
That feature comes in handy for several CG like implementations in time.
Also see skip_first_node.
  )raw_string")
   );


  // DiffOpDt

  m.def("dt", [] (const PyProxyFunction self,py::object comp)
  {
    Array<int> comparr(0);
    if (py::extract<int> (comp).check())
    {
      int c = py::extract<int>(comp)();
      if (c != -1)
      {
        comparr.SetSize(1);
        comparr[0] = c;
      }
    }

    if (py::extract<py::list> (comp).check())
      comparr = makeCArray<int> (py::extract<py::list> (comp)());

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopdt;
    if ( self->GetFESpace()->GetSpatialDimension() == 2) {
        diffopdt = make_shared<T_DifferentialOperator<DiffOpDt<2>>> ();
    }
    else if( self->GetFESpace()->GetSpatialDimension() == 3) {
        diffopdt = make_shared<T_DifferentialOperator<DiffOpDt<3>>> ();
    }
    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopdt = make_shared<CompoundDifferentialOperator> (diffopdt, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->GetFESpace(), self->IsTestFunction(), self->IsComplex(),diffopdt, nullptr, nullptr, nullptr, nullptr, nullptr);
    
    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
    },
          py::arg("proxy"),
          py::arg("comp") = -1,
        docu_string(R"raw_string(
dt is the differential operator in time. This is the variant for a proxy function.

Parameters

proxy : ngsolve.ProxyFunction
  Function to differentiate
  
comp : int or list
  ??
  
)raw_string")
          );

  m.def("dt", [](PyGF self) -> PyCF
  {
    shared_ptr<DifferentialOperator> diffopdt;
    if ( self->GetFESpace()->GetSpatialDimension() == 2)
        diffopdt = make_shared<T_DifferentialOperator<DiffOpDt<2>>> ();
    else if ( self->GetFESpace()->GetSpatialDimension() == 3)
        diffopdt = make_shared<T_DifferentialOperator<DiffOpDt<3>>> ();

    return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopdt));
  }, docu_string(R"raw_string(
dt is the differential operator in time. For a given GridFunction gfu,
dt (gfu) will be its time derivative

Parameters

self : ngsolve.GridFunction
  Function to differentiate

)raw_string")
);

  typedef shared_ptr<TimeVariableCoefficientFunction> PyTimeVariableCF;

  py::class_<TimeVariableCoefficientFunction, PyTimeVariableCF, CoefficientFunction>(m, "TimeVariableCoefficientFunction")
          .def("__init__", [] () -> PyTimeVariableCF { return make_shared<TimeVariableCoefficientFunction>(); })
          .def("FixTime", &TimeVariableCoefficientFunction::FixTime)
          .def("UnfixTime", &TimeVariableCoefficientFunction::UnfixTime)
          .def("IsFixed", [] (PyTimeVariableCF self) -> bool {
            try {
                self->Evaluate(BaseMappedIntegrationPoint());
            } catch (...) {
                 return false;
            }
            return true;
          });

   m.def("ReferenceTimeVariable", []() -> PyTimeVariableCF
   {
     return make_shared<TimeVariableCoefficientFunction> ();
   }, docu_string(R"raw_string(
This is the time variable. Call tref = ReferenceTimeVariable() to have a symbolic variable
for the time like x,y,z for space. That can be used e.g. in lset functions for unfitted methods.
Note that one would typically use tref in [0,1] as one time slab, leading to a call like
t = told + delta_t * tref, when tref is our ReferenceTimeVariable.
)raw_string")
);

   // DiffOpDtVec

   m.def("dt_vec", [] (const PyProxyFunction self,py::object comp)
   {
     Array<int> comparr(0);
     if (py::extract<int> (comp).check())
     {
       int c = py::extract<int>(comp)();
       if (c != -1)
       {
         comparr.SetSize(1);
         comparr[0] = c;
       }
     }

     if (py::extract<py::list> (comp).check())
       comparr = makeCArray<int> (py::extract<py::list> (comp)());

     if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
     {
       throw Exception("cannot work with compounddiffops, prescribe comp != -1");
     }

     shared_ptr<DifferentialOperator> diffopdtvec;
     const int SpaceD = self->GetFESpace()->GetSpatialDimension();
     switch (self->Dimension())
     {
       case 1 : {
         if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 1>>> ();
         else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 1>>> ();
         break;
     }
       case 2 : {
         if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 2>>> ();
         else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 2>>> ();
         break;
     }
       case 3 : {
         if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 3>>> ();
         else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 3>>> ();
         break;
     }
       default : throw Exception("Diffop dt only implemented for dim <= 3 so far.");
     }

     for (int i = comparr.Size() - 1; i >= 0; --i)
     {
       diffopdtvec = make_shared<CompoundDifferentialOperator> (diffopdtvec, comparr[i]);
     }

     auto adddiffop = make_shared<ProxyFunction> (self->GetFESpace(), self->IsTestFunction(), self->IsComplex(), diffopdtvec, nullptr, nullptr, nullptr, nullptr, nullptr);

     if (self->IsOther())
       adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

     return PyProxyFunction(adddiffop);
     },
           py::arg("proxy"),
           py::arg("comp") = -1
           );

   m.def("dt_vec", [](PyGF self) -> PyCF
   {
       shared_ptr<DifferentialOperator> diffopdtvec;
       const int SpaceD = self->GetFESpace()->GetSpatialDimension();
       switch (self->Dimension())
       {
         case 1 : {
           if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 1>>> ();
           else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 1>>> ();
           break;
       }
         case 2 : {
           if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 2>>> ();
           else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 2>>> ();
           break;
       }
         case 3 : {
           if(SpaceD == 2) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<2, 3>>> ();
           else if(SpaceD == 3) diffopdtvec = make_shared<T_DifferentialOperator<DiffOpDtVec<3, 3>>> ();
           break;
       }
         default : throw Exception("Diffop dt only implemented for dim <= 3 so far.");
       }

     return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopdtvec,nullptr,nullptr,0));
   });

   // DiffOpFixt

  m.def("fix_t", [] (const PyProxyFunction self, double time, py::object comp, bool use_FixAnyTime )
  {
    Array<int> comparr(0);
    if (py::extract<int> (comp).check())
    {
      int c = py::extract<int>(comp)();
      if (c != -1)
      {
        comparr.SetSize(1);
        comparr[0] = c;
      }
    }

    if (py::extract<py::list> (comp).check())
      comparr = makeCArray<int> (py::extract<py::list> (comp)());

    if (comparr.Size()== 0 && dynamic_pointer_cast<CompoundDifferentialOperator>(self->Evaluator()))
    {
      throw Exception("cannot work with compounddiffops, prescribe comp != -1");
    }

    shared_ptr<DifferentialOperator> diffopfixt;
    const int SpaceD = self->GetFESpace()->GetSpatialDimension();
    if(!use_FixAnyTime && (time == 0.0 || time == 1.0))
    {
      switch (int(time))
      {
        case 0 : {
          if(SpaceD == 2) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, 0>>> ();
          else if(SpaceD == 3) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, 0>>> ();
          break;
      }
        case 1 : {
          if(SpaceD == 2) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, 1>>> ();
          else if(SpaceD == 3) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, 1>>> ();
          break;
      }
        default : throw Exception("Requested time not implemented yet.");
      }
    }
    else {
      cout << "Calling DiffOpFixAnyTime" << endl;
      if(SpaceD == 2) diffopfixt = make_shared<DiffOpFixAnyTime<2>> (time);
      else if(SpaceD == 3) diffopfixt = make_shared<DiffOpFixAnyTime<3>> (time);
    }


    for (int i = comparr.Size() - 1; i >= 0; --i)
    {
      diffopfixt = make_shared<CompoundDifferentialOperator> (diffopfixt, comparr[i]);
    }

    auto adddiffop = make_shared<ProxyFunction> (self->GetFESpace(), self->IsTestFunction(), self->IsComplex(), diffopfixt, nullptr, nullptr, nullptr, nullptr, nullptr);

    if (self->IsOther())
      adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

    return PyProxyFunction(adddiffop);
    },
          py::arg("proxy"),
          py::arg("time"),
          py::arg("comp") = -1,
          py::arg("use_FixAnyTime") = false
          );

   m.def("fix_t", [](PyGF self, double time, bool use_FixAnyTime) -> PyCF
   {
     shared_ptr<DifferentialOperator> diffopfixt;
     const int SpaceD = self->GetFESpace()->GetSpatialDimension();
     if(!use_FixAnyTime && (time == 0.0 || time == 1.0))
     {
       switch (int(time))
       {
         case 0 : {
           if(SpaceD == 2) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, 0>>> ();
           else if(SpaceD == 3) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, 0>>> ();
           break;
       }
         case 1 : {
           if(SpaceD == 2) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<2, 1>>> ();
           else if(SpaceD == 3) diffopfixt = make_shared<T_DifferentialOperator<DiffOpFixt<3, 1>>> ();
           break;
       }
         default : throw Exception("Requested time not implemented yet.");
       }
     }
     else {
       cout << "Calling DiffOpFixAnyTime" << endl;
       if(SpaceD == 2) diffopfixt = make_shared<DiffOpFixAnyTime<2>> (time);
       else if(SpaceD == 3) diffopfixt = make_shared<DiffOpFixAnyTime<3>> (time);
     }


     return PyCF(make_shared<GridFunctionCoefficientFunction> (self, diffopfixt));
   },
          py::arg("gridfunction"),
          py::arg("time"),
          py::arg("use_FixAnyTime") = false
          ,
    docu_string(R"raw_string(
fix_t fixes the time (ReferenceTimeVariable) of a given expression.
This is the variant for a gridfunction.

Parameters

self: ngsolve.GridFunction
  Gridfunction in which the time should be fixed
  
time: double
  Value the time should become
  
use_FixAnyTime: bool
  Bool flag to control whether the time value should be expected to be
  a node of the time scalar finite element or interpolation should be used.
  use_FixAnyTime = True means interpolation is used. use_FixAnyTime = False
  currently only supports time 0 and 1.

)raw_string")
);

   m.def("CreateTimeRestrictedGF", [](PyGF st_GF,double time) -> PyGF
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     return st_FES->CreateRestrictedGF(st_GF,time);
   },
   py::arg("gf"),
   py::arg("reference_time") = 0.0,
   "Create spatial-only Gridfunction corresponding to a fixed time.");

   m.def("RestrictGFInTime", [](PyGF st_GF,double time,PyGF s_GF)
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     switch (st_FES->GetDimension())
     {
       case 1:
         st_FES->RestrictGFInTime<double>(st_GF,time,s_GF);
         break;
       case 2:
         st_FES->RestrictGFInTime<Vec<2>>(st_GF,time,s_GF);
         break;
       case 3:
         st_FES->RestrictGFInTime<Vec<3>>(st_GF,time,s_GF);
         break;
       default:
         throw Exception("cannot handle GridFunction type (dimension too large?).");
         break;
     }
   }, 
   py::arg("spacetime_gf"),
   py::arg("reference_time") = 0.0,
   py::arg("space_gf"),
   "Extract Gridfunction corresponding to a fixed time from a space-time GridFunction.");

   m.def("SpaceTimeInterpolateToP1", [](PyCF st_CF, PyCF tref, PyGF st_GF)
   {
     FESpace* raw_FE = (st_GF->GetFESpace()).get();
     SpaceTimeFESpace * st_FES = dynamic_cast<SpaceTimeFESpace*>(raw_FE);
     if (!st_FES) throw Exception("not a spacetime gridfunction");
     st_FES->InterpolateToP1(st_CF,tref,st_GF);
   }, 
   py::arg("spacetime_cf"),
   py::arg("time"),
   py::arg("spacetime_gf"),
   "Interpolate nodal in time (possible high order) and nodal in space (P1).");


}

PYBIND11_MODULE(ngsxfem_py, m)
{
  cout << "importing ngs-xfem" << NGSXFEM_VERSION << endl;
  ExportNgsx(m);
  ExportNgsx_cutint(m);
  ExportNgsx_utils(m);
  ExportNgsx_xfem(m);
}
