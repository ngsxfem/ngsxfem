#include <regex>
#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

//using namespace ngcomp;
#include "../utils/ngsxstd.hpp"
#include "../cutint/straightcutrule.hpp"
#include "../cutint/xintegration.hpp"
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
          static Timer t ("IntegrateX"); RegionTimer reg(t);
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

               const IntegrationRule * ir;
               Array<double> wei_arr;
               tie (ir, wei_arr) = CreateCutIntegrationRule(cf_lset, gf_lset, trafo, dt, order, time_order, lh, subdivlvl, quad_dir_policy);

               if (ir != nullptr)
               {
                 BaseMappedIntegrationRule & mir = trafo(*ir, lh);
                 FlatMatrix<> val(mir.Size(), 1, lh);

                 cf -> Evaluate (mir, val);

                 double lsum = 0.0;
                 for (int i = 0; i < mir.Size(); i++)
                     lsum += mir[i].GetMeasure()*wei_arr[i]*val(i,0);
                   //lsum += mir[i].GetWeight()*val(i,0);

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

  m.def("InterpolateToP1",  [] (PyGF gf_ho, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(gf_ho, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("gf_ho")=NULL,py::arg("gf_p1")=NULL,
        py::arg("eps_perturbation")=1e-14,py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Takes the vertex values of a GridFunction (also possible with a CoefficentFunction) and puts them
into a piecewise (multi-) linear function.

Parameters

gf_ho : ngsolve.GridFunction
  Function to interpolate

gf_p1 : ngsolve.GridFunction
  Function to interpolate to (should be P1)

eps_perturbation : float
  If the absolute value if the function is smaller than eps_perturbation, it will be set to
  eps_perturbation. Thereby, exact and close-to zeros at vertices are avoided (Useful to reduce cut
  configurations for level set based methods).

heapsize : int
  heapsize of local computations.
)raw_string")
    )
    ;

  m.def("InterpolateToP1",  [] (PyCF coef, PyGF gf_p1, double eps_perturbation, int heapsize)
        {
          InterpolateP1 interpol(coef, gf_p1);
          LocalHeap lh (heapsize, "InterpolateP1-Heap");
          interpol.Do(lh,eps_perturbation);
        } ,
        py::arg("coef"),py::arg("gf"),
        py::arg("eps_perturbation")=1e-14,py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Takes the vertex values of a CoefficentFunction) and puts them into a piecewise (multi-) linear
function.

Parameters

coef : ngsolve.CoefficientFunction
  Function to interpolate

gf_p1 : ngsolve.GridFunction
  Function to interpolate to (should be P1)

eps_perturbation : float
  If the absolute value if the function is smaller than eps_perturbation, it will be set to
  eps_perturbation. Thereby, exact and close-to zeros at vertices are avoided (Useful to reduce cut
  configurations for level set based methods).

heapsize : int
  heapsize of local computations.
)raw_string")
    )
    ;

  
  typedef shared_ptr<BitArray> PyBA;

  m.def("RestrictedBilinearForm",
        [](shared_ptr<FESpace> fes,
           const string & aname,
           py::object ael_restriction,
           py::object afac_restriction,
           bool check_unused,
           py::dict bpflags)
        {
          Flags flags = py::extract<Flags> (bpflags)();

          shared_ptr<BitArray> el_restriction = nullptr;
          shared_ptr<BitArray> fac_restriction = nullptr;
          if (py::extract<PyBA> (ael_restriction).check())
            el_restriction = py::extract<PyBA>(ael_restriction)();

          if (py::extract<PyBA> (afac_restriction).check())
            fac_restriction = py::extract<PyBA>(afac_restriction)();

          if (fes->IsComplex())
            throw Exception("RestrictedBilinearForm not implemented for complex fespace");

          shared_ptr<BilinearForm> biform = make_shared<RestrictedBilinearForm> (fes, aname, el_restriction, fac_restriction, flags);
          biform -> SetCheckUnused (check_unused);
          return biform;
        },
        py::arg("space"),
        py::arg("name") = "bfa",
        py::arg("element_restriction") = DummyArgument(),
        py::arg("facet_restriction") = DummyArgument(),
        py::arg("check_unused") = true,
        py::arg("flags") = py::dict(),
        docu_string(R"raw_string(
A restricted bilinear form is a (so far real-valued) bilinear form with a reduced MatrixGraph
compared to the usual BilinearForm. BitArray(s) define on which elements/facets entries will be
created.

Use cases:

 * ghost penalty type stabilization:
    Facet-stabilization that are introduced only act on a few facets in the mesh. By providing the
    information on the corresponding facets, these additional couplings will only be introduced
    where necessary.

 * fictitious domain methods:
    When PDE problems are only solved on a part of a domain while a finite element space is used
    that is still defined on the whole domain, a BitArray can be used to mark the 'active' part of
    the mesh.

Parameters

space : ngsolve.FESpace
  finite element space on which the bilinear form is defined.

name : string
  name of the bilinear form

element_restriction : ngsolve.BitArray
  BitArray defining the 'active mesh' element-wise

facet_restriction : ngsolve.BitArray
  BitArray defining the 'active facets'. This is only relevant if FESpace has DG-terms (dgjumps=True)

check_unused : boolean
  Check if some degrees of freedoms are not considered during assembly

flags : ngsolve.Flags
  additional bilinear form flags
)raw_string"));

  m.def("CompoundBitArray",
        [] (py::list balist)
        {
          size_t cnt = 0;
          for( auto aba : balist )
          {
            shared_ptr<BitArray> ba = py::extract<PyBA>(aba)();
            cnt += ba->Size();
          }
          shared_ptr<BitArray> res = make_shared<BitArray>(cnt);
          res->Clear();
          size_t offset = 0;
          for( auto aba : balist )
          {
            shared_ptr<BitArray> ba = py::extract<PyBA>(aba)();
            for (size_t i = 0; i < ba->Size(); ++i)
            {
              if (ba->Test(i))
                res->Set(offset+i);
            }
            offset += ba->Size();
          }
          return res;
        } ,
        py::arg("balist"),
        docu_string(R"raw_string(
Takes a list of BitArrays and merges them to one larger BitArray. Can be useful for
CompoundFESpaces.
)raw_string")
    );



  typedef shared_ptr<BitArrayCoefficientFunction> PyBACF;
  py::class_<BitArrayCoefficientFunction, PyBACF, CoefficientFunction>
    (m, "BitArrayCF",
        docu_string(R"raw_string(
CoefficientFunction that evaluates a BitArray. On elements with an index i where the BitArray
evaluates to true the CoefficientFunction will evaluate as 1, otherwise as 0.

Similar functionality (also for facets) can be obtained with IndicatorCF.
)raw_string"))
    .def("__init__",
         [](BitArrayCoefficientFunction *instance, shared_ptr<BitArray> ba)
         {
           new (instance) BitArrayCoefficientFunction (ba);
         },
         py::arg("bitarray")
      );
   
  
  typedef shared_ptr<SFESpace> PySFES;
  typedef shared_ptr<CutInformation> PyCI;
  
  m.def("SFESpace", [](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
        -> PyFES
        {
          Flags flags = py::extract<Flags> (bpflags)();
          shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset, order, flags);
          LocalHeap lh (1000000, "SFESpace::Update-heap", true);
          ret->Update(lh);
          ret->FinalizeUpdate(lh);
          return ret;
        },
        docu_string(R"raw_string(
This is a special finite elemetn space which is a 1D polynomial along the zero level of the linearly
approximated level set function lset and constantly extended in normal directions to this.
)raw_string"));

  typedef shared_ptr<XFESpace> PyXFES;

  m.def("XToNegPos",  [] (PyGF gfx, PyGF gfnegpos) {
      XFESpace::XToNegPos(gfx,gfnegpos);
    }, docu_string(R"raw_string(
Takes a GridFunction of an extended FESpace, i.e. a compound space of V and VX = XFESpace(V) and
interpretes it as a function in the CompoundFESpace of V and V. Updates the values of the vector of
the corresponding second GridFunction.
)raw_string") );

  py::class_<CutInformation, shared_ptr<CutInformation>>
    (m, "CutInfo",R"raw(
A CutInfo stores and organizes cut informations in the mesh with respect to a level set function. 
Elements (BND / VOL) and facets can be either cut elements or in the positive (POS) or negative
(NEG) part of the domain. A CutInfo provides information about the cut configuration in terms of
BitArrays and Vectors of Ratios. (Internally also domain_types for different mesh nodes are stored.)
)raw")
    .def("__init__",  [] (CutInformation *instance,
                          shared_ptr<MeshAccess> ma,
                          py::object lset,
                          int time_order,
                          int heapsize)
         {
           new (instance) CutInformation (ma);
           if (py::extract<PyCF> (lset).check())
           {
             PyCF cflset = py::extract<PyCF>(lset)();
             LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
             instance->Update(cflset, time_order, lh);
           }
         },
         py::arg("mesh"),
         py::arg("levelset") = DummyArgument(),
         py::arg("time_order") = -1,
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Creates a CutInfo based on a level set function and a mesh.

Parameters

mesh : Mesh

levelset : ngsolve.CoefficientFunction / None
  level set funciton w.r.t. which the CutInfo is created

time_order : int
  order in time that is used in the integration in time to check for cuts and the ratios. This is
  only relevant for space-time discretizations.
)raw_string")
      )
    .def("Update", [](CutInformation & self,
                      PyCF lset,
                      int time_order,
                      int heapsize)
         {
           LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
           self.Update(lset,time_order,lh);
         },
         py::arg("levelset"),
         py::arg("time_order") = -1,
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Updates a CutInfo based on a level set function.

Parameters

levelset : ngsolve.CoefficientFunction
  level set function w.r.t. which the CutInfo is generated

time_order : int
  order in time that is used in the integration in time to check for cuts and the ratios. This is
  only relevant for space-time discretizations.


)raw_string")
      )
    .def("Mesh", [](CutInformation & self)
         {
           return self.GetMesh();
         },docu_string(R"raw_string(
Returns mesh of CutInfo)raw_string")
      )
    .def("GetElementsOfType", [](CutInformation & self,
                                 py::object dt,
                                 VorB vb)
         {
           COMBINED_DOMAIN_TYPE cdt = CDOM_NO;
           if (py::extract<COMBINED_DOMAIN_TYPE> (dt).check())
             cdt = py::extract<COMBINED_DOMAIN_TYPE>(dt)();
           else if (py::extract<DOMAIN_TYPE> (dt).check())
             cdt = TO_CDT(py::extract<DOMAIN_TYPE>(dt)());
           else
             throw Exception(" unknown type for dt ");
           return self.GetElementsOfDomainType(cdt,vb);
         },
         py::arg("domain_type") = IF,
         py::arg("VOL_or_BND") = VOL,docu_string(R"raw_string(
Returns BitArray that is true for every element that has the 
corresponding combined domain type 
(NO/NEG/POS/UNCUT/IF/HASNEG/HASPOS/ANY))raw_string")
      )
    .def("GetFacetsOfType", [](CutInformation & self,
                               py::object dt)
         {
           COMBINED_DOMAIN_TYPE cdt = CDOM_NO;
           if (py::extract<COMBINED_DOMAIN_TYPE> (dt).check())
             cdt = py::extract<COMBINED_DOMAIN_TYPE>(dt)();
           else if (py::extract<DOMAIN_TYPE> (dt).check())
             cdt = TO_CDT(py::extract<DOMAIN_TYPE>(dt)());
           else
             throw Exception(" unknown type for dt ");
           return self.GetFacetsOfDomainType(cdt);
         },
         py::arg("domain_type") = IF,docu_string(R"raw_string(
Returns BitArray that is true for every facet that has the 
corresponding combined domain type 
(NO/NEG/POS/UNCUT/IF/HASNEG/HASPOS/ANY))raw_string")
      )
    .def("GetCutRatios", [](CutInformation & self,
                            VorB vb)
         {
           return self.GetCutRatios(vb);
         },
         py::arg("VOL_or_BND") = VOL,docu_string(R"raw_string(
Returns Vector of the ratios between the measure of the NEG domain on a (boundary) element and the
full (boundary) element
)raw_string"))
    ;


  m.def("GetFacetsWithNeighborTypes",
        [] (shared_ptr<MeshAccess> ma,
            shared_ptr<BitArray> a,
            bool bv_a,
            bool bv_b,
            bool use_and,
            py::object bb,
            int heapsize)
        {
          LocalHeap lh (heapsize, "FacetsWithNeighborTypes-heap", true);
          shared_ptr<BitArray> b = nullptr;
          if (py::extract<PyBA> (bb).check())
            b = py::extract<PyBA>(bb)();
          else
            b = a;
          return GetFacetsWithNeighborTypes(ma,a,b,bv_a,bv_b,use_and,lh);
        } ,
        py::arg("mesh"),
        py::arg("a"),
        py::arg("bnd_val_a") = true,
        py::arg("bnd_val_b") = true,
        py::arg("use_and") = true,
        py::arg("b") = DummyArgument(),
        py::arg("heapsize") = 1000000, docu_string(R"raw_string(
Given a mesh and two BitArrays (if only one is provided these are set to be equal) facets will be
marked (in terms of BitArrays) depending on the BitArray-values on the neighboring elements. The
BitArrays are complemented with flags for potential boundary values for the BitArrays. The decision
on every facet is now based on the values a and b (left and right) where a or b can also be obtained
from the BitArray boundary values.
The result is:
  result =    (a(left) and b(right)) 
           or (b(left) and a(right)) 
or 
  result =    (a(left) or b(right)) 
           or (b(left) or a(right)) 

Parameters:

mesh : 
  mesh

a : ngsolve.BitArray
  first BitArray 

b : ngsolve.BitArray / None
  second BitArray. If None, b=a

bnd_val_a : boolean
  BitArray-replacement for a if a(left) or a(right) is not valid (at the boundary)

bnd_val_a : boolean
  BitArray-replacement for b if b(left) or b(right) is not valid (at the boundary)

use_and : boolean
  use 'and'-relation to evaluate the result. Otherwise use 'or'-relation 

heapsize : int
  heapsize of local computations.
)raw_string")
    );

  m.def("GetElementsWithNeighborFacets",
        [] (shared_ptr<MeshAccess> ma,
            shared_ptr<BitArray> a,
            int heapsize)
        {
          LocalHeap lh (heapsize, "GetElementsWithNeighborFacets-heap", true);
          return GetElementsWithNeighborFacets(ma,a,lh);
        } ,
        py::arg("mesh"),
        py::arg("a"),
        py::arg("heapsize") = 1000000,
        docu_string(R"raw_string(
Given a BitArray marking some facets extract
a BitArray of elements that are neighboring
these facets

Parameters:

mesh : 
  mesh

a : ngsolve.BitArray
  BitArray for marked facets

heapsize : int
  heapsize of local computations.
)raw_string")
    );

  m.def("GetDofsOfElements",
        [] (PyFES fes,
            PyBA a,
            int heapsize)
        {
          LocalHeap lh (heapsize, "GetDofsOfElements-heap", true);
          return GetDofsOfElements(fes,a,lh);
        } ,
        py::arg("space"),
        py::arg("a"),
        py::arg("heapsize") = 1000000,
        docu_string(R"raw_string(
Given a BitArray marking some elements in a
mesh extract all unknowns that are supported
on these elements as a BitArray.

Parameters:

space : ngsolve.FESpace
  finite element space from which the 
  corresponding dofs should be extracted

a : ngsolve.BitArray
  BitArray for marked elements

heapsize : int
  heapsize of local computations.
)raw_string")

    );

  m.def("GetDofsOfFacets",
        [] (PyFES fes,
            PyBA a,
            int heapsize)
        {
          LocalHeap lh (heapsize, "GetDofsOfFacets-heap", true);
          return GetDofsOfFacets(fes,a,lh);
        } ,
        py::arg("space"),
        py::arg("a"),
        py::arg("heapsize") = 1000000,
        docu_string(R"raw_string(
Given a BitArray marking some facets in a
mesh extract all unknowns that are associated
to these facets as a BitArray.

Parameters:

space : ngsolve.FESpace
  finite element space from which the 
  corresponding dofs should be extracted

a : ngsolve.BitArray
  BitArray for marked Facets

heapsize : int
  heapsize of local computations.
)raw_string")

    );


//   .def("__init__",  [] (XFESpace *instance,
  m.def("XFESpace", [] (
          PyFES basefes,
          py::object acutinfo,
          py::object alset,
          py::dict bpflags,
          int heapsize)
        {
          shared_ptr<CoefficientFunction> cf_lset = nullptr;
          shared_ptr<CutInformation> cutinfo = nullptr;
          if (py::extract<PyCI> (acutinfo).check())
            cutinfo = py::extract<PyCI>(acutinfo)();
          if (py::extract<PyCF> (acutinfo).check())
            cf_lset = py::extract<PyCF>(acutinfo)();
          if (py::extract<PyCF> (alset).check())
            cf_lset = py::extract<PyCF>(alset)();


          Flags flags = py::extract<Flags> (bpflags)();
          shared_ptr<XFESpace> ret = nullptr;
          shared_ptr<MeshAccess> ma = basefes->GetMeshAccess();
          if (cutinfo)
          {
            if (ma->GetDimension()==2)
              ret = make_shared<T_XFESpace<2>> (ma, basefes, cutinfo, flags);
            else
              ret = make_shared<T_XFESpace<3>> (ma, basefes, cutinfo, flags);
          }
          else if (cf_lset)
          {
            if (ma->GetDimension()==2)
              ret = make_shared<T_XFESpace<2>> (ma, basefes, cf_lset, flags);
            else
              ret = make_shared<T_XFESpace<3>> (ma, basefes, cf_lset, flags);
          }
          else
            throw Exception("levelset and cutinfo are invalid");
          LocalHeap lh (heapsize, "XFESpace::Update-heap", true);
          ret->Update(lh);
          return ret;
        },
        py::arg("basefes"),
        py::arg("cutinfo") = DummyArgument(),
        py::arg("lset") = DummyArgument(),
        py::arg("flags") = py::dict(),
        py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Constructor for XFESpace [For documentation of XFESpace-class see help(CXFESpace)]:

Extended finite element space. Takes a basis FESpace and creates an enrichment space based on cut
information. The cut information is provided by a CutInfo object or - if a level set function is
only provided - a CutInfo object is created. The enrichment doubles the unknowns on all cut elements
and assigns to them a sign (NEG/POS). One of the differential operators neg(...) or pos(...)
evaluates like the basis function of the origin space, the other one as zero for every basis
function. Away from cut elements no basis function is supported.

Parameters

basefes : ngsolve.FESpace
  basic FESpace to be extended

cutinfo : xfem.CutInfo / None
  Information on the cut configurations (cut elements, sign of vertices....)

lset : ngsolve.CoefficientFunction / None
  level set function to construct own CutInfo (if no CutInfo is provided)

flags : Flags
  additional FESpace-flags

heapsize : int
  heapsize of local computations.
)raw_string"));
         

  py::class_<XFESpace, PyXFES, FESpace>
    (m, "CXFESpace",docu_string(R"raw_string(
XFESpace-class [For documentation of the XFESpace-constructor see help(XFESpace)]:

Extended finite element space. Takes a basis FESpace and creates an enrichment space based on cut
information.  The cut information is provided by a CutInfo object or - if a level set function is
only provided - a CutInfo object is created. The enrichment doubles the unknowns on all cut elements
and assigns to them a sign (NEG/POS). One of the differential operators neg(...) or pos(...)
evaluates like the basis function of the origin space, the other one as zero for every basis
function. Away from cut elements no basis function is supported.
)raw_string"))
    .def("GetCutInfo", [](PyXFES self)
         {
           return self->GetCutInfo();
         },
         "Get Information of cut geometry")
    .def("BaseDofOfXDof", [](PyXFES self, int i)
         {
           return self->GetBaseDofOfXDof(i);
         },docu_string(R"raw_string(
To an unknown of the extended space, get the corresponding unknown of the base FESpace.

Parameters

i : int
  degree of freedom 
)raw_string"))
    .def("GetDomainOfDof", [](PyXFES self, int i)
         {
           return self->GetDomainOfDof(i);
         },docu_string(R"raw_string(
Get Domain (NEG/POS) of a degree of freedom of the extended FESpace.

Parameters

i : int
  degree of freedom 
)raw_string"))
    .def("GetDomainNrs",   [] (PyXFES self, int elnr) {
        Array<DOMAIN_TYPE> domnums;
        self->GetDomainNrs( elnr, domnums );
        return domnums;
      },docu_string(R"raw_string(
Get Array of Domains (Array of NEG/POS) of degrees of freedom of the extended FESpace on one element.

Parameters

elnr : int
  element number
)raw_string"))
    ;

  typedef shared_ptr<BilinearFormIntegrator> PyBFI;
  typedef shared_ptr<LinearFormIntegrator> PyLFI;

  m.def("SymbolicCutBFI", [](PyCF lset,
                             DOMAIN_TYPE dt,
                             int order,
                             int time_order,
                             int subdivlvl,
                             SWAP_DIMENSIONS_POLICY quad_dir_pol,
                             PyCF cf,
                             VorB vb,
                             bool element_boundary,
                             bool skeleton,
                             py::object definedon,
                             py::object definedonelem,
                             py::object deformation)
        -> PyBFI
        {

          py::extract<Region> defon_region(definedon);
          if (defon_region.check())
            vb = VorB(defon_region());

          // check for DG terms
          bool has_other = false;
          cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                            {
                              if (dynamic_cast<ProxyFunction*> (&cf))
                                if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                  has_other = true;
                            });
          if (has_other && !element_boundary && !skeleton)
            throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");
          if (element_boundary)
            throw Exception("No Facet BFI with Symbolic cuts..");

          shared_ptr<BilinearFormIntegrator> bfi;
          if (!has_other && !skeleton)
          {
            auto bfime = make_shared<SymbolicCutBilinearFormIntegrator> (lset, cf, dt, order, subdivlvl,quad_dir_pol,vb);
            bfime->SetTimeIntegrationOrder(time_order);
            bfi = bfime;
          }
          else
          {
            if (time_order >= 0)
              throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for time_order >= 0..");
            if (vb == BND)
              throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for boundaries..");
            bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (lset, cf, dt, order, subdivlvl);
          }
          if (py::extract<py::list> (definedon).check())
            bfi -> SetDefinedOn (makeCArray<int> (definedon));

          if (defon_region.check())
          {
            cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
            bfi->SetDefinedOn(defon_region().Mask());
          }

          if (! py::extract<DummyArgument> (definedonelem).check())
            bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! py::extract<DummyArgument> (deformation).check())
            bfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyBFI(bfi);
        },
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("subdivlvl")=0,
        py::arg("quad_dir_policy")=FIND_OPTIMAL,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=DummyArgument(),
        py::arg("definedonelements")=DummyArgument(),
        py::arg("deformation")=DummyArgument(),
        docu_string(R"raw_string(
see documentation of SymbolicBFI (which is a wrapper))raw_string")
    );

  m.def("SymbolicFacetPatchBFI", [](PyCF cf,
                                    int order,
                                    int time_order,
                                    bool skeleton,
                                    py::object definedonelem,
                                    py::object deformation)
        -> PyBFI
        {
          // check for DG terms
          bool has_other = false;
          cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                            {
                              if (dynamic_cast<ProxyFunction*> (&cf))
                                if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                  has_other = true;
                            });
          if (!has_other)
            cout << " no Other() used?!" << endl;

  
          shared_ptr<BilinearFormIntegrator> bfi;
          if (skeleton)
          {
            auto bfime = make_shared<SymbolicFacetBilinearFormIntegrator2> (cf, order);
            bfime->SetTimeIntegrationOrder(time_order);
            bfi = bfime;
          }
          else
          {
            // throw Exception("Patch facet blf not implemented yet: TODO(2)!");
            auto bfime = make_shared<SymbolicFacetPatchBilinearFormIntegrator> (cf, order);
            bfime->SetTimeIntegrationOrder(time_order);
            bfi = bfime;
          }

          if (! py::extract<DummyArgument> (definedonelem).check())
            bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! py::extract<DummyArgument> (deformation).check())
            bfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyBFI(bfi);
        },
        py::arg("form"),
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("skeleton") = true,
        py::arg("definedonelements")=DummyArgument(),
        py::arg("deformation")=DummyArgument(),
        docu_string(R"raw_string(
Integrator on facet patches. Two versions are possible:
* Either (skeleton=False) an integration on the element patch consisting of two neighboring elements is applied, 
* or (skeleton=True) the integration is applied on the facet. 

Parameters

form : ngsolve.CoefficientFunction
  var form to integrate

force_intorder : int
  (only active in the facet patch case (skeleton=False)) use this integration order in the integration

skeleton : boolean
  decider on facet patch vs facet integration

definedonelements : ngsolve.BitArray/None
  array which decides on which facets the integrator should be applied

time_order : int
  order in time that is used in the space-time integration. time_order=-1 means that no space-time
  rule will be applied. This is only relevant for space-time discretizations.
)raw_string")
    );

  m.def("SymbolicCutLFI", [](PyCF lset,
                             DOMAIN_TYPE dt,
                             int order,
                             int time_order,
                             int subdivlvl,
                             SWAP_DIMENSIONS_POLICY quad_dir_pol,
                             PyCF cf,
                             VorB vb,
                             bool element_boundary,
                             bool skeleton,
                             py::object definedon,
                             py::object definedonelem,
                             py::object deformation)
        -> PyLFI
        {

          py::extract<Region> defon_region(definedon);
          if (defon_region.check())
            vb = VorB(defon_region());

          // if (vb == BND)
          //   throw Exception("Symbolic cuts not yet (tested) for boundaries..");

          if (element_boundary || skeleton)
            throw Exception("No Facet LFI with Symbolic cuts..");

          auto lfime  = make_shared<SymbolicCutLinearFormIntegrator> (lset, cf, dt, order, subdivlvl, quad_dir_pol,vb);
          lfime->SetTimeIntegrationOrder(time_order);
          shared_ptr<LinearFormIntegrator> lfi = lfime;

          if (py::extract<py::list> (definedon).check())
            lfi -> SetDefinedOn (makeCArray<int> (definedon));

          if (defon_region.check())
          {
            cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
            lfi->SetDefinedOn(defon_region().Mask());
          }

          if (! py::extract<DummyArgument> (definedonelem).check())
            lfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! py::extract<DummyArgument> (deformation).check())
            lfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyLFI(lfi);
        },
        py::arg("lset"),
        py::arg("domain_type")=NEG,
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("subdivlvl")=0,
        py::arg("quad_dir_policy")=FIND_OPTIMAL,
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=py::bool_(false),
        py::arg("skeleton")=py::bool_(false),
        py::arg("definedon")=DummyArgument(),
        py::arg("definedonelements")=DummyArgument(),
        py::arg("deformation")=DummyArgument(),
        docu_string(R"raw_string(
see documentation of SymbolicLFI (which is a wrapper))raw_string")
    );

  typedef shared_ptr<ProxyFunction> PyProxyFunction;
  m.def("dn", [] (const PyProxyFunction self, int order, py::object comp, int dim_space, bool hdiv)
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

          shared_ptr<DifferentialOperator> diffopdudnk;
          if (! hdiv)
          {
            if (dim_space == 2)
            {
              switch (order)
              {
              case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,1>>> (); break;
              case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,2>>> (); break;
              case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,3>>> (); break;
              case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
              case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,5>>> (); break;
              case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,6>>> (); break;
              case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,7>>> (); break;
              case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,8>>> (); break;
              default : throw Exception("no order higher than 8 implemented yet");
              }
            }
            else
            {
              switch (order)
              {
              case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,1>>> (); break;
              case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,2>>> (); break;
              case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,3>>> (); break;
              case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,4>>> (); break;
              case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,5>>> (); break;
              case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,6>>> (); break;
              case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,7>>> (); break;
              case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<3,8>>> (); break;
              default : throw Exception("no order higher than 8 implemented yet");
              }
            }
          }
          else
            switch (order)
            {
            case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,1>>> (); break;
            case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,2>>> (); break;
            case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,3>>> (); break;
            case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,4>>> (); break;
            case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,5>>> (); break;
            case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,6>>> (); break;
            case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,7>>> (); break;
            case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnkHDiv<2,8>>> (); break;
            default : throw Exception("no order higher than 8 implemented yet");
            }

          for (int i = comparr.Size() - 1; i >= 0; --i)
          {
            diffopdudnk = make_shared<CompoundDifferentialOperator> (diffopdudnk, comparr[i]);
          }

          auto adddiffop = make_shared<ProxyFunction> (self->GetFESpace(),self->IsTestFunction(), self->IsComplex(),
                                                       diffopdudnk, nullptr, nullptr, nullptr, nullptr, nullptr);

          if (self->IsOther())
            adddiffop = adddiffop->Other(make_shared<ConstantCoefficientFunction>(0.0));

          return PyProxyFunction(adddiffop);
        },
        py::arg("proxy"),
        py::arg("order"),
        py::arg("comp") = -1,
        py::arg("dim_space") = 2,
        py::arg("hdiv") = false,
        docu_string(R"raw_string(
Normal derivative of higher order. This is evaluated via numerical differentiation which offers only
limited accuracy (~ 1e-7).

Parameters

proxy : ngsolve.ProxyFunction
  test / trialfunction to the the normal derivative of

order : int
  order of derivative (in normal direction)

comp : int
  component of proxy if test / trialfunction is a component of a compound or vector test / trialfunction

dim_space : int
  dimension of the space

hdiv : boolean
  assumes scalar FEs if false, otherwise assumes hdiv
)raw_string")
    );

  m.def("dn", [](PyGF gf, int order) -> PyCF
        {
          shared_ptr<DifferentialOperator> diffopdudnk;
          switch (order)
          {
          case 1 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,1>>> (); break;
          case 2 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,2>>> (); break;
          case 3 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,3>>> (); break;
          case 4 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,4>>> (); break;
          case 5 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,5>>> (); break;
          case 6 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,6>>> (); break;
          case 7 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,7>>> (); break;
          case 8 : diffopdudnk = make_shared<T_DifferentialOperator<DiffOpDuDnk<2,8>>> (); break;
          default : throw Exception("no order higher than 8 implemented yet");
          }
          return PyCF(make_shared<GridFunctionCoefficientFunction> (gf, diffopdudnk));
        },
        py::arg("gf"),
        py::arg("order"),
        docu_string(R"raw_string(
Normal derivative of higher order for a GridFunction. This is evaluated via numerical
differentiation which offers only limited accuracy (~ 1e-7).

Parameters

gf : ngsolve.GridFunction
  (scalar) GridFunction to the the normal derivative of

order : int
  order of derivative (in normal direction)
)raw_string")
);


  typedef shared_ptr<P1Prolongation> PyP1P;
  py::class_<P1Prolongation, PyP1P, Prolongation>
    (m, "P1Prolongation",
        docu_string(R"raw_string(
Prolongation for P1-type spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P1Prolongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P1Prolongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P1Prolongation> p1p, shared_ptr<FESpace> fes)
         {
           p1p -> Update(*fes);
         },
         py::arg("space")
      );

      typedef shared_ptr<P2Prolongation> PyP2P;
  py::class_<P2Prolongation, PyP2P, Prolongation>
    (m, "P2Prolongation",
        docu_string(R"raw_string(
Prolongation for P2 spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P2Prolongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P2Prolongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P2Prolongation> p2p, shared_ptr<FESpace> fes)
         {
           p2p -> Update(*fes);
         },
         py::arg("space")
      );


      typedef shared_ptr<P2CutProlongation> PyP2CutP;
  py::class_<P2CutProlongation, PyP2CutP, Prolongation>
    (m, "P2CutProlongation",
        docu_string(R"raw_string(
Prolongation for P2 spaces (with possibly inactive dofs) --- 
As is asks the fespace for dofs to vertices at several occasions the 
current implementation is not very fast and should be primarily used
for prototype and testing...
)raw_string"))
    .def("__init__",
         [](P2CutProlongation *instance, shared_ptr<MeshAccess> ma)
         {
           new (instance) P2CutProlongation (ma);
         },
         py::arg("mesh"))
    .def("Update",
         [](shared_ptr<P2CutProlongation> p2p, shared_ptr<FESpace> fes)
         {
           p2p -> Update(*fes);
         },
         py::arg("space")
      );

    typedef shared_ptr<CompoundProlongation> PyCProl;
    py::class_< CompoundProlongation, PyCProl, Prolongation>
    ( m, "CompoundProlongation", 
     docu_string(R"raw_string(prolongation for compound spaces)raw_string"))
     .def("__init__",
          [](CompoundProlongation *instance, const FESpace *fes)
          {
            new (instance) CompoundProlongation( dynamic_cast<const CompoundFESpace*>(fes) );
          },
          py::arg("compoundFESpace"))
      .def("Update", &CompoundProlongation::Update, py::arg("fespace"))
      .def ("Prolongate", &CompoundProlongation::ProlongateInline, py::arg("finelevel"), py::arg("vec"))
      .def ("Restrict", &CompoundProlongation::RestrictInline, py::arg("finelevel"), py::arg("vec"))
      .def ("AddProlongation", [](shared_ptr<CompoundProlongation> cprol, shared_ptr<Prolongation> prol )
            { cprol -> AddProlongation( prol ); }, py::arg("p1prol")
          );
      //.def ("AddProlongation" &CompoundProlongation::AddProlongation, py::arg("prol"));
  

}

PYBIND11_MODULE(ngsxfem_py, m)
{
  cout << "importing ngs-xfem" << NGSXFEM_VERSION << endl;
  ExportNgsx(m);
}
