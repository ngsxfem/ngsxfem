#include <python_ngstd.hpp>
#include "../xfem/sFESpace.hpp"
#include "../xfem/cutinfo.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/ghostpenalty.hpp"

using namespace ngcomp;

void ExportNgsx_xfem(py::module &m)
{


  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<SFESpace> PySFES;

  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  typedef shared_ptr<CutInformation> PyCI;
  typedef shared_ptr<BitArray> PyBA;
  
  
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
                             py::object definedonelem)
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
        docu_string(R"raw_string(
see documentation of SymbolicBFI (which is a wrapper))raw_string")
    );

  m.def("SymbolicFacetPatchBFI", [](PyCF cf,
                                    int order,
                                    int time_order,
                                    bool skeleton,
                                    py::object definedonelem)
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

          return PyBFI(bfi);
        },
        py::arg("form"),
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("skeleton") = true,
        py::arg("definedonelements")=DummyArgument(),
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
                             py::object definedonelem)
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
  
}

PYBIND11_MODULE(ngsxfem_xfem_py,m)
{
  cout << "importing ngsxfem-xfem lib" << endl;
  ExportNgsx_xfem(m);
}
