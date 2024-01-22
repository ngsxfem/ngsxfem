#include <python_ngstd.hpp>
#include "../xfem/sFESpace.hpp"
#include "../xfem/cutinfo.hpp"
#include "../xfem/aggregates.hpp"
#include "../xfem/xFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../xfem/ghostpenalty.hpp"

#include <typeinfo>

using namespace ngcomp;

void ExportNgsx_xfem(py::module &m)
{


  typedef shared_ptr<CoefficientFunction> PyCF;
  typedef shared_ptr<FESpace> PyFES;
  typedef shared_ptr<SFESpace> PySFES;

  typedef GridFunction GF;
  typedef shared_ptr<GF> PyGF;
  typedef shared_ptr<CutInformation> PyCI;
  typedef shared_ptr<MultiLevelsetCutInformation> PyMLCI;
  typedef shared_ptr<BitArray> PyBA;
  
  
  m.def("SFESpace", [](shared_ptr<MeshAccess> ma, PyCF lset, int order, py::dict bpflags)
        -> PyFES
        {
          Flags flags = py::extract<Flags> (bpflags)();
          shared_ptr<FESpace> ret = make_shared<SFESpace> (ma, lset, order, flags);
          LocalHeap lh (1000000, "SFESpace::Update-heap", true);
          ret->Update();
          ret->FinalizeUpdate();
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
                          int subdivlvl,
                          int time_order,
                          int heapsize)
         {
           new (instance) CutInformation (ma);
           if ((!lset.is_none()) && py::extract<PyCF> (lset).check())
           {
             PyCF cflset = py::extract<PyCF>(lset)();
             LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
             instance->Update(cflset, subdivlvl, time_order, lh);
           }
         },
         py::arg("mesh"),
         py::arg("levelset") = py::none(),
         py::arg("subdivlvl") = 0,
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
                      int subdivlvl,
                      int time_order,
                      int heapsize)
         {
           LocalHeap lh (heapsize, "CutInfo::Update-heap", true);
           self.Update(lset,subdivlvl,time_order,lh);
         },
         py::arg("levelset"),
         py::arg("subdivlvl") = 0,
         py::arg("time_order") = -1,
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Updates a CutInfo based on a level set function.

Parameters

levelset : ngsolve.CoefficientFunction
  level set function w.r.t. which the CutInfo is generated

subdivlvl : int
  subdivision for numerical integration

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
    .def("GetElementsWithThresholdContribution", [](CutInformation & self,
                                                    py::object dt,
                                                    double threshold,
                                                    VorB vb)
         {
           DOMAIN_TYPE _dt = NEG;
           if (py::extract<DOMAIN_TYPE> (dt).check() && py::extract<DOMAIN_TYPE> (dt)() != IF)
              _dt = py::extract<DOMAIN_TYPE>(dt)();
           else
              throw Exception("Unknown/Invalid type for dt: Only POS, NEG are implemented a.t.m.");
           return self.GetElementsWithThresholdContribution(_dt, threshold, vb);
         },
         py::arg("domain_type") = NEG,
         py::arg("threshold") = 1.0,
         py::arg("VOL_or_BND") = VOL, docu_string(R"raw_string(
Returns BitArray marking the elements where the cut ratio is greater or equal to the given 
threshold.

Parameters

domain_type : ENUM
    Check POS or NEG elements.

threshold : float
    Mark elements with cut ratio (volume of domain_type / volume background mesh) greater or equal to threshold.

VOL_or_BND : ngsolve.comp.VorB
    input VOL, BND, ..

)raw_string"))
    ;


py::class_<ElementAggregation, shared_ptr<ElementAggregation>>
    (m, "ElementAggregation",R"raw(
ElementAggregation does the following:  
  It collects patches of elements that allow to stabilize bad cut elements by at least one
  good element (the root element).
)
)raw")
    .def("__init__",  [] (ElementAggregation *instance,
                          shared_ptr<MeshAccess> ma,
                          py::object proot, 
                          py::object pbad,
                          int heapsize                          
                          )
         {
           auto self = new (instance) ElementAggregation (ma);
           PyBA root = nullptr, bad= nullptr;
           if (!proot.is_none() && py::extract<PyBA> (proot).check())
             root = py::extract<PyBA>(proot)();
           if (!pbad.is_none() && py::extract<PyBA> (pbad).check())
             bad = py::extract<PyBA>(pbad)();

           if (root && bad)
           {
             LocalHeap lh (heapsize, "ElementAggregation::Update-heap", true);
             self->Update(root, bad, lh);
           }
         },
         py::arg("mesh"),
         py::arg("root_elements") = py::none(),
         py::arg("bad_elements") = py::none(),
         py::arg("heapsize") = 1000000,
         docu_string(R"raw_string(
Creates a ElementAggregation based on a mesh, a list of root and a list of bad elements.
)raw_string")
      )
        .def("Update", [](ElementAggregation & self,
                          PyBA root, 
                          PyBA bad,
                          int heapsize)
         {
           LocalHeap lh (heapsize, "ElementAggregation::Update-heap", true);
           self.Update(root,bad,lh);
         },
         py::arg("root_elements"),
         py::arg("bad_elements"),
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Updates a Element Aggregation based ...
)raw_string")
      )
      .def_property_readonly("patch_interior_facets",
                [](ElementAggregation & self) { return self.GetInnerPatchFacets(); },
               "BitArray that is true for every facet that is *inside* an aggregation cluster")      
      .def_property_readonly("els_in_trivial_patch",
                [](ElementAggregation & self) { return self.GetElsInTrivialPatch(); },
               "BitArray that is true for every element that is not part of a (non-trivial) patch")      
      .def_property_readonly("els_in_nontrivial_patch",
                [](ElementAggregation & self) { return self.GetElsInNontrivialPatch(); },
               "BitArray that is true for every element that is part of a (non-trivial) patch")      
      .def_property_readonly("element_to_patch",
                [](ElementAggregation & self) { 
                  py::list ret;
                  for (int i : self.GetElementToPatch())
                    ret.append(i);
                  return ret; 
                  },
               "vector mapping elements to (non-trivial) patches")      
      .def_property_readonly("facet_to_patch",
                [](ElementAggregation & self) { 
                  py::list ret;
                  for (int i : self.GetFacetToPatch())
                    ret.append(i);
                  return ret; 
                  },
               "vector mapping facets to (non-trivial) patches")      
    ;

  m.def("PatchwiseSolve", [](shared_ptr<ElementAggregation> elagg, 
                              shared_ptr<FESpace> fes,
                              shared_ptr<SumOfIntegrals> bf,
                              shared_ptr<SumOfIntegrals> lf,
                              int heapsize
                            )
  {
    VVector<double> hvec(fes->GetNDof());
    shared_ptr<BaseVector> vec = make_shared<VVector<double>>(hvec); 
    LocalHeap lh(heapsize, "Patchwisesolve-heap", true);
    PatchwiseSolve(elagg, fes, bf, lf, vec, lh);
    return vec;
  },
    py::arg("elagg"),
    py::arg("fes"),
    py::arg("bf"),
    py::arg("lf"),
    py::arg("heapsize") = 1000000, docu_string(R"raw_string(
Solve patch-wise problem based on the patches provided by the element aggregation input.

Parameters

elagg: ElementAggregatetion
  The instance defining the patches

fes : FESpace
  The finite element space on which the local solve is performed.

bf : SumOfIntegrals
  Integrators defining the matrix problem.

lf : SumOfIntegrals
  Integrators defining the right-hand side.
)raw_string")
  );


  m.def("ExtensionEmbedding", [] (shared_ptr<ElementAggregation> elagg, 
                                 shared_ptr<FESpace> fes,
                                 shared_ptr<SumOfIntegrals> bf,
                                 int heapsize
                                )
        { 
          LocalHeap lh(heapsize, "ExtensionEmbedding-heap", true);
          return SetupExtensionEmbedding(elagg, fes, bf, lh); 
        },
        py::arg("elagg"),
        py::arg("fes"),
        py::arg("bf"),
        py::arg("heapsize") = 1000000,
        docu_string(R"raw_string(
Computes the embedding matrix for an extension minimizing a described energy on 
patched (followed by averaging if some dofs are shared by multiple patches)

Parameters:

elagg : ElementAggregation
  ElementAggregation instace defining the patches which are aggregated into a single element

fes : ngsolve.FESpace
  The finite element space which is aggregated. 

bf : ngsolve.SumOfIntegrals
  The bilinear form describing the energy to be minized

)raw_string")
        );

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
          if ((!bb.is_none()) && py::extract<PyBA> (bb).check())
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
        py::arg("b")=py::none(),
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

  m.def("GetElementsWithSharedVertex",
      [] (shared_ptr<MeshAccess> ma,
         shared_ptr<BitArray> a,
         int heapsize)
      {
          LocalHeap lh (heapsize, "GetElementsWithNeighborFacets-heap", true);
          return GetElementsWithSharedVertex(ma,a,lh);
      } ,
      py::arg("mesh"),
      py::arg("a"),
      py::arg("heapsize") = 1000000);

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


  py::class_<MultiLevelsetCutInformation, shared_ptr<MultiLevelsetCutInformation>>
    (m, "MultiLevelsetCutInfo",R"raw(
A minimal version of a CutInfo that allows for several levelsets and a list of tuples of domain_types.
)raw")
    .def("__init__",  [] (MultiLevelsetCutInformation *instance,
                          shared_ptr<MeshAccess> ma,
                          py::object lsets_in)
         {
           py::extract<py::list> lsets_(lsets_in);
           if (!lsets_.check())
             throw Exception("levelset not compatible.");
           auto lsets = lsets_();
           for (int i = 0; i < py::len(lsets); i++)
             if (!(py::extract<shared_ptr<GridFunction>>(lsets[i]).check()))
               throw Exception("all lsets need to be GridFunctions!");
           Array<shared_ptr<GridFunction>> lset_a = makeCArray<shared_ptr<GridFunction>> (lsets);
           Array<shared_ptr<GridFunction>> lset_b;
           for (int i = 0; i < py::len(lsets); i++)
           {
            lset_b.Append(CreateGridFunction(lset_a[i]->GetFESpace(), "lset_p1", Flags()));
            lset_b[i]->Update();
            lset_b[i]->GetVectorPtr()->Set(1.0, lset_a[i]->GetVector());
           }
           new (instance) MultiLevelsetCutInformation (ma, lset_b);
         },
         py::arg("mesh"),
         py::arg("levelset"),
         docu_string(R"raw_string(
Creates a MultiLevelsetCutInfo based on a mesh and a tuple of levelsets.

Parameters

mesh : 
  mesh

levelsets : tuple(ngsolve.GridFunction)
  tuple of GridFunctions w.r.t. which elements are marked 
)raw_string")
      )
    .def("Mesh", [](MultiLevelsetCutInformation & self)
         {
           return self.GetMesh();
         },docu_string(R"raw_string(
Returns mesh of CutInfo.
)raw_string")
      )
    .def("Update", [] (MultiLevelsetCutInformation & self,
                       py::object lsets_in,
                       int heapsize)
         {
           LocalHeap lh (heapsize, "MultiLevelsetCutInfo-heap", true);

           py::extract<py::list> lsets_(lsets_in);
           if (!lsets_.check())
             throw Exception("levelset not compatible.");
           auto lsets = lsets_();
           for (int i = 0; i < py::len(lsets); i++)
             if (!(py::extract<shared_ptr<GridFunction>>(lsets[i]).check()))
               throw Exception("all lsets need to be GridFunctions!");
           if (py::len(lsets) != self.GetLen())
             throw Exception("New levelset tuple must have the same length as the original!");

           Array<shared_ptr<GridFunction>> lsets_a = makeCArray<shared_ptr<GridFunction>> (lsets);
           self.Update(lsets_a, lh);
         },
         py::arg("levelsets"),
         py::arg("heapsize") = 1000000,
         docu_string(R"raw_string(
Updates the tuple of levelsets behind the MultiLevelsetCutInfo and 
recomputes any element marker arrays which have been created with this
instance.

Parameters

levelsets : tuple(ngsolve.GridFunction)
  tuple of GridFunctions w.r.t. which elements are marked.

heapsize : int = 1000000
  heapsize of local computations)raw_string")
        )
    .def("GetElementsOfType", [](MultiLevelsetCutInformation & self,
                                 py::object dt_in,
                                 VorB vb,
                                 int heapsize)
         {

           LocalHeap lh (heapsize, "MultiLevelsetCutInfo-heap", true);

           if (py::isinstance<py::tuple>(dt_in))
           {
             if (py::len(dt_in) != self.GetLen())
               throw Exception("Number of domains does not match number of levelsets");
             Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (py::extract<py::tuple>(dt_in)());
             Array<Array<DOMAIN_TYPE>> dts(1);
             dts[0] = dts_;
             return self.GetElementsOfDomainType(dts, vb, lh);
           }           

           py::list dts_list;
           if (py::hasattr(dt_in, "as_list") && py::isinstance<py::list>(dt_in.attr("as_list")))
             dts_list = dt_in.attr("as_list");      
           else if (py::isinstance<py::list>(dt_in))
             dts_list = dt_in;
           else
             throw Exception("domain_type is neither a tuple nor a list nor a DomainTypeArray.");

           Array<Array<DOMAIN_TYPE>> cdts_aa(py::len(dts_list));
           int common_length = -1; //not a list
           for (int i = 0; i < py::len(dts_list); i++)
           {
             auto dta(dts_list[i]);

             // Check that input is valid
             if (!py::isinstance<py::tuple>(dta))
               throw Exception("domain_type arrays are incompatible. Maybe you used a list instead of a tuple?");
             else
             {
               if ((i>0) && (common_length != py::len(dta)))
                 throw Exception("domain_type arrays have different length");
               else
                 common_length = py::len(dta);
             }

             // Valid input. Add domain to pass to self.GetElementsOfDomainType
             cdts_aa[i] = makeCArray<DOMAIN_TYPE> (dta);
           }
           if (common_length != self.GetLen())
             throw Exception("Number of domains does not match number of levelsets");
           
           return self.GetElementsOfDomainType(cdts_aa, vb, lh);
         },
         py::arg("domain_type"),
         py::arg("VOL_or_BND") = VOL,
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Returns BitArray that is true for every element that has the 
corresponding domain type. This BitArray remains attached to the mlci class
instance and is updated on mlci.Update(lsets).

Parameters

domain_type : {tuple(ENUM), list(tuple(ENUM)), DomainTypeArray}
  Description of the domain.

heapsize : int = 1000000
  heapsize of local computations.
)raw_string")
    )
    .def("GetElementsWithContribution", [](MultiLevelsetCutInformation & self,
                                           py::object dt_in,
                                           VorB vb,
                                           int heapsize)
         {

           LocalHeap lh (heapsize, "MultiLevelsetCutInfo-heap", true);

           if (py::isinstance<py::tuple>(dt_in))
           {
             if (py::len(dt_in) != self.GetLen())
               throw Exception("Number of domains does not match number of levelsets");
             Array<DOMAIN_TYPE> dts_ = makeCArray<DOMAIN_TYPE> (py::extract<py::tuple>(dt_in)());
             Array<Array<DOMAIN_TYPE>> dts(1);
             dts[0] = dts_;
             return self.GetElementsWithContribution(dts, vb, lh);
           }

           py::list dts_list;
           if (py::hasattr(dt_in, "as_list") && py::isinstance<py::list>(dt_in.attr("as_list")))
             dts_list = dt_in.attr("as_list");      
           else if (py::isinstance<py::list>(dt_in))
             dts_list = dt_in;
           else
             throw Exception("domain_type is neither a tuple nor a list nor a DomainTypeArray.");

           Array<Array<DOMAIN_TYPE>> cdts_aa(py::len(dts_list));
           int common_length = -1; //not a list
           for (int i = 0; i < py::len(dts_list); i++)
           {
             auto dta(dts_list[i]);

             // Check valid input
             if (!py::isinstance<py::tuple>(dta))
               throw Exception("domain_type arrays are incompatible. Maybe you used a list instead of a tuple?");
             else
             {
               if ((i>0) && (common_length != py::len(dta)))
                 throw Exception("domain_type arrays have different length");
               else
                 common_length = py::len(dta);
             }

             // Valid input: Add domain to pass to self.GetElementsWithContribution
             cdts_aa[i] = makeCArray<DOMAIN_TYPE> (dta);
           }
           if (common_length != self.GetLen())
             throw Exception("Number of domains does not match number of levelsets");
           
           return self.GetElementsWithContribution(cdts_aa, vb, lh);
         },
         py::arg("domain_type"),
         py::arg("VOL_or_BND") = VOL,
         py::arg("heapsize") = 1000000,docu_string(R"raw_string(
Returns BitArray that is true for every element that has the 
a contribution to the corresponding level set domain. This BitArray 
remains attached to the mlci class instance and is updated on 
mlci.Update(lsets).

Parameters

domain_type : {tuple(ENUM), list(tuple(ENUM)), DomainTypeArray}
  Description of the domain.

heapsize : int = 1000000
  heapsize of local computations.
)raw_string"));


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
          if ((!acutinfo.is_none()) && py::extract<PyCI> (acutinfo).check())
            cutinfo = py::extract<PyCI>(acutinfo)();
          if ((!acutinfo.is_none()) && py::extract<PyCF> (acutinfo).check())
            cf_lset = py::extract<PyCF>(acutinfo)();
          if ((!alset.is_none()) && py::extract<PyCF> (alset).check())
            cf_lset = py::extract<PyCF>(alset)();


          Flags flags = py::extract<Flags> (bpflags)();

          if (basefes->IsComplex())
            flags.SetFlag("complex");

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
          ret->Update();
          return ret;
        },
        py::arg("basefes"),
        py::arg("cutinfo")=py::none(),
        py::arg("lset")=py::none(),
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

  m.def("SymbolicCutBFI", [](py::dict lsetdom,
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
          if ((!definedon.is_none() ) && defon_region.check())
            vb = VorB(py::extract<Region>(definedon)());

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

          VorB element_vb;
          if (element_boundary) element_vb = BND;
          else element_vb = VOL;

          shared_ptr<LevelsetIntegrationDomain> lsetintdom = PyDict2LevelsetIntegrationDomain(lsetdom);
          shared_ptr<BilinearFormIntegrator> bfi;
          if (!has_other && !skeleton)
          {
            bfi  = make_shared<SymbolicCutBilinearFormIntegrator> (*lsetintdom, cf, vb, element_vb);
          }
          else
          {
            if (lsetintdom->GetTimeIntegrationOrder() >= 0)
              throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for time_order >= 0..");
            if (vb == BND)
              throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for boundaries..");
            bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (*lsetintdom, cf);
          }

          if (!definedon.is_none() )
          {
            if (py::extract<py::list> (definedon).check())
              bfi -> SetDefinedOn (makeCArray<int> (definedon));
            
            if (defon_region.check())
            {
              cout << IM(3) << "definedon = " << defon_region().Mask() << endl;
              bfi->SetDefinedOn(defon_region().Mask());
            }
          }

          if (! definedonelem.is_none() )
            bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! deformation.is_none() )
            bfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyBFI(bfi);
        },
        py::arg("levelset_domain"),
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=false,
        py::arg("skeleton")=false,
        py::arg("definedon")=py::none(),
        py::arg("definedonelements")=py::none(),
        py::arg("deformation")=py::none(),
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
            cout << IM(2) << " no Other() used?!" << endl;

  
          shared_ptr<BilinearFormIntegrator> bfi;
          if (skeleton)
          {
            auto bfime = make_shared<SymbolicFacetBilinearFormIntegrator2> (cf);
            bfime->SetTimeIntegrationOrder(time_order);
            bfi = bfime;
          }
          else
          {
            // throw Exception("Patch facet blf not implemented yet: TODO(2)!");
            auto bfime = make_shared<SymbolicFacetPatchBilinearFormIntegrator> (cf);
            bfime->SetTimeIntegrationOrder(time_order);
            bfi = bfime;
          }

          if (! definedonelem.is_none())
            bfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! deformation.is_none())
            bfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyBFI(bfi);
        },
        py::arg("form"),
        py::arg("force_intorder")=-1,
        py::arg("time_order")=-1,
        py::arg("skeleton") = true,
        py::arg("definedonelements")=py::none(),
        py::arg("deformation")=py::none(),
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

  
  m.def("SymbolicCutLFI", [](py::dict lsetdom,
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
          if ((!definedon.is_none()) && defon_region.check())
            vb = VorB(defon_region());

          // if (vb == BND)
          //   throw Exception("Symbolic cuts not yet (tested) for boundaries..");

          if (element_boundary || skeleton)
            throw Exception("No Facet LFI with Symbolic cuts..");

          shared_ptr<LevelsetIntegrationDomain> lsetintdom = PyDict2LevelsetIntegrationDomain(lsetdom);
          auto lfi  = make_shared<SymbolicCutLinearFormIntegrator> (*lsetintdom, cf, vb);

          if ((!definedon.is_none()) && py::extract<py::list> (definedon).check())
            lfi -> SetDefinedOn (makeCArray<int> (definedon));

          if (defon_region.check())
          {
            cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
            lfi->SetDefinedOn(defon_region().Mask());
          }

          if (! definedonelem.is_none())
            lfi -> SetDefinedOnElements (py::extract<PyBA>(definedonelem)());

          if (! deformation.is_none())
            lfi->SetDeformation(py::extract<PyGF>(deformation)());

          return PyLFI(lfi);
        },
        py::arg("levelset_domain"),
        py::arg("form"),
        py::arg("VOL_or_BND")=VOL,
        py::arg("element_boundary")=py::bool_(false),
        py::arg("skeleton")=py::bool_(false),
        py::arg("definedon")=py::none(),
        py::arg("definedonelements")=py::none(),
        py::arg("deformation")=py::none(),
        docu_string(R"raw_string(
see documentation of SymbolicLFI (which is a wrapper))raw_string")
    );

  typedef shared_ptr<ProxyFunction> PyProxyFunction;
  m.def("dn", [] (const PyProxyFunction self, int order, py::object comp, bool hdiv)
        {

          const int dim_space = self->GetFESpace()->GetSpatialDimension();
          
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
          else if (py::extract<py::list> (comp).check())
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
