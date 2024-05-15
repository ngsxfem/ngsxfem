#include <python_ngstd.hpp>

/// from ngsolve
#include <solve.hpp>
#include <comp.hpp>
#include <fem.hpp>

#include "../cutint/straightcutrule.hpp"
#include "../cutint/xintegration.hpp"
#include "../cutint/mlsetintegration.hpp"
#include "../cutint/cutintegral.hpp"

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
           py::object gfdeformation,
           py::object ip_container,
           bool element_wise,
           int heapsize)
        {
          static Timer timer("IntegrateX");
          RegionTimer reg (timer);
          shared_ptr<py::list> ip_cont = nullptr;
          py::extract<py::object> ip_cont_(ip_container);
          if (!ip_cont_().is_none())
          {
            py::extract<py::list> ip_cont_as_list(ip_container);
            if (ip_cont_as_list.check())
              ip_cont = make_shared<py::list>(ip_cont_as_list());
          }
          shared_ptr<GridFunction> deformation = nullptr;
          py::extract<py::object> deformation_(gfdeformation);
          if (!deformation_().is_none())
          {
            py::extract<shared_ptr<GridFunction>> deformation_as_gf(gfdeformation);
            if (deformation_as_gf.check())
              deformation = deformation_as_gf();
          }

          shared_ptr<LevelsetIntegrationDomain> lsetintdom = PyDict2LevelsetIntegrationDomain(lsetdom);
          bool space_time = lsetintdom->GetTimeIntegrationOrder() >= 0;
          LocalHeap lh(heapsize, "lh-IntegrateX");


          // int DIM = ma->GetDimension();

          int cfdim = cf->Dimension();
          if(element_wise && cfdim != 1)
            throw Exception("element_wise only implemented for 1 dimensional coefficientfunctions");

          MyMutex mutex_ip_cont;
          //MyLock lock_ip_cont(mutex_ip_cont);
					
          Vector<> element_sum(element_wise ? ma->GetNE(VOL) : 0);
          element_sum = 0.0;
					double sum = 0.0;
					if(globxvar.SIMD_EVAL && ip_cont==nullptr) {
            try {
						  ma->IterateElements
              (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
              {
                HeapReset hr(lh);
						  	auto & trafo1 = ma->GetTrafo (el, lh);
         		  	auto & trafo = trafo1.AddDeformation(deformation.get(), lh);

                const IntegrationRule *ns_ir;
						  	Array<double> ns_wei_arr;
						  	tie(ns_ir, ns_wei_arr) = CreateCutIntegrationRule(*lsetintdom, trafo, lh);

						  	if (ns_ir == nullptr)
						  		return;        

						  	SIMD_IntegrationRule simd_ir(*ns_ir, lh);
						  	FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(ns_wei_arr, lh);
  
						  	SIMD_BaseMappedIntegrationRule &simd_mir = trafo(simd_ir, lh);
						  	FlatMatrix<SIMD<double>> val(simd_mir.Size(), 1, lh);
  
						  	cf -> Evaluate(simd_mir, val);
  
						  	SIMD<double> lsum = 0.0;
						  	for (int i = 0; i < simd_mir.Size(); i++)
						  		lsum += simd_mir[i].GetMeasure()*simd_wei_arr[i]*val(i,0);

                AtomicAdd(sum,HSum(lsum));
                if (element_wise)
                  element_sum(el.Nr()) = HSum(lsum); // problem ?

						  });
						  py::object result;
              if (element_wise)
                result = py::cast(element_sum);
              else
  						  result = py::cast(ma->GetCommunicator().AllReduce(sum, NG_MPI_SUM));
	
  					  return result;
            } catch (ExceptionNOSIMD e) {
              cout << IM(6) << e.What() << endl
                   << "switching to non-SIMD evaluation" << endl;
            }

					}

					if(globxvar.SIMD_EVAL && ip_cont!=nullptr) 
            cout << IM(6) << "Switching to non-SIMD evaluation due to IntegrationPointContainer" << endl;
          
					ma->IterateElements
            (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
             {
               HeapReset hr(lh);
               auto & trafo1 = ma->GetTrafo (el, lh);
               auto & trafo = trafo1.AddDeformation(deformation.get(), lh);

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

                 if (ip_cont != nullptr)
                   for (int i = 0; i < mir.Size(); i++)
                   {
                     MyLock lock_ip_cont(mutex_ip_cont);
                     if (space_time)
                         ip_cont->append(tuple ( MeshPoint{mir[i].IP()(0), mir[i].IP()(1), mir[i].IP()(2),
                                                ma.get(), VOL, static_cast<int>(el.Nr())}, (*ir)[i].Weight()));
                     else
                        ip_cont->append(MeshPoint{mir[i].IP()(0), mir[i].IP()(1), mir[i].IP()(2),
                                               ma.get(), VOL, static_cast<int>(el.Nr())});
                   }
                 if (element_wise)
                   element_sum(el.Nr()) = lsum;
                 
                 AtomicAdd(sum,lsum);
               }
             });

					py::object result;
          if (element_wise) {
            result = py::cast(element_sum);
          } else {
            result = py::cast(ma->GetCommunicator().AllReduce(sum, NG_MPI_SUM));
          }
          return result;
        },
        py::arg("levelset_domain"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("deformation")=py::none(),
        py::arg("ip_container")=py::none(),
        py::arg("element_wise")=false,
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

deformation : gridfunction (or None)
  deformation of the mesh

ip_container : list (or None)
  a list to store integration points (for debugging or visualization purposes)

element_wise : bool
  result will return the integral w.r.t. each element individually.

heapsize : int
  heapsize for local computations.
)raw_string"));

    m.def("IntegrationPointExtrema",
        [](py::dict lsetdom,
           shared_ptr<MeshAccess> ma,
           PyCF cf,
           int heapsize)
        {
          static Timer timer("IntegrationPointExtrema"); RegionTimer reg (timer);

          shared_ptr<LevelsetIntegrationDomain> lsetintdom = PyDict2LevelsetIntegrationDomain(lsetdom);
          // bool space_time = lsetintdom->GetTimeIntegrationOrder() >= 0;
          LocalHeap lh(heapsize, "lh-IntegrationPointExtrema");

          double min = 1e99;
          double max = -1e99;

          // int DIM = ma->GetDimension();
          // int cfdim = cf->Dimension();

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

                 double lmin = 1e99;
                 double lmax = -1e99;
                 for (int i = 0; i < mir.Size(); i++)
                 {
                     if (val(i,0) > lmax) 
                       lmax = val(i,0);
                     if (val(i,0) < lmin) 
                       lmin = val(i,0);
                 }

                 AtomicMax(max,lmax);
                 AtomicMin(min,lmin);
               }
             });
          py::list minmax;
          minmax.append(min);
          minmax.append(max);
          return py::tuple(minmax);
        },
        py::arg("levelset_domain"),
        py::arg("mesh"),
        py::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
        py::arg("heapsize")=1000000,
        docu_string(R"raw_string(
Determine minimum and maximum on integration points on a level set domain. The sampling uses the same
integration rule as in Integrate and is determined by 'order' w.r.t. a (multi-)linear
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


  py::class_<CutDifferentialSymbol,DifferentialSymbol>(m, "CutDifferentialSymbol",
docu_string(R"raw_string(
CutDifferentialSymbol that allows to formulate linear, bilinear forms and integrals on
level set domains in an intuitive form:

Example use case:

  dCut = CutDifferentialSymbol(VOL)
  dx = dCut(lset,NEG)
  a = BilinearForm(...)
  a += u * v * dx

Note that the most important options are set in the second line when the basic
CutDifferentialSymbol is further specified.
)raw_string")  
  )
    .def(py::init<>(), docu_string(R"raw_string(
Constructor of CutDifferentialSymbol.

  Argument: none
)raw_string"))
    .def(py::init<VorB>(), docu_string(R"raw_string(
Constructor of CutDifferentialSymbol.

  Argument: VOL_or_BND (boundary or volume form?).
)raw_string"))
    .def("__call__", [](CutDifferentialSymbol & self,
                        py::dict lsetdom,
                        optional<variant<Region,string>> definedon,
                        VorB vb, 
                        bool element_boundary,
                        VorB element_vb, bool skeleton,
                        shared_ptr<GridFunction> deformation,
                        shared_ptr<BitArray> definedonelements)
         {
           if (element_boundary) element_vb = BND;
           auto dx = CutDifferentialSymbol(PyDict2LevelsetIntegrationDomain(lsetdom), 
                                           vb, element_vb, skeleton);
           if (definedon)
             {
               if (auto definedon_region = get_if<Region>(&*definedon); definedon_region)
                 {
                   dx.definedon = definedon_region->Mask();
                   dx.vb = VorB(*definedon_region);
                 }
               if (auto definedon_string = get_if<string>(&*definedon); definedon_string)
                 dx.definedon = *definedon_string;
             }
           dx.deformation = deformation;
           dx.definedonelements = definedonelements;
           return dx;
         },
         py::arg("levelset_domain"),
         py::arg("definedon")=nullptr,
         py::arg("vb")=VOL,
         py::arg("element_boundary")=false,
         py::arg("element_vb")=VOL,
         py::arg("skeleton")=false,
         py::arg("deformation")=nullptr,
         py::arg("definedonelements")=nullptr,
docu_string(R"raw_string(
The call of a CutDifferentialSymbol allows to specify what is needed to specify the 
integration domain. It returns a new CutDifferentialSymbol.

Parameters:

levelset_domain (dict) : specifies the level set domain.
definedon (Region or Array) : specifies on which part of the mesh (in terms of regions)
  the current form shall be defined.
vb (VOL/BND/BBND/BBBND) : Where does the integral take place from point of view 
  of the mesh.
element_boundary (bool) : Does the integral take place on the boundary of an element-
element_vb (VOL/BND/BBND/BBBND) : Where does the integral take place from point of view
  of an element.
skeleton (bool) : is it an integral on facets (the skeleton)?
deformation (GridFunction) : which mesh deformation shall be applied (default : None)
definedonelements (BitArray) : Set of elements or facets where the integral shall be
  defined.
)raw_string"))
    .def("__rmul__", [](CutDifferentialSymbol & self, double x)
    {
      return CutDifferentialSymbol(self, x );
    })

    .def("order", [](CutDifferentialSymbol & self, int order)
    {
      auto _cds = CutDifferentialSymbol(self);
      _cds.lsetintdom->SetIntegrationOrder(order);
      return _cds;
    },
    py::arg("order"))
    .def_property("vb",
                  [](CutDifferentialSymbol & self) { return self.vb; },
                  [](CutDifferentialSymbol & self, VorB vb) { self.vb = vb; return self.vb;},
                  "Volume of boundary?")
    ;
    
  py::class_<FacetPatchDifferentialSymbol,DifferentialSymbol>(m, "FacetPatchDifferentialSymbol",
docu_string(R"raw_string(
FacetPatchDifferentialSymbol that allows to formulate integrals on facet patches.
Example use case:

  dFacetPatch = FacetPatchDifferentialSymbol(VOL)
  dw = dFacetPatch(definedonelements = ...)
  a = BilinearForm(...)
  a += (u-u.Other()) * (v-v.Other()) * dw

)raw_string")  
  )
    .def(py::init<VorB>(), docu_string(R"raw_string(
Constructor of FacetPatchDifferentialSymbol.

  Argument: VOL_or_BND (boundary or volume form?).
)raw_string"))
    .def("__call__", [](FacetPatchDifferentialSymbol & self,
                        optional<variant<Region,string>> definedon,
                        bool element_boundary,
                        VorB element_vb, bool skeleton,
                        shared_ptr<GridFunction> deformation,
                        shared_ptr<BitArray> definedonelements,
                        int time_order, 
                        optional<double> tref)
         {
           if (element_boundary) element_vb = BND;
           auto dx = FacetPatchDifferentialSymbol(self.vb, element_vb, skeleton,time_order,tref);
           if (definedon)
             {
               if (auto definedon_region = get_if<Region>(&*definedon); definedon_region)
                 {
                   dx.definedon = definedon_region->Mask();
                   dx.vb = VorB(*definedon_region);
                 }
               if (auto definedon_string = get_if<string>(&*definedon); definedon_string)
                 dx.definedon = *definedon_string;
             }
           dx.deformation = deformation;
           dx.definedonelements = definedonelements;
           return dx;
         },
         py::arg("definedon")=nullptr,
         py::arg("element_boundary")=false,
         py::arg("element_vb")=VOL,
         py::arg("skeleton")=false,
         py::arg("deformation")=nullptr,
         py::arg("definedonelements")=nullptr,
         py::arg("time_order")=-1,
         py::arg("tref")=nullopt,
docu_string(R"raw_string(
The call of a FacetPatchDifferentialSymbol allows to specify what is needed to specify the 
integration domain of an integral that runs over the volume patch of each facet. 
It returns a new CutDifferentialSymbol.

Parameters:

definedon (Region or Array) : specifies on which part of the mesh (in terms of regions)
  the current form shall be defined.
element_boundary (bool) : Does the integral take place on the boundary of an element-
element_vb (VOL/BND/BBND/BBBND) : Where does the integral take place from point of view
  of an element.
skeleton (bool) : is it an integral on facets (the skeleton)?
deformation (GridFunction) : which mesh deformation shall be applied (default : None)
definedonelements (BitArray) : Set of elements or facets where the integral shall be
  defined.
time_order (int) : integration order in time (for space-time) (default : -1).
tref (float) : turn space integral into space-time integral with fixed time tref.
)raw_string"         
         ))
    .def("__rmul__", [](FacetPatchDifferentialSymbol & self, double x)
    {
      return FacetPatchDifferentialSymbol(self, x );
    })
    ;



  
}
