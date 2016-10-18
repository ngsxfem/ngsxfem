#ifdef NGSX_PYTHON
//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../xfem/xFESpace.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"
#include "../stokes/xstokesspace.hpp"
#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../utils/error.hpp"

#include "../cutint/straightcutrule.hpp"

//using namespace ngcomp;

void ExportNgsx() 
{
  std::string nested_name = "xfem";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".xfem");
  
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

  cout << "exporting xfem as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("xfem") = module ;

  bp::scope local_scope(module);

  typedef PyWrapper<FESpace> PyFES;
  typedef PyWrapper<CoefficientFunction> PyCF;
  typedef GridFunction GF;
  typedef PyWrapperDerived<GF, CoefficientFunction> PyGF;
  
  bp::enum_<DOMAIN_TYPE>("DOMAIN_TYPE")
    .value("POS", POS)
    .value("NEG", NEG)
    .value("IF", IF)
    .export_values()
    ;

  // typedef PyWrapperDerived<CompoundFESpace, FESpace> PyCompFES;
  
  typedef PyWrapperDerived<XFESpace, FESpace> PyXFES;
  typedef PyWrapperDerived<XStdFESpace, FESpace> PyXStdFES;
  typedef PyWrapperDerived<XStokesFESpace, FESpace> PyXStokesFES;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XStdFESpace>);
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XFESpace>);
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XStokesFESpace>);
  
  bp::def("CastToXFESpace", FunctionPointer( [] (PyFES fes) -> PyXFES { return PyXFES(dynamic_pointer_cast<XFESpace>(fes.Get())); } ) );
  bp::def("CastToXStdFESpace", FunctionPointer( [] (PyFES fes) -> PyXStdFES { return PyXStdFES(dynamic_pointer_cast<XStdFESpace>(fes.Get())); } ) );
  bp::def("CastToXStokesFESpace", FunctionPointer( [] (PyFES fes) -> PyXStokesFES { return PyXStokesFES(dynamic_pointer_cast<XStokesFESpace>(fes.Get())); } ) );

  bp::def("XToNegPos", FunctionPointer( [] (PyGF gfx, PyGF gfnegpos) { XFESpace::XToNegPos(gfx.Get(),gfnegpos.Get()); } ) );
  

  bp::class_<PyXFES, bp::bases<PyFES>>
    ("XFESpace", bp::no_init)
    .def("SetLevelSet", FunctionPointer ([](PyXFES self, PyCF cf) 
                                         { self.Get()->SetLevelSet(cf.Get()); }),
         "Update information on level set function")
    .def("SetLevelSet", FunctionPointer ([](PyXFES self, PyGF gf) 
                                         { self.Get()->SetLevelSet(gf.Get()); }),
         "Update information on level set function")
    .def("SetBaseFESpace", FunctionPointer ([](PyXFES self, PyFES fes) 
                                            { self.Get()->SetBaseFESpace(fes.Get()); }),
         "Update information on base FESpace")
    .def("BaseDofOfXDof", FunctionPointer ([](PyXFES self, int i) 
                                         { return self.Get()->GetBaseDofOfXDof(i); }),
         "get corresponding dof of base FESpace")
    .def("GetNVertexDofs", FunctionPointer ([](PyXFES self) 
                                            { return self.Get()->GetNVertexDof(); }),
         "get number of x dofs at vertices")
    .def("CutElements", FunctionPointer ([](PyXFES self) 
                                         { return self.Get()->CutElements(); }),
         "get BitArray of cut elements")
    .def("CutSurfaceElements", FunctionPointer ([](PyXFES self) 
                                         { return self.Get()->CutSurfaceElements(); }),
         "get BitArray of cut surface elements")
    .def("GetDomainOfDof", FunctionPointer ([](PyXFES self, int i) 
                                         { return self.Get()->GetDomainOfDof(i); }),
         "get domain_type of degree of freedom")
    .def("GetDomainOfElement", FunctionPointer ([](PyXFES self, int i) 
                                         { return self.Get()->GetDomainOfElement(i); }),
         "get domain_type of element")
    .def("GetDomainNrs",  FunctionPointer( [] (PyXFES self, int elnr) {
               Array<DOMAIN_TYPE> domnums;
               self.Get()->GetDomainNrs( elnr, domnums );
               return domnums;
            }))
    ;


  bp::class_<PyXStdFES, bp::bases<PyFES>>
    ("XStdFESpace", bp::no_init)
    .add_property("XFESpace", FunctionPointer ([](const PyXStdFES self) 
                                               { return PyXFES(dynamic_pointer_cast<XFESpace> ((*self.Get())[1])); }
                    ),
         "return XFESpace part of XStdFESpace")
    .add_property("StdFESpace", FunctionPointer ([](const PyXStdFES self) 
                                                 { return PyFES((*self.Get())[0]); }
                    ),
         "return 'standard' FESpace part of XStdFESpace")
    ;

  bp::class_<PyXStokesFES, bp::bases<PyFES>>
    ("XStokesFESpace", bp::no_init)
    .def("SetLevelSet", FunctionPointer ([](PyXStokesFES & self, PyCF cf) 
                                         { self.Get()->SetLevelSet(cf.Get()); }),
         "Update information on level set function")
    .def("SetLevelSet", FunctionPointer ([](PyXStokesFES & self, PyGF gf) 
                                         { self.Get()->SetLevelSet(gf.Get()); }),
         "Update information on level set function")
    ;

  bp::def("InterpolateToP1", FunctionPointer( [] (PyGF gf_ho, PyGF gf_p1, int heapsize)
                                              {
                                                InterpolateP1 interpol(gf_ho.Get(), gf_p1.Get());
                                                LocalHeap lh (heapsize, "InterpolateP1-Heap");
                                                interpol.Do(lh);
                                              } ),
          (bp::arg("gf_ho")=NULL,bp::arg("gf_p1")=NULL,bp::arg("heapsize")=1000000))
    ;

  bp::def("InterpolateToP1", FunctionPointer( [] (PyCF coef, PyGF gf_p1, int heapsize)
                                              {
                                                InterpolateP1 interpol(coef.Get(), gf_p1.Get());
                                                LocalHeap lh (heapsize, "InterpolateP1-Heap");
                                                interpol.Do(lh);
                                              } ),
          (bp::arg("coef"),bp::arg("heapsize")=1000000))
    ;

  bp::class_<StatisticContainer, shared_ptr<StatisticContainer>,  boost::noncopyable>("StatisticContainer", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([]()
                           {
                             return make_shared<StatisticContainer> ();
                           })))
    .def("Print", FunctionPointer ([](StatisticContainer & self, string label, string select)
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
                           }),
         (bp::arg("self")=NULL,bp::arg("label")="something",bp::arg("select")="all")
      )
    ;

  bp::def("CalcMaxDistance", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, int heapsize)
                                              {
                                                StatisticContainer dummy;
                                                LocalHeap lh (heapsize, "CalcDistance-Heap");
                                                if (lset_p1->GetMeshAccess()->GetDimension()==2)
                                                  CalcDistances<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  dummy, lh, -1.0, false);
                                                else
                                                  CalcDistances<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  dummy, lh, -1.0, false);
                                                return (double) dummy.ErrorMaxNorm[dummy.ErrorMaxNorm.Size()-1];
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("heapsize")=1000000))
    ;

  bp::def("CalcDistances", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
                                              {
                                                LocalHeap lh (heapsize, "CalcDistance-Heap");
                                                if (lset_p1.Get()->GetMeshAccess()->GetDimension()==2)
                                                  CalcDistances<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  stats, lh, refine_threshold, absolute);
                                                else
                                                  CalcDistances<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(),  stats, lh, refine_threshold, absolute);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("stats")=NULL,bp::arg("heapsize")=1000000,bp::arg("refine_threshold")=-1.0,bp::arg("absolute")=false))
    ;

  bp::def("CalcDeformationError", FunctionPointer( [] (PyCF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, StatisticContainer & stats, double lower, double upper, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "CalcDeformationError-Heap");
                                                if (lset_p1.Get()->GetMeshAccess()->GetDimension()==2)
                                                  CalcDeformationError<2>(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), stats, lh, lower, upper);
                                                else
                                                  CalcDeformationError<3>(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), stats, lh, lower, upper);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("qn")=NULL,bp::arg("stats")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("heapsize")=1000000))
    ;

  bp::def("ProjectShift", FunctionPointer( [] (PyGF lset_ho, PyGF lset_p1, PyGF deform, PyCF qn, double lower, double upper, double threshold, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "ProjectShift-Heap");
                                                ProjectShift(lset_ho.Get(), lset_p1.Get(), deform.Get(), qn.Get(), lower, upper, threshold, lh);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("qn")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("threshold")=1.0,bp::arg("heapsize")=1000000))
    ;

// ProjectShift

  
  bp::def("RefineAtLevelSet", FunctionPointer( [] (PyGF lset_p1, double lower, double upper, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
                                                RefineAtLevelSet(lset_p1.Get(), lower, upper, lh);
                                              } ),
          (bp::arg("lset_p1")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("heapsize")=1000000))
    ;


  bp::def("CalcTraceDiff", FunctionPointer( [] (PyGF gf, PyCF coef, int intorder, int heapsize)
                                              {
                                                Array<double> errors;
                                                LocalHeap lh (heapsize, "CalcTraceDiff-Heap");
                                                if (gf.Get()->GetMeshAccess()->GetDimension() == 2)
                                                  CalcTraceDiff<2>(gf.Get(),coef.Get(),intorder,errors,lh);
                                                else 
                                                  CalcTraceDiff<3>(gf.Get(),coef.Get(),intorder,errors,lh);
                                                return errors;
                                              } ),
          (bp::arg("gf"),bp::arg("coef"),bp::arg("intorder")=6,bp::arg("heapsize")=1000000))
    ;


  bp::def("RefineAtLevelSet", FunctionPointer( [] (PyGF gf, double lower_lset_bound, double upper_lset_bound, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
                                                RefineAtLevelSet(gf.Get(),lower_lset_bound,upper_lset_bound,lh);
                                              } ),
          (bp::arg("gf"),bp::arg("lower_lset_bound")=0.0,bp::arg("upper_lset_bound")=0.0,bp::arg("heapsize")=10000000))
    ;


  bp::def("IntegrateX", 
          FunctionPointer([](PyCF lset,
                             shared_ptr<MeshAccess> ma, 
                             PyCF cf_neg,
                             PyCF cf_pos,
                             PyCF cf_interface,
                             int order, int subdivlvl, bp::dict domains, int heapsize)
                          {
                            static Timer timer ("IntegrateX");
                            static Timer timercutgeom ("IntegrateX::MakeCutGeom");
                            static Timer timerevalcoef ("IntegrateX::EvalCoef");
                            RegionTimer reg (timer);
                            LocalHeap lh(heapsize, "lh-Integrate");
                            
                            Flags flags = bp::extract<Flags> (domains)();
                            
                            Array<bool> tointon(3); tointon = false;
                            Array<shared_ptr<CoefficientFunction>> cf(3); cf = nullptr;
                            if (flags.GetDefineFlag("negdomain")){
                              tointon[int(NEG)] = true;
                              if (cf_neg.Get() == nullptr)
                                throw Exception("no coef for neg domain given");
                              else
                                cf[int(NEG)] = cf_neg.Get();
                            }
                            if (flags.GetDefineFlag("posdomain")){
                              tointon[int(POS)] = true;
                              if (cf_pos.Get() == nullptr)
                                throw Exception("no coef for pos domain given");
                              else
                                cf[int(POS)] = cf_pos.Get();
                            }
                            if (flags.GetDefineFlag("interface")){
                              tointon[int(IF)] = true;
                              if (cf_interface.Get() == nullptr)
                                throw Exception("no coef for interface domain given");
                              else
                                cf[int(IF)] = cf_interface.Get();
                            }

                            Vector<> domain_sum(3); // [val_neg,val_pos,val_interface]
                            domain_sum = 0.0;
                            int DIM = ma->GetDimension();
                            ma->IterateElements
                              (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
                               {
                                 auto & trafo = ma->GetTrafo (el, lh);
                                 auto lset_eval
                                   = ScalarFieldEvaluator::Create(DIM,*(lset.Get()),trafo,lh);
                                 ELEMENT_TYPE eltype = el.GetType();
                                 timercutgeom.Start();
                                 shared_ptr<XLocalGeometryInformation> xgeom = nullptr;

                                 CompositeQuadratureRule<2> cquad2d;
                                 CompositeQuadratureRule<3> cquad3d;
                                 if (DIM == 2)
                                   xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT,
                                                                             *lset_eval, cquad2d, lh,
                                                                             order, 0, subdivlvl, 0);
                                 else
                                   xgeom = XLocalGeometryInformation::Create(eltype, ET_POINT,
                                                                             *lset_eval, cquad3d, lh,
                                                                             order, 0, subdivlvl, 0);
                                 DOMAIN_TYPE element_domain = xgeom->MakeQuadRule();
                                 timercutgeom.Stop();
                                 if (element_domain == IF)
                                 {
                                   for (auto domtype : {NEG,POS})
                                   {
                                     if( ! tointon[int(domtype)] ) continue;
                                     double hsum = 0.0;
                                     // double value = 0.0;

                                     if (DIM == 2)
                                     {
                                       const QuadratureRule<2> & domain_quad = cquad2d.GetRule(domtype);
                                       IntegrationRule ir_domain (domain_quad.Size(),lh);
                                       for (int i = 0; i < ir_domain.Size(); ++i)
                                         ir_domain[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
                                       
                                       MappedIntegrationRule<2,2> mir_domain(ir_domain, trafo,lh);
                                       FlatMatrix<> values(mir_domain.Size(), 1, lh);
                                       timerevalcoef.Start();
                                       cf[int(domtype)] -> Evaluate (mir_domain, values);
                                       timerevalcoef.Stop();

                                       for (int i = 0; i < domain_quad.Size(); ++i)
                                         hsum += values(i,0) * mir_domain[i].GetWeight(); 
                                     }
                                     else
                                     {
                                       const QuadratureRule<3> & domain_quad = cquad3d.GetRule(domtype);
                                       IntegrationRule ir_domain (domain_quad.Size(),lh);
                                       for (int i = 0; i < ir_domain.Size(); ++i)
                                         ir_domain[i] = IntegrationPoint (&domain_quad.points[i](0),domain_quad.weights[i]);
                                       
                                       MappedIntegrationRule<3,3> mir_domain(ir_domain, trafo,lh);
                                       FlatMatrix<> values(mir_domain.Size(), 1, lh);
                                       timerevalcoef.Start();
                                       cf[int(domtype)] -> Evaluate (mir_domain, values);
                                       timerevalcoef.Stop();

                                       for (int i = 0; i < domain_quad.Size(); ++i)
                                         hsum += values(i,0) * mir_domain[i].GetWeight(); 
                                     }


                                     double & rsum = domain_sum(int(domtype));
                                     AsAtomic(rsum) += hsum;
                                   }

                                   if (tointon[int(IF)])
                                   {
                                     double hsum_if = 0.0;
                                     double value_if = 0.0;
                                     if (DIM == 2)
                                     {
                                       const QuadratureRuleCoDim1<2> & interface_quad(cquad2d.GetInterfaceRule());
                                       IntegrationRule ir_interface (interface_quad.Size(),lh);
                                       for (int i = 0; i < ir_interface.Size(); ++i)
                                         ir_interface[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i]);
                                       
                                       MappedIntegrationRule<2,2> mir_interface(ir_interface, trafo,lh);
                                       FlatMatrix<> values(mir_interface.Size(), 1, lh);
                                       timerevalcoef.Start();
                                       cf[int(IF)] -> Evaluate (mir_interface, values);
                                       timerevalcoef.Stop();
                                       
                                       for (int i = 0; i < interface_quad.Size(); ++i)
                                       {
                                         MappedIntegrationPoint<2,2> & mip(mir_interface[i]);

                                         Mat<2,2> Finv = mip.GetJacobianInverse();
                                         const double absdet = mip.GetMeasure();

                                         Vec<2> nref = interface_quad.normals[i];
                                         Vec<2> normal = absdet * Trans(Finv) * nref ;
                                         double len = L2Norm(normal);
                                         const double weight = interface_quad.weights[i] * len;

                                         hsum_if += values(i,0) * weight; 
                                       }
                                     }
                                     else
                                     {
                                       
                                       const QuadratureRuleCoDim1<3> & interface_quad(cquad3d.GetInterfaceRule());
                                       IntegrationRule ir_interface (interface_quad.Size(),lh);
                                       for (int i = 0; i < ir_interface.Size(); ++i)
                                         ir_interface[i] = IntegrationPoint (&interface_quad.points[i](0),interface_quad.weights[i]);
                                       
                                       MappedIntegrationRule<3,3> mir_interface(ir_interface, trafo,lh);
                                       FlatMatrix<> values(mir_interface.Size(), 1, lh);
                                       timerevalcoef.Start();
                                       cf[int(IF)] -> Evaluate (mir_interface, values);
                                       timerevalcoef.Stop();
                                       
                                       for (int i = 0; i < interface_quad.Size(); ++i)
                                       {
                                         MappedIntegrationPoint<3,3> & mip(mir_interface[i]);

                                         Mat<3,3> Finv = mip.GetJacobianInverse();
                                         const double absdet = mip.GetMeasure();

                                         Vec<3> nref = interface_quad.normals[i];
                                         Vec<3> normal = absdet * Trans(Finv) * nref ;
                                         double len = L2Norm(normal);
                                         const double weight = interface_quad.weights[i] * len;
                                         
                                         hsum_if += values(i,0) * weight; 
                                       }
                                     }

                                     double & rsum_if = domain_sum(int(IF));
                                     AsAtomic(rsum_if) += hsum_if;
                                   }
                                 }
                                 else if( tointon[int(element_domain)] )
                                 {
                                   double hsum = 0.0;
                                   IntegrationRule ir(trafo.GetElementType(), order);
                                   BaseMappedIntegrationRule & mir = trafo(ir, lh);

                                   FlatMatrix<> values(ir.Size(), 1, lh);
                                   timerevalcoef.Start();
                                   cf[int(element_domain)] -> Evaluate (mir, values);
                                   timerevalcoef.Stop();
                                   for (int i = 0; i < values.Height(); i++)
                                     hsum += mir[i].GetWeight() * values(i,0);

                                   double & rsum = domain_sum(int(element_domain));
                                   AsAtomic(rsum) += hsum;
                                 }
                                   
                               });

                            bp::dict resdict;
                            if (tointon[int(NEG)])
                              resdict["negdomain"] = domain_sum(int(NEG));
                            if (tointon[int(POS)])
                              resdict["posdomain"] = domain_sum(int(POS));
                            if (tointon[int(IF)])
                              resdict["interface"] = domain_sum(int(IF));

                            return resdict;
                            // bp::object result;
                            // return  bp::list(bp::object(result_vec));
                          }),
          (bp::arg("lset"), bp::arg("mesh"), 
           bp::arg("cf_neg")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)), 
           bp::arg("cf_pos")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
           bp::arg("cf_interface")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)),
           bp::arg("order")=5, bp::arg("subdivlvl")=0, bp::arg("domains")=bp::dict(), bp::arg("heapsize")=1000000))
    ;


  typedef PyWrapper<BilinearFormIntegrator> PyBFI;
  typedef PyWrapper<LinearFormIntegrator> PyLFI;
  
  bp::def("SymbolicCutBFI", FunctionPointer
          ([](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              bp::object definedon)
           -> PyBFI
           {
             
             bp::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             if (vb == BND)
               throw Exception("Symbolic cuts not yet (tested) for boundaries..");

             // check for DG terms
             bool has_other = false;
             cf.Get()->TraverseTree ([&has_other] (CoefficientFunction & cf)
                               {
                                 if (dynamic_cast<ProxyFunction*> (&cf))
                                   if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                                     has_other = true;
                               });
               
             
             if (has_other || element_boundary || skeleton)
               throw Exception("No Facet BFI with Symbolic cuts..");
             
             shared_ptr<BilinearFormIntegrator> bfi
               = make_shared<SymbolicCutBilinearFormIntegrator> (lset.Get(), cf.Get(), dt, order, subdivlvl);
             
             if (bp::extract<bp::list> (definedon).check())
               bfi -> SetDefinedOn (makeCArray<int> (definedon));

             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 bfi->SetDefinedOn(defon_region().Mask());
               }
             
             return PyBFI(bfi);
           }),
          (bp::args("lset"),
           bp::args("domain_type")=NEG,
           bp::args("force_intorder")=-1,
           bp::args("subdivlvl")=0,
           bp::args("form"),
           bp::args("VOL_or_BND")=VOL,
           bp::args("element_boundary")=false,
           bp::args("skeleton")=false,
           bp::arg("definedon")=bp::object())
          );
  
  
  bp::def("SymbolicCutLFI", FunctionPointer
          ([](PyCF lset,
              DOMAIN_TYPE dt,
              int order,
              int subdivlvl,
              PyCF cf,
              VorB vb,
              bool element_boundary,
              bool skeleton,
              bp::object definedon)
           -> PyLFI
           {
             
             bp::extract<Region> defon_region(definedon);
             if (defon_region.check())
               vb = VorB(defon_region());

             if (vb == BND)
               throw Exception("Symbolic cuts not yet (tested) for boundaries..");
               
             if (element_boundary || skeleton)
               throw Exception("No Facet LFI with Symbolic cuts..");
             
             shared_ptr<LinearFormIntegrator> lfi
               = make_shared<SymbolicCutLinearFormIntegrator> (lset.Get(), cf.Get(), dt, order, subdivlvl);
             
             if (bp::extract<bp::list> (definedon).check())
               lfi -> SetDefinedOn (makeCArray<int> (definedon));

             if (defon_region.check())
               {
                 cout << IM(3) << "defineon = " << defon_region().Mask() << endl;
                 lfi->SetDefinedOn(defon_region().Mask());
               }
             
             return PyLFI(lfi);
           }),
          (bp::args("lset"),
           bp::args("domain_type")=NEG,
           bp::args("force_intorder")=-1,
           bp::args("subdivlvl")=0,
           bp::args("form"),
           bp::args("VOL_or_BND")=VOL,
           bp::args("element_boundary")=false,
           bp::args("skeleton")=false,
           bp::arg("definedon")=bp::object())
          );
  
  
  // bp::def("GFCoeff", FunctionPointer( [] (shared_ptr<GridFunction> in) { return dynamic_pointer_cast<CoefficientFunction>(make_shared<GridFunctionCoefficientFunction>(in)); } ) );

  // bp::implicitly_convertible 
  //   <shared_ptr<GridFunctionCoefficientFunction>, 
  //   shared_ptr<CoefficientFunction> >(); 
  

  // void RefineAtLevelSet (PyGF gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh){

 
  // bp::docstring_options local_docstring_options(true, true, false);
  
  // std::string nested_name = "comp";
  // if( bp::scope() )
  //   nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");


  // new implementation: only straight cuts - start with triangles only for a start!

  bp::def("NewIntegrateX",
          FunctionPointer([](bp::object lset,
                             shared_ptr<MeshAccess> ma, 
                             PyCF cf,
                             int order, DOMAIN_TYPE dt, int heapsize)
                          {
                            static Timer timer ("NewIntegrateX");
                            static Timer timercutgeom ("NewIntegrateX::MakeCutGeom");
                            static Timer timerevalcoef ("NewIntegrateX::EvalCoef");
                            static Timer timeradding("NewIntegrateX::Adding");
                            static Timer timermapir("NewIntegrateX::MapingIntergrRule");

                            cout << "before problems ?" << endl;
                            bp::extract<PyGF> bpgf(lset);
                            if (!bpgf.check())
                              throw Exception("cast failed... need new candidates..");
                            shared_ptr<GridFunction> gf_lset = bpgf().Get();
                            cout << "lset is a gf of type :" << gf_lset->GetFESpace()->GetName() << endl;
                            cout << "lset is a gf of order:" << gf_lset->GetFESpace()->GetOrder() << endl;

                            RegionTimer reg (timer);
                            LocalHeap lh(heapsize, "lh-New-Integrate");

                            double sum = 0.0;
                            int DIM = ma->GetDimension();

                            ma->IterateElements
                              (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
                               {
                                 auto & trafo = ma->GetTrafo (el, lh);

                                 Array<int> dnums;
                                 gf_lset->GetFESpace()->GetDofNrs(el.Nr(),dnums);
                                 FlatVector<> elvec(dnums.Size(),lh);
                                 gf_lset->GetVector().GetIndirect(dnums,elvec);

                                 cout << "elvec: " << elvec << endl;

                                 timercutgeom.Start();
                                 const IntegrationRule * ir = StraightCutIntegrationRule(gf_lset, trafo, dt, order, lh);
                                 timercutgeom.Stop();

                                 if (ir != nullptr)
                                 {
                                   timermapir.Start();
                                   BaseMappedIntegrationRule & mir = trafo(*ir, lh);
                                   FlatMatrix<> val(mir.Size(), 1, lh);
                                   timermapir.Stop();

                                   timerevalcoef.Start();
                                   cf -> Evaluate (mir, val);
                                   timerevalcoef.Stop();

                                   timeradding.Start();
                                   double lsum = 0.0;
                                   for (int i = 0; i < mir.Size(); i++)
                                     lsum += mir[i].GetWeight()*val(i,0);

                                   AsAtomic(sum) += lsum;
                                   timeradding.Stop();
                                 }
                               });

                            return sum;
                          }),
          (bp::arg("lset"), bp::arg("mesh"), 
           bp::arg("cf")=PyCF(make_shared<ConstantCoefficientFunction>(0.0)), 
           bp::arg("order")=5, bp::arg("domain_type")=IF, bp::arg("heapsize")=1000000));

}

BOOST_PYTHON_MODULE(libngsxfem_py) 
{
  ExportNgsx();
}
#endif // NGSX_PYTHON
