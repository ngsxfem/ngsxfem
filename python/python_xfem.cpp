#ifdef NGSX_PYTHON
//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../xfem/xFESpace.hpp"
#include "../stokes/xstokesspace.hpp"
#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"
#include "../utils/error.hpp"

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

  bp::enum_<DOMAIN_TYPE>("DOMAIN_TYPE")
    .value("POS", POS)
    .value("NEG", NEG)
    .value("IF", IF)
    .export_values()
    ;

  bp::def("CastToXStokesFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<XStokesFESpace>(fes); } ) );

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XStokesFESpace>);
  
  bp::def("GFCoeff", FunctionPointer( [] (shared_ptr<GridFunction> in) { return dynamic_pointer_cast<CoefficientFunction>(make_shared<GridFunctionCoefficientFunction>(in)); } ) );

  // bp::implicitly_convertible 
  //   <shared_ptr<GridFunctionCoefficientFunction>, 
  //   shared_ptr<CoefficientFunction> >(); 

  bp::class_<XStokesFESpace, shared_ptr<XStokesFESpace>, bp::bases<FESpace>, boost::noncopyable>
    ("XStokesFESpace", bp::no_init)
    .def("SetLevelSet", FunctionPointer ([](XStokesFESpace & self, shared_ptr<CoefficientFunction> cf) 
                                         { self.SetLevelSet(cf); }),
         "Update information on level set function")
    .def("SetLevelSet", FunctionPointer ([](XStokesFESpace & self, shared_ptr<GridFunction> gf) 
                                         { self.SetLevelSet(gf); }),
         "Update information on level set function")
    ;


  bp::class_<XFESpace, shared_ptr<XFESpace>, bp::bases<FESpace>, boost::noncopyable>
    ("XFESpace", bp::no_init)
    .def("SetLevelSet", FunctionPointer ([](XFESpace & self, shared_ptr<CoefficientFunction> cf) 
                                         { self.SetLevelSet(cf); }),
         "Update information on level set function")
    .def("SetLevelSet", FunctionPointer ([](XFESpace & self, shared_ptr<GridFunction> gf) 
                                         { self.SetLevelSet(gf); }),
         "Update information on level set function")
    .def("SetBaseFESpace", FunctionPointer ([](XFESpace & self, shared_ptr<FESpace> fes) 
                                         { self.SetBaseFESpace(fes); }),
         "Update information on base FESpace")
    .def("BaseDofOfXDof", FunctionPointer ([](XFESpace & self, int i) 
                                         { return self.GetBaseDofOfXDof(i); }),
         "get corresponding dof of base FESpace")
    ;
  
  bp::def("CastToXStdFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<XStdFESpace>(fes); } ) );
  bp::def("CastToXFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<XFESpace>(fes); } ) );
  bp::def("XToNegPos", FunctionPointer( [] (shared_ptr<GridFunction> gfx, shared_ptr<GridFunction> gfnegpos) { XFESpace::XToNegPos(gfx,gfnegpos); } ) );

  // bp::def("CastToFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<FESpace>(fes); } ) );

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XStdFESpace>);
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<XFESpace>);

  bp::class_<XStdFESpace, shared_ptr<XStdFESpace>, bp::bases<CompoundFESpace>, boost::noncopyable>
    ("XStdFESpace", bp::no_init)
    .add_property("XFESpace", FunctionPointer ([](XStdFESpace & self) 
                                               { return dynamic_pointer_cast<XFESpace> (self[1]); }
                    ),
         // bp::return_value_policy<bp::reference_existing_object>(),
         "return XFESpace part of XStdFESpace")
    .add_property("StdFESpace", FunctionPointer ([](XStdFESpace & self) 
                                               { return self[0]; }
                    ),
         // bp::return_value_policy<bp::reference_existing_object>(),
         "return 'standard' FESpace part of XStdFESpace")
    ;

  bp::implicitly_convertible 
    <shared_ptr<XFESpace>, shared_ptr<FESpace> >(); 
  bp::implicitly_convertible 
    <shared_ptr<XStokesFESpace>, shared_ptr<FESpace> >(); 
  bp::implicitly_convertible 
    <shared_ptr<XStdFESpace>, shared_ptr<FESpace> >(); 
  
  bp::def("InterpolateToP1", FunctionPointer( [] (shared_ptr<GridFunction> gf_ho, shared_ptr<GridFunction> gf_p1, int heapsize)
                                              {
                                                InterpolateP1 interpol(gf_ho, gf_p1);
                                                LocalHeap lh (heapsize, "InterpolateP1-Heap");
                                                interpol.Do(lh);
                                              } ),
          (bp::arg("gf_ho")=NULL,bp::arg("gf_p1")=NULL,bp::arg("heapsize")=1000000))
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

  bp::def("CalcMaxDistance", FunctionPointer( [] (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> lset_p1, shared_ptr<GridFunction> deform, int heapsize)
                                              {
                                                StatisticContainer dummy;
                                                LocalHeap lh (heapsize, "CalcDistance-Heap");
                                                if (lset_p1->GetMeshAccess()->GetDimension()==2)
                                                  CalcDistances<2>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
                                                else
                                                  CalcDistances<3>(lset_ho, lset_p1, deform,  dummy, lh, -1.0, false);
                                                return (double) dummy.ErrorMaxNorm[dummy.ErrorMaxNorm.Size()-1];
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("heapsize")=1000000))
    ;

  bp::def("CalcDistances", FunctionPointer( [] (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> lset_p1, shared_ptr<GridFunction> deform, StatisticContainer & stats, int heapsize, double refine_threshold, bool absolute)
                                              {
                                                LocalHeap lh (heapsize, "CalcDistance-Heap");
                                                if (lset_p1->GetMeshAccess()->GetDimension()==2)
                                                  CalcDistances<2>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
                                                else
                                                  CalcDistances<3>(lset_ho, lset_p1, deform,  stats, lh, refine_threshold, absolute);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("stats")=NULL,bp::arg("heapsize")=1000000,bp::arg("refine_threshold")=-1.0,bp::arg("absolute")=false))
    ;

  bp::def("CalcDeformationError", FunctionPointer( [] (shared_ptr<CoefficientFunction> lset_ho, shared_ptr<GridFunction> lset_p1, shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn, StatisticContainer & stats, double lower, double upper, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "CalcDeformationError-Heap");
                                                if (lset_p1->GetMeshAccess()->GetDimension()==2)
                                                  CalcDeformationError<2>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
                                                else
                                                  CalcDeformationError<3>(lset_ho, lset_p1, deform, qn, stats, lh, lower, upper);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("qn")=NULL,bp::arg("stats")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("heapsize")=1000000))
    ;

  bp::def("ProjectShift", FunctionPointer( [] (shared_ptr<GridFunction> lset_ho, shared_ptr<GridFunction> lset_p1, shared_ptr<GridFunction> deform, shared_ptr<CoefficientFunction> qn, double lower, double upper, double threshold, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "ProjectShift-Heap");
                                                ProjectShift(lset_ho, lset_p1, deform, qn, lower, upper, threshold, lh);
                                              } ),
          (bp::arg("lset_ho")=NULL,bp::arg("lset_p1")=NULL,bp::arg("deform")=NULL,bp::arg("qn")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("threshold")=1.0,bp::arg("heapsize")=1000000))
    ;

// ProjectShift

  
  bp::def("RefineAtLevelSet", FunctionPointer( [] (shared_ptr<GridFunction> lset_p1, double lower, double upper, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
                                                RefineAtLevelSet(lset_p1, lower, upper, lh);
                                              } ),
          (bp::arg("lset_p1")=NULL,bp::arg("lower")=0.0,bp::arg("upper")=0.0,bp::arg("heapsize")=1000000))
    ;


  bp::def("CalcTraceDiff", FunctionPointer( [] (shared_ptr<GridFunction> gf, shared_ptr<CoefficientFunction> coef, int intorder, int heapsize)
                                              {
                                                Array<double> errors;
                                                LocalHeap lh (heapsize, "CalcTraceDiff-Heap");
                                                if (gf->GetMeshAccess()->GetDimension() == 2)
                                                  CalcTraceDiff<2>(gf,coef,intorder,errors,lh);
                                                else 
                                                  CalcTraceDiff<3>(gf,coef,intorder,errors,lh);
                                                return errors;
                                              } ),
          (bp::arg("gf"),bp::arg("coef"),bp::arg("intorder")=6,bp::arg("heapsize")=1000000))
    ;


  bp::def("RefineAtLevelSet", FunctionPointer( [] (shared_ptr<GridFunction> gf, double lower_lset_bound, double upper_lset_bound, int heapsize)
                                              {
                                                LocalHeap lh (heapsize, "RefineAtLevelSet-Heap");
                                                RefineAtLevelSet(gf,lower_lset_bound,upper_lset_bound,lh);
                                              } ),
          (bp::arg("gf"),bp::arg("lower_lset_bound")=0.0,bp::arg("upper_lset_bound")=0.0,bp::arg("heapsize")=10000000))
    ;


  bp::def("IntegrateX", 
          FunctionPointer([](shared_ptr<CoefficientFunction> lset,
                             shared_ptr<MeshAccess> ma, 
                             shared_ptr<CoefficientFunction> cf_neg,
                             shared_ptr<CoefficientFunction> cf_pos,
                             shared_ptr<CoefficientFunction> cf_interface,
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
                              if (cf_neg == nullptr)
                                throw Exception("no coef for neg domain given");
                              else
                                cf[int(NEG)] = cf_neg;
                            }
                            if (flags.GetDefineFlag("posdomain")){
                              tointon[int(POS)] = true;
                              if (cf_pos == nullptr)
                                throw Exception("no coef for pos domain given");
                              else
                                cf[int(POS)] = cf_pos;
                            }
                            if (flags.GetDefineFlag("interface")){
                              tointon[int(IF)] = true;
                              if (cf_interface == nullptr)
                                throw Exception("no coef for interface domain given");
                              else
                                cf[int(IF)] = cf_interface;
                            }

                            Vector<> domain_sum(3); // [val_neg,val_pos,val_interface]
                            domain_sum = 0.0;
                            int DIM = ma->GetDimension();
                            ma->IterateElements
                              (VOL, lh, [&] (Ngs_Element el, LocalHeap & lh)
                               {
                                 auto & trafo = ma->GetTrafo (el, lh);
                                 auto lset_eval
                                   = ScalarFieldEvaluator::Create(DIM,*lset,trafo,lh);
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
                                     if( ! tointon[domtype] ) continue;
                                     double hsum = 0.0;
                                     double value = 0.0;

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
                                 else if( tointon[element_domain] )
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
           bp::arg("cf_neg")=make_shared<ConstantCoefficientFunction>(0.0), 
           bp::arg("cf_pos")=make_shared<ConstantCoefficientFunction>(0.0),
           bp::arg("cf_interface")=make_shared<ConstantCoefficientFunction>(0.0),
           bp::arg("order")=5, bp::arg("subdivlvl")=0, bp::arg("domains")=bp::dict(), bp::arg("heapsize")=1000000))
    ;
  
  
  
  

  // void RefineAtLevelSet (shared_ptr<GridFunction> gf_lset_p1, double lower_lset_bound, double upper_lset_bound, LocalHeap & lh){

 
  // bp::docstring_options local_docstring_options(true, true, false);
  
  // std::string nested_name = "comp";
  // if( bp::scope() )
  //   nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".comp");
}

BOOST_PYTHON_MODULE(libngsxfem_py) 
{
  ExportNgsx();
}
#endif // NGSX_PYTHON
