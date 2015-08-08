#ifdef NGSX_PYTHON
//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../xfem/xFESpace.hpp"
#include "../stokes/xstokesspace.hpp"
#include "../lsetcurving/p1interpol.hpp"
#include "../lsetcurving/calcgeomerrors.hpp"
#include "../lsetcurving/lsetrefine.hpp"
#include "../lsetcurving/projshift.hpp"

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
    ;

  bp::def("CastToXStdFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<XStdFESpace>(fes); } ) );
  bp::def("CastToXFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<XFESpace>(fes); } ) );
  bp::def("XToNegPos", FunctionPointer( [] (shared_ptr<GridFunction> gfx, shared_ptr<GridFunction> gfnegpos) { XFESpace::XToNegPos(gfx,gfnegpos); } ) );

  // bp::def("CastToFESpace", FunctionPointer( [] (shared_ptr<FESpace> fes) { return dynamic_pointer_cast<FESpace>(fes); } ) );


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
