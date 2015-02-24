#ifdef NGSX_PYTHON
//#include "../ngstd/python_ngstd.hpp"
#include <python_ngstd.hpp>
#include "../xfem/xFESpace.hpp"

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
