/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
// #include "xintegration.hpp"
#include "../spacetime/spacetimefespace.hpp"
#include "../xfem/xfemIntegrators.hpp"
#include "../xfem/stxfemIntegrators.hpp"
#include "../xfem/setvaluesx.hpp"
#include "../utils/error.hpp"
#include "../utils/output.hpp"
#include "../utils/calccond.hpp"
#include "../xfem/xFESpace.hpp"

using namespace ngsolve;
// using namespace xintegration;
using namespace ngfem;

namespace ngcomp
{ 

  class ValueField : public Array<double>
  {
    string name = "none";
  public:
    void SetName(string aname){ name = aname; }
    string Name(){ return name;}
  };
  
/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  template <int D> 
  class VTKOutput
  {
  protected:

    shared_ptr<MeshAccess> ma = nullptr;
    
    Array<shared_ptr<GridFunction>> gfus;
    Array<shared_ptr<CoefficientFunction>> coefs;
    
    int subdivision;
    bool onlygrid;
    Flags myflags;
    
    shared_ptr<ofstream> fileout;

    Array<Vec<D>> points;
    Array<INT<D+1>> cells;

    int n_coef_fields;
    int n_gf_fields;
    
    Array<shared_ptr<ValueField>> value_field;

    string filename;
    
  public:

    VTKOutput (const Array<shared_ptr<CoefficientFunction>> & a_coef_lset,
               const Array<shared_ptr<GridFunction>> & a_gfus,
               const Flags & flags,
               shared_ptr<MeshAccess> ama = nullptr);
    
    void ResetArrays();
    
    void FillReferenceData2D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_trigs);    
    void FillReferenceData3D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets);
    void PrintPoints();
    void PrintCells();
    void PrintCellTypes();
    void PrintFieldData();    

    void Do (LocalHeap & lh);
  };


  class NumProcVTKOutput : public NumProc
  {
  protected:
    shared_ptr<VTKOutput<2>> vtkout2 = nullptr;
    shared_ptr<VTKOutput<3>> vtkout3 = nullptr;
  public:
    NumProcVTKOutput (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcVTKOutput() { }

    virtual string GetClassName () const
    {
      return "NumProcVTKOutput";
    }
    
    virtual void Do (LocalHeap & lh);
  };
  
}
