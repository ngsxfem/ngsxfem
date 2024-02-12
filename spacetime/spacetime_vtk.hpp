#pragma once

/*********************************************************************/
/* File:   spacetime_vtk.hpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   28. September 2019 (based on vtkout.*pp from NGSOlve)     */
/*********************************************************************/
#include <comp.hpp>
#include "vtkoutput.hpp"

namespace ngcomp
{ 

  class SpaceTimeVTKOutput
  {
  protected:
    shared_ptr<MeshAccess> ma = nullptr;
    Array<shared_ptr<CoefficientFunction>> coefs;
    Array<string> fieldnames;
    string filename;
    int subdivisionx, subdivisiont;
    int only_element = -1;

    Array<shared_ptr<ValueField>> value_field;
    Array<Vec<3>> points;
    Array<IVec<ELEMENT_MAXPOINTS+1>> cells;

    int output_cnt = 0;
    
    shared_ptr<ofstream> fileout;
    
  public:

    SpaceTimeVTKOutput (const Array<shared_ptr<CoefficientFunction>> &,
               const Flags &,shared_ptr<MeshAccess>);

    SpaceTimeVTKOutput (shared_ptr<MeshAccess>, const Array<shared_ptr<CoefficientFunction>> &,
                        const Array<string> &, string, int, int, int);
    virtual ~SpaceTimeVTKOutput() { ; }
    
    void ResetArrays();
    
    void FillReferenceHex(Array<IntegrationPoint> & ref_coords,Array<IVec<ELEMENT_MAXPOINTS+1>> & ref_elems); 
    void FillReferencePrism(Array<IntegrationPoint> & ref_coords,Array<IVec<ELEMENT_MAXPOINTS+1>> & ref_elems);    

    void PrintPoints();
    void PrintCells();
    void PrintCellTypes(VorB vb, const BitArray * drawelems=nullptr);
    void PrintFieldData();    

    virtual void Do (LocalHeap & lh, VorB vb = VOL, const BitArray * drawelems = 0, double t_start = 0, double t_end = 1);
  };

}


