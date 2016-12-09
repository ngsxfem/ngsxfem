/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
// #include "xintegration.hpp"
// #include "../spacetime/spacetimefespace.hpp"
#include "../xfem/xfemIntegrators.hpp"
// #include "../xfem/stxfemIntegrators.hpp"
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

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
  // template <int D> 
  class NumProcTraceOutput : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    shared_ptr<CoefficientFunction> coef_lset;
    int subdivision;
    bool onlygrid;
    Flags myflags;
    shared_ptr<ofstream> fileout;

    Array<Vec<3>> points;
    Array<INT<4>> cells;
    Array<double> values_lset;
    Array<double> values_tracesol;
    
    bool breset;
    bool btime;

    string filename;
  public:


    NumProcTraceOutput (shared_ptr<PDE> apde, const Flags & flags)
        : NumProc (apde), myflags(flags)
    { 
      gfu  = apde->GetGridFunction (flags.GetStringFlag ("gridfunction"));
      coef_lset  = apde->GetCoefficientFunction (flags.GetStringFlag ("levelset"));
      subdivision = (int) flags.GetNumFlag ( "subdivision", 0);
      onlygrid = flags.GetDefineFlag ("onlymesh");
      filename = flags.GetStringFlag ("filename","tracefem");
      btime = flags.GetDefineFlag ("instat");
      breset = flags.GetDefineFlag ("reset");
			
    }

    void ResetArrays()
    {
      points.SetSize(0);
      cells.SetSize(0);
      values_lset.SetSize(0);
      values_tracesol.SetSize(0);
    }






    virtual ~NumProcTraceOutput()
    {
    }

    virtual string GetClassName () const
    {
      return "NumProcTraceOutput";
    }

    void FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<4>> & ref_tets)
    {
      if (subdivision == 0)
      {
        ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
        ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
        ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
        ref_coords.Append(IntegrationPoint(0.0,0.0,1.0));
        ref_tets.Append(INT<4>(0,1,2,3));
      }
      else
      {
        const int r = pow(2,subdivision);
        const int s = r + 1;

        // std::cout << " r = " << r << std::endl;
        // std::cout << " s = " << s << std::endl;

        Array<INT<3>> pidx_to_ijk ( (r+1)*(r+2)*(r+3) / 6);

        
        const double h = 1.0/r;

        int pidx = 0;
        for (int i = 0; i <= r; ++i)
          for (int j = 0; i+j <= r; ++j)
            for (int k = 0; i+j+k <= r; ++k)
            {
              ref_coords.Append(IntegrationPoint(i*h,j*h,k*h));
              pidx_to_ijk[pidx++] = INT<3>(i,j,k);
            }

        pidx = 0;
        for (int i = 0; i <= r; ++i)
          for (int j = 0; i+j <= r; ++j)
            for (int k = 0; i+j+k <= r; ++k, pidx++)
            {
              if (i+j+k == r) continue;
              int pidx_curr = pidx;
              int pidx_incr_k = pidx+1;
              int pidx_incr_j = pidx+s-i-j;
              int pidx_incr_i = pidx+(s-i)*(s+1-i)/2-j;

              int pidx_incr_kj = pidx_incr_j + 1;

              int pidx_incr_ij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j;
              int pidx_incr_ki = pidx+(s-i)*(s+1-i)/2-j + 1;
              int pidx_incr_kij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j + 1;

              ref_tets.Append(INT<4>(pidx,pidx_incr_k,pidx_incr_j,pidx_incr_i));
              if (i+j+k+1 == r)
                continue;

              // std::cout << " i = " << i << std::endl;
              // std::cout << " j = " << j << std::endl;
              // std::cout << " k = " << k << std::endl;
              // std::cout << " pidx = " << pidx << std::endl;
              // std::cout << " pidx_incr_k = " << pidx_incr_k << std::endl;
              // std::cout << " pidx_incr_j = " << pidx_incr_j << std::endl;
              // std::cout << " pidx_incr_i = " << pidx_incr_i << std::endl;
              // std::cout << " pidx_incr_ki = " << pidx_incr_ki << std::endl;
              // std::cout << " pidx_incr_kj = " << pidx_incr_kj << std::endl;
              // std::cout << " pidx_incr_ij = " << pidx_incr_ij << std::endl;
              // std::cout << " pidx_incr_kij = " << pidx_incr_kij << std::endl;

              ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_j,pidx_incr_i));
              ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_ki,pidx_incr_i));

              // ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ki,pidx_incr_kij,pidx_incr_i));

              ref_tets.Append(INT<4>(pidx_incr_j,pidx_incr_i,pidx_incr_kj,pidx_incr_ij));
              ref_tets.Append(INT<4>(pidx_incr_i,pidx_incr_kj,pidx_incr_ij,pidx_incr_ki));
              
              if (i+j+k+2 != r)
                ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ij,pidx_incr_ki,pidx_incr_kij));
            }              
      }
    }

    void PrintPoints()
    {
      *fileout << "POINTS " << points.Size() << " float" << endl;
      for (auto p : points)
        *fileout << p << endl;
    }

    void PrintCells()
    {
      *fileout << "CELLS " << cells.Size() << " " << 5 * cells.Size() << endl;
      for (auto c : cells)
        *fileout << "4 " << c << endl;
    }

    void PrintCellTypes()
    {
      *fileout << "CELL_TYPES " << cells.Size() << endl;
      for (auto c : cells)
        *fileout << "10 " << endl;
      *fileout << "CELL_DATA " << cells.Size() << endl;
      *fileout << "POINT_DATA " << points.Size() << endl;
    }

    void PrintFieldData()
    {
      *fileout << "FIELD FieldData 2"<< endl;

      *fileout << "levelset 1 " << values_lset.Size() << " float" << endl;
      for (auto v : values_lset)
        *fileout << v << " ";
      *fileout << endl;

      *fileout << "solution 1 " << values_tracesol.Size() << " float" << endl;
      for (auto v : values_tracesol)
        *fileout << v << " ";
      *fileout << endl;
    }
    

    virtual void Do (LocalHeap & lh)
    {
      static int refinements = 0;
      static int time = 0;
      if (breset)
      {
        refinements++;
        time = 0;
        cout <<refinements;
      }  
      else
      {

        ostringstream filenamefinal;
        if (!btime)
        {
          filenamefinal << filename << refinements << ".vtk";
          cout << " This is the Do-call on refinement level " << refinements << std::endl;
          refinements++;
        }
        else
        {
          filenamefinal << filename << refinements << ".vtk." << time;
          time++;
        }


        fileout = make_shared<ofstream>(filenamefinal.str());



        ResetArrays();

        Array<IntegrationPoint> ref_vertices(0);
        Array<INT<4>> ref_tets(0);

        FillReferenceData(ref_vertices,ref_tets);
        
        // header:
        *fileout << "# vtk DataFile Version 3.0" << endl;
        *fileout << "vtk output" << endl;
        *fileout << "ASCII" << endl;
        *fileout << "DATASET UNSTRUCTURED_GRID" << endl;


        XFESpace & xfes = * dynamic_pointer_cast<XFESpace>(gfu->GetFESpace());
        int ne = gfu->GetMeshAccess()->GetNE();
        for ( int i : Range(ne))
        {
          if (!xfes.IsElementCut(i)) continue;
          HeapReset hr(lh);

          ElementTransformation & eltrans = gfu->GetMeshAccess()->GetTrafo (i, 0, lh);
          const FiniteElement & fel = xfes.GetFE (i, lh);
          const XFiniteElement & xfe =
            dynamic_cast<const XFiniteElement &> (fel);
          const ScalarFiniteElement<3> & scafe =
            dynamic_cast<const ScalarFiniteElement<3> & > (xfe.GetBaseFE());

          int ndof = scafe.GetNDof();
          FlatVector<> shape(ndof,lh);

          Array<int> dnums (ndof, lh);
          xfes.GetDofNrs (i, dnums);
          FlatVector<> elvec(ndof,lh);
          gfu->GetVector().GetIndirect (dnums, elvec);

          int offset = points.Size();
          for ( auto ip : ref_vertices)
          {
            MappedIntegrationPoint<3,3> mip(ip, eltrans);
            points.Append(mip.GetPoint());
            values_lset.Append(coef_lset->Evaluate(mip));

            scafe.CalcShape (mip.IP(),shape);
            values_tracesol.Append(InnerProduct(shape,elvec));
          }

          for ( auto tet : ref_tets)
          {
            INT<4> new_tet = tet;
            for (int i = 0; i < 4; ++i)
              new_tet[i] += offset;
            cells.Append(new_tet);
          }

        }
        PrintPoints();
        PrintCells();
        PrintCellTypes();
        PrintFieldData();

      }
    }    
    

  };
}

static RegisterNumProc<NumProcTraceOutput> nptraceout("traceoutput");
