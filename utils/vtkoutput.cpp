/*********************************************************************/
/* File:   vtkoutput.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   1. June 2014                                              */
/*********************************************************************/

#include "vtkoutput.hpp"

namespace ngcomp
{ 

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */


  template <int D> 
  VTKOutput<D>::VTKOutput (const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                           const Array<shared_ptr<GridFunction>> & a_gfus,
                           const Flags & flags,
                           shared_ptr<MeshAccess> ama )
    : ma(ama), gfus(a_gfus), coefs(a_coefs), myflags(flags)
  {
    if ((a_gfus.Size()==0) && (ama == nullptr))
      throw Exception("need mesh access");
    if (!ma)
      ma = a_gfus[0]->GetMeshAccess();
    
    subdivision = (int) flags.GetNumFlag ( "subdivision", 0);
    onlygrid = flags.GetDefineFlag ("onlymesh");
    filename = flags.GetStringFlag ("filename","output");

    Array<string> fieldnames(flags.GetStringListFlag ("fieldnames" ));
    
    n_coef_fields = a_coefs.Size();
    n_gf_fields = a_gfus.Size();
    
    value_field.SetSize(n_coef_fields+n_gf_fields);
    for (int i = 0; i < n_coef_fields+n_gf_fields; i++)
    {
      value_field[i] = make_shared<ValueField>();
      if (fieldnames.Size() > i)
        value_field[i]->SetName(fieldnames[i]);
      else 
        value_field[i]->SetName("dummy" + to_string(i));
    }

  }

  template <int D> 
  void VTKOutput<D>::ResetArrays()
  {
    points.SetSize(0);
    cells.SetSize(0);
    for (auto field : value_field)
      field->SetSize(0);
  }
    
  template <int D> 
  void VTKOutput<D>::FillReferenceData2D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_trigs)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_trigs.Append(INT<D+1>(0,1,2));
    }
    else
    {
      throw Exception("not yet implemented");
    }
  }
    
  template <int D> 
  void VTKOutput<D>::FillReferenceData3D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,0.0,1.0));
      ref_tets.Append(INT<D+1>(0,1,2,3));
    }
    else
    {
      throw Exception("not yet implemented");
        
      // const int r = pow(2,subdivision);
      // const int s = r + 1;

      // // std::cout << " r = " << r << std::endl;
      // // std::cout << " s = " << s << std::endl;

      // Array<INT<D+1>> pidx_to_ijk ( (r+1)*(r+2)*(r+3) / 6);

        
      // const double h = 1.0/r;

      // int pidx = 0;
      // for (int i = 0; i <= r; ++i)
      //   for (int j = 0; i+j <= r; ++j)
      //     for (int k = 0; i+j+k <= r; ++k)
      //     {
      //       ref_coords.Append(IntegrationPoint(i*h,j*h,k*h));
      //       pidx_to_ijk[pidx++] = INT<3>(i,j,k);
      //     }

      // pidx = 0;
      // for (int i = 0; i <= r; ++i)
      //   for (int j = 0; i+j <= r; ++j)
      //     for (int k = 0; i+j+k <= r; ++k, pidx++)
      //     {
      //       if (i+j+k == r) continue;
      //       int pidx_curr = pidx;
      //       int pidx_incr_k = pidx+1;
      //       int pidx_incr_j = pidx+s-i-j;
      //       int pidx_incr_i = pidx+(s-i)*(s+1-i)/2-j;

      //       int pidx_incr_kj = pidx_incr_j + 1;

      //       int pidx_incr_ij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j;
      //       int pidx_incr_ki = pidx+(s-i)*(s+1-i)/2-j + 1;
      //       int pidx_incr_kij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j + 1;

      //       ref_tets.Append(INT<4>(pidx,pidx_incr_k,pidx_incr_j,pidx_incr_i));
      //       if (i+j+k+1 == r)
      //         continue;

      //       // std::cout << " i = " << i << std::endl;
      //       // std::cout << " j = " << j << std::endl;
      //       // std::cout << " k = " << k << std::endl;
      //       // std::cout << " pidx = " << pidx << std::endl;
      //       // std::cout << " pidx_incr_k = " << pidx_incr_k << std::endl;
      //       // std::cout << " pidx_incr_j = " << pidx_incr_j << std::endl;
      //       // std::cout << " pidx_incr_i = " << pidx_incr_i << std::endl;
      //       // std::cout << " pidx_incr_ki = " << pidx_incr_ki << std::endl;
      //       // std::cout << " pidx_incr_kj = " << pidx_incr_kj << std::endl;
      //       // std::cout << " pidx_incr_ij = " << pidx_incr_ij << std::endl;
      //       // std::cout << " pidx_incr_kij = " << pidx_incr_kij << std::endl;

      //       ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_j,pidx_incr_i));
      //       ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_ki,pidx_incr_i));

      //       // ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ki,pidx_incr_kij,pidx_incr_i));

      //       ref_tets.Append(INT<4>(pidx_incr_j,pidx_incr_i,pidx_incr_kj,pidx_incr_ij));
      //       ref_tets.Append(INT<4>(pidx_incr_i,pidx_incr_kj,pidx_incr_ij,pidx_incr_ki));
              
      //       if (i+j+k+2 != r)
      //         ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ij,pidx_incr_ki,pidx_incr_kij));
      //     }              
    }
  }

  template <int D> 
  void VTKOutput<D>::PrintPoints()
  {
    *fileout << "POINTS " << points.Size() << " float" << endl;
    for (auto p : points)
    {
      *fileout << p;
      if (D==2)
        *fileout << "\t 0.0";
      *fileout << endl;
    }
  }

  template <int D> 
  void VTKOutput<D>::PrintCells()
  {
    *fileout << "CELLS " << cells.Size() << " " << (D+2) * cells.Size() << endl;
    for (auto c : cells)
      *fileout << D+1 <<" " << c << endl;
  }

  template <int D> 
  void VTKOutput<D>::PrintCellTypes()
  {
    *fileout << "CELL_TYPES " << cells.Size() << endl;
    if (D==3)
      for (auto c : cells)
        *fileout << "10 " << endl;
    else
      for (auto c : cells)
        *fileout << "5 " << endl;
    *fileout << "CELL_DATA " << cells.Size() << endl;
    *fileout << "POINT_DATA " << points.Size() << endl;
  }

  template <int D> 
  void VTKOutput<D>::PrintFieldData()
  {
    *fileout << "FIELD FieldData " << value_field.Size() << endl;

    for (auto field : value_field)
    {
      *fileout << field->Name() << " 1 " << field->Size() << " float" << endl;
      for (auto v : *field)
        *fileout << v << " ";
      *fileout << endl;
    }
    
  }
    

  template <int D> 
  void VTKOutput<D>::Do (LocalHeap & lh)
  {
    static int refinements = 0;
    ostringstream filenamefinal;
    filenamefinal << filename << refinements << ".vtk";
    fileout = make_shared<ofstream>(filenamefinal.str());
    cout << " This is the Do-call on refinement level " << refinements << std::endl;
    refinements++;

    ResetArrays();

    Array<IntegrationPoint> ref_vertices(0);
    Array<INT<D+1>> ref_tets(0);

    if (D==3)
      FillReferenceData3D(ref_vertices,ref_tets);
    else
      FillReferenceData2D(ref_vertices,ref_tets);
      
    // header:
    *fileout << "# vtk DataFile Version 3.0" << endl;
    *fileout << "vtk output" << endl;
    *fileout << "ASCII" << endl;
    *fileout << "DATASET UNSTRUCTURED_GRID" << endl;

    int ne = ma->GetNE();
    for ( int elnr : Range(ne))
    {
      HeapReset hr(lh);

      ElementTransformation & eltrans = ma->GetTrafo (elnr, 0, lh);

      int offset = points.Size();
      for ( auto ip : ref_vertices)
      {
        MappedIntegrationPoint<D,D> mip(ip, eltrans);
        points.Append(mip.GetPoint());
      }
      
      for (int i = 0; i < n_coef_fields; i++)
      {
        for ( auto ip : ref_vertices)
        {
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          value_field[i]->Append(coefs[i]->Evaluate(mip));
        }
      }
      
      for (int i = n_coef_fields; i < n_coef_fields + n_gf_fields; i++)
      {
        int k = i - n_coef_fields;
        const FiniteElement & fel = gfus[k]->GetFESpace()->GetFE (elnr, lh);
        const ScalarFiniteElement<D> & scafe =
          dynamic_cast<const ScalarFiniteElement<D> & > (fel);
        int ndof = scafe.GetNDof();
        FlatVector<> shape(ndof,lh);
        Array<int> dnums (ndof, lh);
        gfus[k]->GetFESpace()->GetDofNrs (elnr, dnums);
        FlatVector<> elvec(ndof,lh);
        gfus[k]->GetVector().GetIndirect (dnums, elvec);
        
        for ( auto ip : ref_vertices)
        {
          scafe.CalcShape (ip,shape);
          value_field[i]->Append(InnerProduct(shape,elvec));
        }
      }

      for ( auto tet : ref_tets)
      {
        INT<D+1> new_tet = tet;
        for (int i = 0; i < D+1; ++i)
          new_tet[i] += offset;
        cells.Append(new_tet);
      }

    }

    PrintPoints();
    PrintCells();
    PrintCellTypes();
    PrintFieldData();
      
  }    

  NumProcVTKOutput::NumProcVTKOutput (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    const Array<string> & coefs_strings = flags.GetStringListFlag ("coefficients");
    const Array<string> & gfs_strings = flags.GetStringListFlag ("gridfunctions");
    
    Array<shared_ptr<CoefficientFunction>> coefs;
    for (int i = 0; i < coefs_strings.Size(); ++i)
      coefs.Append(apde->GetCoefficientFunction (coefs_strings[i]));

    Array<shared_ptr<GridFunction>> gfus;
    for (int i = 0; i < gfs_strings.Size(); ++i)
    {
      gfus.Append(apde->GetGridFunction (gfs_strings[i]));
    }
    
    if (apde->GetMeshAccess()->GetDimension() == 2)
      vtkout2 = make_shared<VTKOutput<2>>(coefs, gfus, flags, apde->GetMeshAccess());
    else 
      vtkout3 = make_shared<VTKOutput<3>>(coefs, gfus, flags, apde->GetMeshAccess());
  }


  void NumProcVTKOutput::Do (LocalHeap & lh)
  {
    if (vtkout2)
      vtkout2->Do(lh);
    if (vtkout3)
      vtkout3->Do(lh);
  }
  
}

static RegisterNumProc<NumProcVTKOutput> npvtkout("vtkoutput");
