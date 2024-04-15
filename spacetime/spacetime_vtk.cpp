/*********************************************************************/
/* File:   spacetime_vtk.hpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   28. September 2019 (based on vtkout.*pp from NGSOlve)     */
/*********************************************************************/

#include "spacetime_vtk.hpp"
#include "../utils/ngsxstd.hpp"

namespace ngcomp
{ 

  SpaceTimeVTKOutput::SpaceTimeVTKOutput (shared_ptr<MeshAccess> ama,
                                          const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                                          const Array<string> & a_field_names,
                                          string a_filename,
                                          int a_subdivision_x, int a_subdivision_t, int a_only_element)
    : ma(ama), coefs(a_coefs), fieldnames(a_field_names),
      filename(a_filename),
      subdivisionx(a_subdivision_x), subdivisiont(a_subdivision_t), only_element(a_only_element)
  {
    value_field.SetSize(a_coefs.Size());
    for (int i = 0; i < a_coefs.Size(); i++)
      if (fieldnames.Size() > i)
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),fieldnames[i]);
      else 
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),"dummy" + to_string(i));
  }


  /// Empty all field 
  void SpaceTimeVTKOutput::ResetArrays()
  {
    points.SetSize(0);
    cells.SetSize(0);
    for (auto field : value_field)
      field->SetSize(0);
  }
    

  /// Fill principil lattices (points and connections on subdivided reference hexahedron) in 3D
  void SpaceTimeVTKOutput::FillReferenceHex(Array<IntegrationPoint> & ref_coords,Array<IVec<ELEMENT_MAXPOINTS+1>> & ref_elems)
  {
    if(subdivisionx == 0 && subdivisiont == 0)
    {
      double p[8][3] = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };
      for (int i = 0; i < 8; i++)
      {
        auto tmp = IntegrationPoint(p[i][0],p[i][1],0.0,p[i][2]);
        MarkAsSpaceTimeIntegrationPoint(tmp);
        ref_coords.Append(tmp);
      }
      IVec<ELEMENT_MAXPOINTS+1> elem;
      elem[0] = 8;
      for(int i=0; i<ElementTopology::GetNVertices(ET_HEX); i++)
        elem[i+1] = i;
      ref_elems.Append(elem);
      
    }
    else
    {
      const int rx = 1<<subdivisionx;
      const int rt = 1<<subdivisiont;
      // const int s = r + 1;

      const double hx = 1.0/rx;
      const double ht = 1.0/rt;

      int pidx = 0;
      for(int i = 0; i <= rx; ++i)
        for(int j = 0; j <= rx; ++j)
          for(int k = 0; k <= rt; ++k)
          {
            auto tmp = IntegrationPoint(i*hx,j*hx,0.0,k*ht);
            MarkAsSpaceTimeIntegrationPoint(tmp);
            ref_coords.Append(tmp);
          }

      for(int i = 0; i < rx; ++i)
      {
        int incr_i = (rx+1)*(rt+1);
        for(int j = 0; j < rx; ++j)
        {
          int incr_j = rt+1;
          pidx = i*incr_i + j*incr_j;
          for(int k = 0; k < rt; ++k, pidx++)
          {
            ref_elems.Append(IVec<ELEMENT_MAXPOINTS+1>(8, pidx, pidx+1, pidx+incr_j+1, pidx+incr_j,
                                                         pidx+incr_i, pidx+incr_i+1, pidx+incr_i+incr_j+1, pidx+incr_j+incr_i));
          }          
        } 
      }
    }
  }

  
  void SpaceTimeVTKOutput::FillReferencePrism(Array<IntegrationPoint> & ref_coords,Array<IVec<ELEMENT_MAXPOINTS+1>> & ref_elems)
  {
    if(subdivisionx == 0 && subdivisiont == 0)
    {

      double p[6][3] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1} };
      for (int i = 0; i < 6; i++)
      {
        auto tmp = IntegrationPoint(p[i][0],p[i][1],0.0,p[i][2]);
        MarkAsSpaceTimeIntegrationPoint(tmp);
        ref_coords.Append(tmp);
      }
      IVec<ELEMENT_MAXPOINTS+1> elem;
      elem[0] = 6;
      for(int i=0; i<ElementTopology::GetNVertices(ET_PRISM); i++)
        elem[i+1] = i;
      ref_elems.Append(elem);
    }
    else
    {
      const int rx = 1<<subdivisionx;
      const int rt = 1<<subdivisiont;
      const int s = rx + 1;

      const double hx = 1.0/rx;
      const double ht = 1.0/rt;

      int pidx = 0;
      for(int k=0; k<=rt; k++)
        for(int i = 0; i <= rx; ++i)
          for(int j = 0; i+j <= rx; ++j)
          {
            auto tmp = IntegrationPoint(j*hx,i*hx,0.0,k*ht);
            MarkAsSpaceTimeIntegrationPoint(tmp);
            ref_coords.Append(tmp);
          }

      pidx = 0;
      for(int k=0; k<rt; k++)
      {
        int incr_k = (rx+2)*(rx+1)/2;
        pidx = k*incr_k;
        for(int i = 0; i <= rx; ++i)
          for(int j = 0; i+j <= rx; ++j,pidx++)
          {
            // int pidx_curr = pidx;
            if(i+j == rx) continue;
            int pidx_incr_i = pidx+1;
            int pidx_incr_j = pidx+s-i;

            ref_elems.Append(IVec<ELEMENT_MAXPOINTS+1>(6, pidx, pidx_incr_i, pidx_incr_j, pidx+incr_k, pidx_incr_i+incr_k, pidx_incr_j+incr_k, 0, 0 ));
              
            int pidx_incr_ij = pidx_incr_j + 1;

            if(i+j+1<rx)
              ref_elems.Append(IVec<ELEMENT_MAXPOINTS+1>(6, pidx_incr_i, pidx_incr_ij, pidx_incr_j, pidx_incr_i+incr_k, pidx_incr_ij+incr_k, pidx_incr_j+incr_k,0,0));
          }
      }
    }
  }

  /// output of data points
  void SpaceTimeVTKOutput::PrintPoints()
  {
    *fileout << "POINTS " << points.Size() << " float" << endl;
    for (auto p : points)
    {
      *fileout << p;
      *fileout << endl;
    }
  }

  /// output of cells in form vertices
  void SpaceTimeVTKOutput::PrintCells()
  {
    // count number of data for cells, one + number of vertices
    int ndata = 0;
    for (auto c : cells)
    {
      ndata++;
      ndata += c[0];
    }
    *fileout << "CELLS " << cells.Size() << " " << ndata << endl;
    for (auto c : cells)
    {
      int nv = c[0];
      *fileout << nv << "\t";
      for (int i=0; i<nv; i++)
        *fileout << c[i+1] << "\t";
      *fileout << endl;
    }
  }

  /// output of cell types (here only simplices)
  void SpaceTimeVTKOutput::PrintCellTypes(VorB vb, const BitArray * drawelems)
  {
    *fileout << "CELL_TYPES " << cells.Size() << endl;
    int factor = (1<<subdivisionx)*(1<<subdivisionx);
    // if (D==3 && vb == VOL) // 
    factor *= (1<<subdivisiont);
    for (auto e : ma->Elements(vb))
    {
      if (drawelems && !(drawelems->Test(e.Nr())))
          continue;

      switch(ma->GetElType(e))
      {
      case ET_QUAD:
        for(int i=0; i<factor; i++)
          *fileout << "12 " << endl;
        break;
      case ET_TRIG:
        for(int i=0; i<factor; i++)
          *fileout << "13 " << endl;
        break;
      default:
        cout << IM(1) << "SpaceTimeVTKOutput Element Type " << ma->GetElType(e) << " not supported!" << endl;
        throw Exception("shouldn't get this element type");
      }
    }
    *fileout << "CELL_DATA " << cells.Size() << endl;
    *fileout << "POINT_DATA " << points.Size() << endl;
  }

  /// output of field data (coefficient values)
  void SpaceTimeVTKOutput::PrintFieldData()
  {
    // *fileout << "FIELD FieldData " << value_field.Size() << endl;

    for (auto field : value_field)
    {
      *fileout << "SCALARS " << field->Name()
               << " float " << field->Dimension() << endl
               << "LOOKUP_TABLE default" << endl;
      for (auto v : *field)
        *fileout << v << " ";
      *fileout << endl;
    }
    
  }
    

  void SpaceTimeVTKOutput::Do (LocalHeap & lh, VorB vb, const BitArray * drawelems,
                               double t_start, double t_end)
  {
    ostringstream filenamefinal;
    filenamefinal << filename;
    // if (output_cnt > 0)
    filenamefinal << "_" << output_cnt;
    filenamefinal << ".vtk";
    fileout = make_shared<ofstream>(filenamefinal.str());
    cout << IM(4) << " Writing SpaceTimeVTK-Output";
    // if (output_cnt > 0)
    cout << IM(4) << " ( " << output_cnt << " )";
    cout << IM(4) << ":" << flush;
    
    output_cnt++;

    ResetArrays();

    Array<IntegrationPoint> ref_vertices_prism(0), ref_vertices_hex(0);
    Array<IVec<ELEMENT_MAXPOINTS+1>> ref_prisms(0), ref_hexes(0);
    FlatArray<IntegrationPoint> ref_vertices;
    FlatArray<IVec<ELEMENT_MAXPOINTS+1>> ref_elems;
    FillReferencePrism(ref_vertices_prism,ref_prisms);
    FillReferenceHex(ref_vertices_hex,ref_hexes);
      
    // header:
    *fileout << "# vtk DataFile Version 3.0" << endl;
    *fileout << "vtk output" << endl;
    *fileout << "ASCII" << endl;
    *fileout << "DATASET UNSTRUCTURED_GRID" << endl;

    int ne = ma->GetNE(vb);

    IntRange range = only_element >= 0 ? IntRange(only_element,only_element+1) : IntRange(ne);
    
    for ( int elnr : range)
    {
      if (drawelems && !(drawelems->Test(elnr)))
          continue;
      
      HeapReset hr(lh);

      ElementId ei(vb, elnr);
      ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
      ELEMENT_TYPE eltype = ma->GetElType(ei);

      switch(eltype)
      {
      case ET_TRIG:
        ref_vertices.Assign(ref_vertices_prism);
        ref_elems.Assign(ref_prisms);
        break;
      case ET_QUAD:
        ref_vertices.Assign(ref_vertices_hex);
        ref_elems.Assign(ref_hexes);
        break;
      default:
        throw Exception("VTK output for element-type"+ToString(eltype)+"not supported");
      }

      int offset = points.Size();
      for ( auto ip : ref_vertices)
      {
        if (vb==VOL)
        {
          MappedIntegrationPoint<2,2> mip(ip, eltrans);
          Vec<2> space_x = mip.GetPoint();
          points.Append(Vec<3>(space_x[0],space_x[1],t_start + ip.Weight() * (t_end - t_start)));
        }
        else
        {
          throw Exception("not VOL-VTK-export not yet supported");
          MappedIntegrationPoint<1,2> mip(ip, eltrans);
          Vec<2> space_x = mip.GetPoint();
          points.Append(Vec<3>(space_x[0],space_x[1],t_start + ip.Weight() * (t_end - t_start)));
        }
      }
      
      for (int i = 0; i < coefs.Size(); i++)
      {
        for ( auto ip : ref_vertices)
        {
          if (vb==VOL)
          {
            MappedIntegrationPoint<2,2> mip(ip, eltrans);

            
            const int dim = coefs[i]->Dimension();
            FlatVector<> tmp(dim,lh);
            coefs[i]->Evaluate(mip,tmp);
            for (int d = 0; d < dim; ++d)
              value_field[i]->Append(tmp(d));
          }
          else
          {
            MappedIntegrationPoint<1,2> mip(ip, eltrans);
            const int dim = coefs[i]->Dimension();
            FlatVector<> tmp(dim,lh);
            coefs[i]->Evaluate(mip,tmp);
            for (int d = 0; d < dim; ++d)
              value_field[i]->Append(tmp(d));
          }
        }
      }
      
      for ( auto elem : ref_elems)
      {
        IVec<ELEMENT_MAXPOINTS+1> new_elem = elem;
        for (int i = 1; i <= new_elem[0]; ++i)
          new_elem[i] += offset;
        cells.Append(new_elem);
      }

    }

    PrintPoints();
    PrintCells();
    PrintCellTypes(vb,drawelems);
    PrintFieldData();
      
    cout << IM(4) << " Done." << endl;
  }    

  
}

