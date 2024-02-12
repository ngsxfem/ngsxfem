#include "restrictedblf.hpp"
#include <comp.hpp>

namespace ngcomp
{
  
  //associate elements / specialelements / facets (DG) to dofs for precomputed sparsity pattern
  Table<int> MeshEntityToDofTable(shared_ptr<FESpace> fes, 
                                  shared_ptr<BitArray> active_elems = nullptr,
                                  shared_ptr<BitArray> active_facets = nullptr,
                                  bool eliminate_internal = false, 
                                  bool eliminate_hidden = false, 
                                  const Array<unique_ptr<SpecialElement>> * specels = nullptr)
  {
    if (eliminate_internal)
      if (fes->UsesDGCoupling())
        cout << IM(4) << "static condensation should work here, but the user should make sure that all local dofs are not involved in dg-couplings." << endl;
    if (eliminate_hidden)
      if (fes->UsesDGCoupling())
        cout << IM(4) << "eliminating hidden should work here, but the user should make sure that the hidden dofs are not involved in dg-couplings." << endl;

    shared_ptr<MeshAccess> ma = fes->GetMeshAccess();
    
    size_t ndof = fes->GetNDof();
    size_t nf = ma->GetNFacets();
    size_t neV = ma->GetNE(VOL);
    size_t neB = ma->GetNE(BND);
    size_t neBB = ma->GetNE(BBND);
    size_t nspe = 0;
    if (specels)
      nspe = specels->Size();
    
    int maxind = neV + neB + neBB + nspe;
    if (fes->UsesDGCoupling()) 
      maxind += nf;

    TableCreator<int> creator(maxind);

    for ( ; !creator.Done(); creator++)
    {
      for(VorB vb : {VOL, BND, BBND})
      {
        size_t shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV+neB);
        bool condensation_allowed = (vb == VOL) || ((neV==0) && (vb == BND));
	      ParallelForRange (ma->GetNE(vb), [&](IntRange r)
	      {
          Array<DofId> dnums;
          for (auto i : r)
          {
            if (vb == VOL)
              if (active_elems && (! active_elems->Test(i)))
                continue;
            auto eid = ElementId(vb,i);
            if (!fes->DefinedOn (vb, ma->GetElIndex(eid))) continue;
            if (condensation_allowed && eliminate_internal)
              fes->GetDofNrs (eid, dnums, EXTERNAL_DOF);
            else if (condensation_allowed && eliminate_hidden)
              fes->GetDofNrs (eid, dnums, VISIBLE_DOF);
            else
              fes->GetDofNrs (eid, dnums);
            for (DofId d : dnums)
              if (IsRegularDof(d)) creator.Add (shift+i, d);
          }
        });
      }
      if (specels) 
      {
        Array<DofId> dnums;
        for (int i = 0; i < specels->Size(); i++)
        {
          (*specels)[i]->GetDofNrs (dnums);
          QuickSort(dnums);
          int last = -1;
          for (int d : dnums)
          {
            if (d!=last && IsRegularDof(d))
              creator.Add (neV+neB+neBB+i, d);
            last = d;
          }
        }
      }

      if (fes->UsesDGCoupling())
      {
	      ParallelForRange (nf, [&](IntRange r)
	      {
          //add dofs of neighbour elements as well
          //by associating those to facets
          Array<DofId> dnums;
          Array<DofId> dnums_dg;
          Array<int> elnums;
          Array<int> elnums_per;
          Array<int> nbelems;
          for (auto i : r)
          {
            if (active_facets && (! active_facets->Test(i)))
              continue;

            nbelems.SetSize(0);
            ma->GetFacetElements(i,elnums);
  
            // timerDG1.Stop();
            if(elnums.Size() < 2)
            {
              int facet2 = ma->GetPeriodicFacet(i);
              if(facet2 > i)
              {
                ma->GetFacetElements (facet2, elnums_per);
                // if the facet is identified across subdomain
                // boundary, we only have the surface element
                // and not the other volume element!
                if (elnums_per.Size())
                  elnums.Append(elnums_per[0]);
              }
            }
                
            for (int k=0; k<elnums.Size(); k++)
              nbelems.Append(elnums[k]);
  
            dnums_dg.SetSize(0);
            for (int k=0;k<nbelems.Size();k++)
            {
              int elnr=nbelems[k];
              if (!fes->DefinedOn (VOL,ma->GetElIndex(ElementId(VOL,elnr)))) continue;
              fes->GetDofNrs (ElementId(VOL,elnr), dnums);
              dnums_dg.Append(dnums);
            }
            QuickSort (dnums_dg);
            for (int j = 0; j < dnums_dg.Size(); j++)
              if (IsRegularDof(dnums_dg[j]) && (j==0 || (dnums_dg[j] != dnums_dg[j-1]) ))
                creator.Add (neV+neB+neBB+nspe+i, dnums_dg[j]);
          }
        });
      }
    }
    return creator.MoveTable();
  }

  
   
  template <class TM, class TV> 
  MatrixGraph RestrictedBilinearForm<TM,TV> :: GetGraph (int level, bool symmetric)
  {
    static Timer timer ("BilinearForm::GetGraph");
    RegionTimer reg (timer);

    size_t ndof = this->fespace->GetNDof();

    auto table = MeshEntityToDofTable(this->fespace, el_restriction, fac_restriction, this->eliminate_internal, this->eliminate_hidden, &(this->specialelements) );
    MatrixGraph * graph;
  
    if (!(this->fespace2))
      graph = new MatrixGraph (ndof, ndof, table, table, symmetric);        
    else
    {
      auto table2 = MeshEntityToDofTable(this->fespace2, el_restriction, fac_restriction, this->eliminate_internal, this->eliminate_hidden, &(this->specialelements));
      size_t ndof2 = this->fespace2->GetNDof();
      graph = new MatrixGraph (ndof2, ndof, table2, table, symmetric);
    }
    
    graph -> FindSameNZE();
    return std::move(*graph);
  }

template class RestrictedBilinearForm<double,double>;
template class RestrictedBilinearForm<Complex,Complex>; 
} 

