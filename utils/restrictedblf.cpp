#include "restrictedblf.hpp"
#include <comp.hpp>

namespace ngcomp
{

  RestrictedBilinearForm :: 
  RestrictedBilinearForm (shared_ptr<FESpace> afespace,
                          const string & aname,
                          shared_ptr<BitArray> ael_restriction,
                          shared_ptr<BitArray> afac_restriction,
                          const Flags & flags)
    : T_BilinearForm<double,double>(afespace, aname, flags),
      el_restriction(ael_restriction),
      fac_restriction(afac_restriction)
  {
    ;
  }
  
  
  
  MatrixGraph * RestrictedBilinearForm :: GetGraph (int level, bool symmetric)
  {
    static Timer timer ("BilinearForm::GetGraph");
    RegionTimer reg (timer);

    int ndof = fespace->GetNDof();
    int nf = ma->GetNFacets();
    int neV = ma->GetNE(VOL);
    int neB = ma->GetNE(BND);
    int neBB = ma->GetNE(BBND);
    const Array<SpecialElement*> & specialelements = fespace->GetSpecialElements();
    int nspe = specialelements.Size();

    Array<DofId> dnums;
    Array<int> fnums; //facets of one element
    Array<int> elnums; //elements neighbouring one facet
    Array<int> nbelems; //neighbour elements


    int maxind = neV + neB + neBB + specialelements.Size();
    if (fespace->UsesDGCoupling()) maxind += nf;

    TableCreator<int> creator(maxind);
    for ( ; !creator.Done(); creator++)
      {
	for(VorB vb : {VOL, BND, BBND})
	  {
	    int nre = ma->GetNE(vb);
	    ParallelForRange (Range(nre), [&](IntRange r)
			      {
				Array<DofId> dnums;
				for (auto i : r)
				  {
                                    // if (vb == VOL)
                                    //   if (el_restriction && (! el_restriction->Test(i)))
                                    //     continue;
                                    if (vb == VOL)
                                      if (el_restriction && (! el_restriction->Test(i)))
                                        continue;
				    auto eid = ElementId(vb,i);
				    if (!fespace->DefinedOn (vb,ma->GetElIndex(eid)))
                                      continue;
				    
				    if (vb == VOL && eliminate_internal)
				      fespace->GetDofNrs (i, dnums, EXTERNAL_DOF);
				    else
				      fespace->GetDofNrs (eid, dnums);
				    int shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV+neB);
				    for (int d : dnums)
				      if (d != -1) creator.Add (shift+i, d);
                              }
			      });
	  }

        
        for (int i = 0; i < specialelements.Size(); i++)
          {
            specialelements[i]->GetDofNrs (dnums);
            for (int d : dnums)
              if (d != -1) creator.Add (neV+neB+neBB+i, d);
          }

        if (fespace->UsesDGCoupling())
        {
          //add dofs of neighbour elements as well
          Array<DofId> dnums_dg;
          Array<int> elnums_per;

          for (int i = 0; i < nf; i++)
            {
              if (fac_restriction && (! fac_restriction->Test(i)))
                continue;
              nbelems.SetSize(0);
              ma->GetFacetElements(i,elnums);
              for (int k=0; k<elnums.Size(); k++)
                nbelems.Append(elnums[k]);

              if(nbelems.Size() < 2)
              {
                int facet2 = ma->GetPeriodicFacet(i);
                if(facet2 != i)
                {
                  ma->GetFacetElements (facet2, elnums_per);
                  nbelems.Append(elnums_per[0]);
                }
              }
              dnums_dg.SetSize(0);
              for (int k=0;k<nbelems.Size();k++){
                int elnr=nbelems[k];
                if (!fespace->DefinedOn (VOL,ma->GetElIndex(ElementId(VOL,elnr)))) continue;
                fespace->GetDofNrs (ElementId(VOL,elnr), dnums);
                dnums_dg.Append(dnums);
              }
              QuickSort (dnums_dg);
              for (int j = 0; j < dnums_dg.Size(); j++)
                if (dnums_dg[j] != -1 && (j==0 || (dnums_dg[j] != dnums_dg[j-1]) ))
                  creator.Add (neV+neB+neBB+nspe+i, dnums_dg[j]);
            }
        }

      }
    
    MatrixGraph * graph;

    if (!fespace2)
      {
        auto table = creator.MoveTable();
        graph = new MatrixGraph (ndof, table, table, symmetric);
      }
    else
      {
        throw Exception("not yet implemented");
      }

    graph -> FindSameNZE();
    return graph;
  }
}
