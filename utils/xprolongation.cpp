#include "../utils/xprolongation.hpp"

using namespace ngcomp;

namespace ngmg
{

  void P1Prolongation :: Update (const FESpace & afes)
  {
    fes = &afes;
    if (ma->GetNLevels() > nvlevel.Size())
    {
      nvlevel.Append (ma->GetNV());
      //ndoflevel.Append (fes->GetNDof());
    }
    else
      return;

    int nv = ma->GetNV();
    shared_ptr<Array<int>> vert2dof = make_shared<Array<int>>(nv);

    tmp_vecs.Append(make_shared<VVector<double>>(fes->GetNDof()));

    Array<int> dnums(1);
    for (auto i : Range(nv))
    {
      fes->GetDofNrs(NodeId(NT_VERTEX,i),dnums);    
      if (dnums.Size() > 0 && IsRegularDof(dnums[0]))
        (*vert2dof)[i] = dnums[0];
      else
        (*vert2dof)[i] = NO_DOF_NR;
    }
    v2d_on_lvl.Append(vert2dof);
    //cout << "vert2dof : " << *vert2dof << endl;
  }


  void P1Prolongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    if (fes == nullptr)
      throw Exception("call Update before prolongating");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");    
    static Timer t("Prolongate"); RegionTimer r(t);
    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);

    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel-1]->FV<double>();
    fw = fv;
    fv = 0.0;
    for (size_t i = 0; i < nc; i++)
      if (IsRegularDof(vert2dof_fine[i]))
        fv(vert2dof_fine[i]) = fw(vert2dof_coarse[i]);

    for (size_t i = nc; i < nf; i++)
      {
        if (!IsRegularDof(vert2dof_fine[i]))
          continue;

        auto parents = ma->GetParentNodes (i);
        for (auto j : Range(2))
          //if (IsRegularDof(vert2dof_coarse[parents[j]])) //not necessary
          fv(vert2dof_fine[i]) += 0.5 * fw(vert2dof_coarse[parents[j]]);
      }

  }


  void P1Prolongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    if (fes == nullptr)
      throw Exception("call Update before restricting");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");
    static Timer t("Restrict"); RegionTimer r(t);

    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);

    size_t nc = nvlevel[finelevel-1];
    size_t nf = nvlevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel]->FV<double>();
    fw = fv;
    fv = 0;
    for (size_t i = 0; i < nc; i++)
    {
      if (!IsRegularDof(vert2dof_coarse[i]))
        continue;
      if (IsRegularDof(vert2dof_fine[i]))
        fv(vert2dof_coarse[i]) = fw(vert2dof_fine[i]);
    }

    for (size_t i = nf; i-- > nc; )
    {
      auto parents = ma->GetParentNodes (i);
      if (!IsRegularDof(vert2dof_fine[i]))
        continue;
      for (auto j : Range(2))
        //if (IsRegularDof(vert2dof_coarse[parents[j]])) // not necessary
        fv(vert2dof_coarse[parents[j]]) += 0.5 * fw(vert2dof_fine[i]);
    } 
  }

  void P2Prolongation :: Update (const FESpace & afes)
  {
    fes = &afes;
    if (ma->GetNLevels() > nVertLevel.Size() )
    {
      nVertLevel.Append (ma->GetNV());
      nEdgeLevel.Append( ma->GetNEdges() );
      
      //ndoflevel.Append (fes->GetNDof());
    }
    else
      return;

    tmp_vecs.Append(make_shared<VVector<double>>(fes->GetNDof()));

    /*
    int nv = ma->GetNV();
    // shared_ptr<Array<int>> vert2dof = make_shared<Array<int>>(nv);

    // tmp_vecs.Append(make_shared<VVector<double>>(fes->GetNDof()));

    Array<int> dnums(1);
    for (auto i : Range(nv))
    {
      fes->GetDofNrs(NodeId(NT_VERTEX,i),dnums);    
      if (dnums.Size() > 0 && IsRegularDof(dnums[0]))
        (*vert2dof)[i] = dnums[0];
      else
        (*vert2dof)[i] = NO_DOF_NR;
    }
    // v2d_on_lvl.Append(vert2dof);
    //cout << "vert2dof : " << *vert2dof << endl;
    */
  }

  void P2Prolongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    if (fes == nullptr)
      throw Exception("call Update before prolongating");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");    
    static Timer t("Prolongate"); RegionTimer r(t);
    


    size_t nvC = nVertLevel[finelevel-1];
    size_t nvF = nVertLevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel-1]->FV<double>();
    fw = fv;
    fv = 0.0;

    cout << "update coarse verts" << endl;
    for (size_t i = 0; i < nvC; i++)
    {
        fv(i) = fw(i);
        cout << "vert i: " << i << ", val c: " << fw(i) << ", val f: " << fv(i) << endl << endl;
    }

    cout << "update fine verts" << endl;
    for (size_t i = nvC; i < nvF; i++)
      {        
        auto parents = ma->GetParentNodes (i);
        cout << "vert i: " << i << ", parents: " << parents[0] << ", " << parents[1] << endl;
        cout << "val p1: " << fw(parents[0]) << ", val p2" << fw(parents[1]) << endl;
        for (auto j : Range(2))
          fv( i ) += 0.5 * fw( parents[j] );
        
        int edgeParent = ma->GetParentEdge( i );
        cout << "   edgeparent: " << edgeParent << endl;
        cout << "   val edgeparent: " << fw(edgeParent) << endl;
        fv( i ) -= 0.125 * fw( edgeParent );
        cout << "updated node: " << fv(i) << endl << endl;
      }

    size_t neC = nEdgeLevel[finelevel-1];
    size_t neF = nEdgeLevel[finelevel];

    size_t nedges = nEdgeLevel[finelevel] - nEdgeLevel[finelevel-1];

    cout << "update fine edges, num: " << nedges << endl;
    for ( size_t i = 0; i < nedges; i++ )
    {
      size_t edgenum = neC + i;
      size_t unk = nvF + i;

      cout << "edge i: " << i << ", global edge num: " << edgenum << ", dof nr: " << unk << endl;

      auto edgeconn = ma->GetEdgeConn( edgenum );
      cout << "edgeconn: " << edgeconn[0] << ", " << edgeconn[1] << ", " << edgeconn[2] << endl;

      if( edgeconn[2] == -1 )
      {
        // cout << "single edge connection" << endl;        
        int parentEdge = ma->GetParentEdge( edgeconn[0] );
        // cout << "parent edge" << parentEdge << endl;
        // fv(unk) = 0.25 * fw( parentEdge );
        fv( unk ) = 0.25 * fw( edgeconn[0] );

        continue;
        cout << endl << endl << endl << "----------------------------------------------------" << endl;
        /*fv( unk ) = 0.25 * fw( edgeconn[0] );*/
        auto everts = ma->GetEdgePNums( edgenum );
        int newvert =1;
        int oldvert =0;
        if ( edgeconn[0] == everts[0] )
        {
          newvert = 0;
          oldvert = 1;
        }
        auto coarseparent = ma->GetParentNodes( edgeconn[0] );
        int commonVert = coarseparent[0];
        int detachedVert = coarseparent[1];

        if( commonVert != everts[oldvert] )
        {
          commonVert = coarseparent[1];
          detachedVert = coarseparent[0];
        }


        //int parentEdge = ma->GetParentEdge( edgeconn[0] );

        cout << "commonVert: " << commonVert << ", detached vert: " << detachedVert
             << ", earlier edge: " << parentEdge << endl;
        

        double val = -3./32 * fw( parentEdge ) + 0.75 * fw( commonVert ) + 0.25 * fw( detachedVert );
        val -= 0.5*( fv( everts[0]) + fv( everts[1] ) );
        fv(unk) = -8 * val;
        cout << "coarse val: " << fw( edgeconn[0] ) <<", fine val: " << fv(unk) << endl;
      }
      else if ( edgeconn[2] == -2 )
      {
        cout << "inner single edge" << endl;
        fv( unk ) = 0.25 * fw( edgeconn[0] );
      }
      else
      {
        double fac[3] = {-0.25,0.5,0.5};
        for (auto j: Range(3) )
        {
          cout << "edgeconn: " << edgeconn[j] << ", fac: " << fac[j] << endl;
          fv( unk ) += fac[j] * fw( edgeconn[j] );
          cout << "val coarse: " << fw( edgeconn[j] ) << endl;
        }
      }
      cout << endl;
    }



  }

}

