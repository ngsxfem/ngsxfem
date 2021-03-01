#include "../utils/xprolongation.hpp"

using namespace ngcomp;

//#define __HASP2PROL__

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

    static int timer = NgProfiler::CreateTimer ("Prolongate");
    NgProfiler::RegionTimer reg (timer);

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
    static int timer = NgProfiler::CreateTimer ("Restrict");
    NgProfiler::RegionTimer reg (timer);

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
        if (IsRegularDof(vert2dof_coarse[parents[j]])) // not necessary... apperently it is!
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
    #ifdef __HASP2PROL__

    if (fes == nullptr)
      throw Exception("call Update before prolongating");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");    
    static int timer = NgProfiler::CreateTimer ("Prolongate");
    NgProfiler::RegionTimer reg (timer);
    


    size_t nvC = nVertLevel[finelevel-1];
    size_t nvF = nVertLevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel-1]->FV<double>();
    fw = fv;
    fv = 0.0;

    // update coarse verts
    for (size_t i = 0; i < nvC; i++)
        fv(i) = fw(i);

    // update fine verts
    for (size_t i = nvC; i < nvF; i++)
      {        
        auto parents = ma->GetParentNodes (i);        
        for (auto j : Range(2))
          fv( i ) += 0.5 * fw( parents[j] );
        
        // p2 edge update
        int edgeParent = ma->GetParentEdge( i );        
        fv( i ) -= 0.125 * fw( edgeParent );
      }

    size_t neC = nEdgeLevel[finelevel-1];
    size_t neF = nEdgeLevel[finelevel];

    size_t nedges = nEdgeLevel[finelevel] - nEdgeLevel[finelevel-1];

    // update new edges
    for ( size_t i = 0; i < nedges; i++ )
    {
      // loca edge number: i
      // global edge number: coarse level edges + i
      // dof num: first vertices, then edges: nvF + i
      size_t edgenum = neC + i;
      size_t unk = nvF + i;      

      auto edgeconn = ma->GetEdgeConn( edgenum );

      if( edgeconn[2] == -1 )
      {
        // cout << "single edge connection" << endl;
        int parentEdge = ma->GetParentEdge( edgeconn[0] );       
        fv( unk ) = 0.25 * fw( edgeconn[0] );
        
      }
      else
      {
        //special case inner edge in uniform refinement
        double fac[3] = {-0.25,0.5,0.5};
        for (auto j: Range(3) )
        {
          fv( unk ) += fac[j] * fw( edgeconn[j] );
        }
      }
      
    }
    #endif
  }

  void P2Prolongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    #ifdef __HASP2PROL__
    size_t nvC = nVertLevel[finelevel-1];
    size_t nvF = nVertLevel[finelevel];

    size_t neC = nEdgeLevel[finelevel-1];
    size_t neF = nEdgeLevel[finelevel];

    size_t nedges = neF - neC;

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel]->FV<double>();

    
    fw = fv;
    fv = 0.0;


    for (size_t i = 0; i < nvC; ++i)
      fv(i) = fw(i);

    for(size_t i = nvC; i < nvF; ++i)
    {
      auto parents = ma->GetParentNodes(i);
      int edgeParent = ma->GetParentEdge( i ); 

      for (auto j : Range(2))
        fv( parents[j] ) += 0.5 * fw( i );
      
      fv( edgeParent ) -= 0.125 * fw( i );
    }


    // edge factors
    double fac[3] = {-0.25,0.5,0.5};
    // update new edges
    for ( size_t i = 0; i < nedges; i++ )
    {

      size_t edgenum = neC + i;
      size_t unk = nvF + i;

      auto edgeconn = ma->GetEdgeConn( edgenum );

      if( edgeconn[2] == -1 )
      {
        fv( edgeconn[0] ) += 0.25 * fw( unk );
      }
      else
        for (auto j: Range(3) )
        {
          fv( edgeconn[j] ) += fac[j] * fw( unk );
        }

    }
    #endif

  }

   void P2CutProlongation :: Update (const FESpace & afes)
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


    int nv = ma->GetNV();
    int nlvl = nEdgeLevel.Size();
    cout << "nlvl: " << nlvl << endl;
    size_t neC;

    if ( nlvl == 1 )    
      neC = 0;
    else
      neC = nEdgeLevel[nlvl-2];

    size_t neF = nEdgeLevel[nlvl-1];

    size_t nedges = neF - neC;

    int ndofs = nv + nedges;

    shared_ptr<Array<int>> vert2dof = make_shared<Array<int>>(ndofs);

    tmp_vecs.Append(make_shared<VVector<double>>(fes->GetNDof()));

    Array<int> dnums(1);

    // vertices
    for (auto i : Range(nv))
    {
      fes->GetDofNrs(NodeId(NT_VERTEX,i),dnums);    
      if (dnums.Size() > 0 && IsRegularDof(dnums[0]))
        (*vert2dof)[i] = dnums[0];
      else
        (*vert2dof)[i] = NO_DOF_NR;
      
    }
    // edges
    for ( size_t i = 0; i < nedges; i++ )
    {

      size_t edgenum = neC + i;
      size_t unk = nv + i;
      fes->GetDofNrs(NodeId(NT_EDGE,edgenum),dnums);    
      if (dnums.Size() > 0 && IsRegularDof(dnums[0]))
        (*vert2dof)[unk] = dnums[0];
      else
        (*vert2dof)[unk] = NO_DOF_NR;      
    }

    v2d_on_lvl.Append(vert2dof);
    // cout << "vert2dof : " << *vert2dof << endl;
  }

  void P2CutProlongation :: ProlongateInline (int finelevel, BaseVector & v) const
  {
    #ifdef __HASP2PROL__
    if (fes == nullptr)
      throw Exception("call Update before prolongating");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");    
    static int timer = NgProfiler::CreateTimer ("Prolongate");
    NgProfiler::RegionTimer reg (timer);



    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);    


    size_t nvC = nVertLevel[finelevel-1];
    size_t nvF = nVertLevel[finelevel];

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel-1]->FV<double>();
    fw = fv;
    fv = 0.0;

    // update coarse verts
    for (size_t i = 0; i < nvC; i++)
        if (IsRegularDof(vert2dof_fine[i]))
          fv(vert2dof_fine[i]) = fw(vert2dof_coarse[i]);

    // update fine verts
    for (size_t i = nvC; i < nvF; i++)
      {
        if (!IsRegularDof(vert2dof_fine[i]))
          continue;

        auto parents = ma->GetParentNodes (i);        
        for (auto j : Range(2))
          fv( vert2dof_fine[i] ) += 0.5 * fw( vert2dof_coarse[parents[j]] );
        
        // p2 edge update
        int edgeParent = vert2dof_coarse[ ma->GetParentEdge( i ) ];
        fv(  vert2dof_fine[i] ) -= 0.125 * fw( edgeParent );
      }

    size_t neC = nEdgeLevel[finelevel-1];
    size_t neF = nEdgeLevel[finelevel];

    size_t nedges = nEdgeLevel[finelevel] - nEdgeLevel[finelevel-1];

    // update new edges
    double fac[3] = {-0.25,0.5,0.5};
    for ( size_t i = 0; i < nedges; i++ )
    {
      // loca edge number: i
      // global edge number: coarse level edges + i
      // dof num: first vertices, then edges: nvF + i
      size_t edgenum = neC + i;
      size_t unk = nvF + i;

      if (!IsRegularDof(vert2dof_fine[unk]))
          continue;

      auto edgeconn = ma->GetEdgeConn( edgenum );

      if( edgeconn[2] == -1 )
        fv( vert2dof_fine[unk] ) = 0.25 * fw( vert2dof_coarse[ edgeconn[0] ] );
      //special case inner edge in uniform refinement
      else
        for (auto j: Range(3) )        
          fv( vert2dof_fine[unk] ) += fac[j] * fw( vert2dof_coarse[ edgeconn[j] ] );
    }
    #endif

  }

void P2CutProlongation :: RestrictInline (int finelevel, BaseVector & v) const
  {
    #ifdef __HASP2PROL__
    if (fes == nullptr)
      throw Exception("call Update before restricting");
    if (v.EntrySize() > 1)
      throw Exception("no dim>1 yet");
    static int timer = NgProfiler::CreateTimer ("Restrict");
    NgProfiler::RegionTimer reg (timer);

    size_t nvC = nVertLevel[finelevel-1];
    size_t nvF = nVertLevel[finelevel];

    size_t neC = nEdgeLevel[finelevel-1];
    size_t neF = nEdgeLevel[finelevel];

    size_t nedges = neF - neC;

    FlatVector<> fv = v.FV<double>();
    FlatVector<> fw = tmp_vecs[finelevel]->FV<double>();

    Array<int> & vert2dof_fine = *(v2d_on_lvl[finelevel]);
    Array<int> & vert2dof_coarse = *(v2d_on_lvl[finelevel-1]);
    
    fw = fv;
    fv = 0.0;


    for (size_t i = 0; i < nvC; ++i)
    {
      if (!IsRegularDof(vert2dof_coarse[i]))
        continue;
      if (IsRegularDof(vert2dof_fine[i]))
        fv(vert2dof_coarse[i]) = fw(vert2dof_fine[i]);
    }

    for(size_t i = nvC; i < nvF; ++i)
    {
      if (!IsRegularDof(vert2dof_fine[i]))
        continue;

      auto parents = ma->GetParentNodes(i);
      int edgeParent = ma->GetParentEdge( i ); 

      for (auto j : Range(2))
        fv( vert2dof_coarse[ parents[j] ] ) += 0.5 * fw( vert2dof_fine[i] );
      
      fv( vert2dof_coarse[edgeParent] ) -= 0.125 * fw( vert2dof_fine[i] );
    }


    // edge factors
    double fac[3] = {-0.25,0.5,0.5};
    // update new edges
    for ( size_t i = 0; i < nedges; i++ )
    {
      size_t edgenum = neC + i;
      size_t unk = nvF + i;

      if (!IsRegularDof(vert2dof_fine[unk]))
        continue;

      auto edgeconn = ma->GetEdgeConn( edgenum );

      if( edgeconn[2] == -1 )
      {
        fv( vert2dof_coarse[ edgeconn[0] ] ) += 0.25 * fw( vert2dof_fine[unk] );
      }
      else
        for (auto j: Range(3) )
        {
          fv( vert2dof_coarse[ edgeconn[j] ] ) += fac[j] * fw( vert2dof_fine[unk] );
        }

    }
    #endif
  }
  

}

