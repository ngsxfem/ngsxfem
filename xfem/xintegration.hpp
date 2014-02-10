#ifndef FILE_XINTEGRATION_HPP
#define FILE_XINTEGRATION_HPP

#include "../common/spacetimefe.hpp"   // for ScalarSpaceTimeFiniteElement

#include <set>
#include <vector>

using namespace ngfem;

namespace xintegration
{
  /// struct which defines the relation a < b for Point4DCL (in order to use std::set-features)
  template< int SD>
  struct Pointless {
    bool operator() (const Vec<SD> & a, const Vec<SD> & b) const {
      // if you want to merge point which are "the same up 14 digits..."
      const double EPS = 0.0; // 1e-14*abs(b[i])
      for (int i=0; i<SD; i++)
      {
        if (a[i] < b[i] - EPS)
          return true;
        if (a[i] > b[i] + EPS)
          return false;
      }
      return false;
    }
    typedef Vec<SD> first_argument_type;
    typedef Vec<SD> second_argument_type;
    typedef bool result_type;
  };


  /// Container set constitutes a collection of PointXDs 
  /// main feature: the operator()(const PointXDCL & p)
  /// The points in the container are owned and later 
  /// released by Point4DContainer
  template<int SD>
  class PointContainer
  {
    typedef  std::set<Vec<SD>, Pointless<SD> > SetOfPoints;
  protected:
    SetOfPoints pset;
#ifdef DEBUG
    size_t k;
#endif
  public: 
    PointContainer()
    {
#ifdef DEBUG
      k=0;
#endif
      pset.clear();
    };
    
    /// Access operator to points
    /// Either point is already in the Container, 
    ///   then return pointer to that point
    /// or point is not in the Container yet,
    ///   then add point to container and return pointer to new Point4D
    /// The return value (pointer ) points to a Point4D which is owned
    /// and later released by PointContainer
    const Vec<SD>* operator()(const Vec<SD> & p)
    {
      typename SetOfPoints::iterator it;
      it = pset.find(p);
      if (it == pset.end())
      {
        pset.insert(p);
        return &(*pset.find(p));
      }
      else
      {
#ifdef DEBUG
        k++;
#endif
        return &(*it);
      }
    }

    void Report(std::ostream & out) const
    {
      out << " PointContainer stored " << pset.size() << " points.\n";
#ifdef DEBUG
      out << " PointContainer rejected " << k << " points.\n";
#endif
    }

    ~PointContainer(){};
  };

  enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};

  ostream & operator<<(ostream & s, DOMAIN_TYPE dt)
  {
    switch (dt)
	{
	case NEG: 
      s << "NEG";
      break;
	case POS: 
      s << "POS";
      break;
	case IF: 
      s << "IF";
      break;
	default:
      ;
	}
  }

  template < int SD >
  struct QuadratureRule
  {
    Array < Vec<SD> > points;
    Array < double > weights;
    int Size() { return points.Size(); }
  };

  template < int SD >
  struct CompositeQuadratureRule
  {
    QuadratureRule<SD> quadrule_pos;
    QuadratureRule<SD> quadrule_neg;
    QuadratureRule<SD> quadrule_if;
    
    QuadratureRule<SD> & GetRule(DOMAIN_TYPE dt)
    {
      switch (dt)
      {
      case NEG: 
        return quadrule_neg;
        break;
      case POS: 
        return quadrule_pos;
        break;
      case IF: 
        return quadrule_if;
        break;
      default:
        throw Exception(" DOMAIN_TYPE not known ");
        return quadrule_if;
      }
    }
  };

  // does only make sense for D==3
  template < int D >
  inline double Determinant (const Vec<D> & col1,
			     const Vec<D> & col2,
			     const Vec<D> & col3)
  {
    return
      col1[0] * ( col2[1] * col3[2] - col2[2] * col3[1]) +
      col1[1] * ( col2[2] * col3[0] - col2[0] * col3[2]) +
      col1[2] * ( col2[0] * col3[1] - col2[1] * col3[0]);
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  class NumericalIntegrationStrategy
  {
  public:

    enum { D = ET_trait<ET_SPACE>::DIM };
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM };

    const ScalarSpaceTimeFEEvaluator<D> & lset;
    PointContainer<SD> & pc;

    Array< Vec<D> > verts_space;
    Array< double > verts_time;

    int ref_level_space = 0;
    int ref_level_time = 0;
    int int_order_space = 0;
    int int_order_time = 0;

    CompositeQuadratureRule<SD> & compquadrule;

    int GetIntegrationOrderMax() const{ return max(int_order_space, int_order_time); }

    
    NumericalIntegrationStrategy(const NumericalIntegrationStrategy & a, int reduce_ref_space = 0, int reduce_ref_time = 0)
      : lset(a.lset), pc(a.pc), 
        int_order_space(a.int_order_space), int_order_time(a.int_order_time),
        ref_level_space(a.ref_level_space-reduce_ref_space), ref_level_time(a.ref_level_time-reduce_ref_time),
        compquadrule(a.compquadrule)
    {
    }

    NumericalIntegrationStrategy(const ScalarSpaceTimeFEEvaluator<D> & a_lset, PointContainer<SD> & a_pc, 
                                 CompositeQuadratureRule<SD> & a_compquadrule,
                                 int a_int_order_space = 2, int a_int_order_time = 2, 
                                 int a_ref_level_space = 0, int a_ref_level_time = 0 )
      : lset(a_lset), pc(a_pc), 
        ref_level_space(a_ref_level_space), ref_level_time(a_ref_level_time),
        int_order_space(a_int_order_space), int_order_time(a_int_order_time),
      compquadrule(a_compquadrule)
    {
      ;
    }
    
    void SetVerticesSpace(const Array<Vec<D> > & verts)
    {
      verts_space.SetSize(verts.Size());
      for (int i = 0; i < verts.Size(); ++i)
        verts_space[i] = verts[i];
    }

    void SetVerticesSpace()
    {
      const POINT3D * verts = ElementTopology::GetVertices(ET_SPACE);
      const int nv = ElementTopology::GetNVertices(ET_SPACE);

      verts_space.SetSize(nv);
      for (int i = 0; i < nv; ++i)
      {
        Vec<D> newvert;
        for (int d = 0; d < D; ++d)
          verts_space[i][d] = verts[i][d];
      }
    }

    void SetVerticesTime(const Array<double> & verts)
    {
      verts_time.SetSize(verts.Size());
      for (int i = 0; i < verts.Size(); ++i)
        verts_time[i] = verts[i];
    }

    void SetVerticesTime()
    {
      const int np1dt = pow(2,ref_level_time);
      const double ht = 1.0 / np1dt;
      verts_time.SetSize(np1dt+1);
      for (int i = 0; i < np1dt + 1; ++i)
      {
        verts_time[i] = i * ht;
      }
    }

    void SetVerticesTimeFromUpperHalf(const Array< double >& verts_t)
    {
      const int newsize = (verts_t.Size()+1)/2;
      const int offset = (verts_t.Size()-1)/2;
      verts_time.SetSize(newsize);
      for (int i = 0; i < newsize; ++i)
      {
        verts_time[i] = verts_t[offset+i];
      }
    }

    void SetVerticesTimeFromLowerHalf(const Array< double >& verts_t)
    {
      const int newsize = (verts_t.Size()+1)/2;
      verts_time.SetSize(newsize);
      for (int i = 0; i < newsize; ++i)
      {
        verts_time[i] = verts_t[i];
      }
    }



  };

  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE CheckIfCut(const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {

    bool haspos = false;
    bool hasneg = false;

    int np1ds = pow(2,numint.ref_level_space);
    int np1dt = pow(2,numint.ref_level_time);

    double dx_scalar = 1.0 / np1ds;
    
    Array<double> time (0);
    double ht = 1.0;

    switch (ET_SPACE)
    {
    case ET_TRIG:
    case ET_TET:
    {
      int sum = 0;
      INT< ET_trait<ET_SPACE>::DIM > I;
      Vec< ET_trait<ET_SPACE>::DIM + 1> position;
      for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
        I[i] = 0;

      // cout << " index = ";
      // for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
      //   cout << I[i] << ", \t";
      // cout << endl;

      bool finish = false;
      while (finish == false)
      {
        // do something

        //calculate all points corresponding to the current space position
        // loop over time points
        for (int i = 0; i < np1dt + 1; ++i)
        {
          for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
            position[d] = numint.verts_space[0][d];

          for (int j = 0; j < ET_trait<ET_SPACE>::DIM; ++j)
          {
            for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
            {
               position[d] += I[j] * dx_scalar * (numint.verts_space[j+1][d] - numint.verts_space[0][d]);
            }
          }
          // for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
          //   cout << position[d] << ",\t";

          if (ET_TIME == ET_SEGM)
          {
            position[ET_trait<ET_SPACE>::DIM] = numint.verts_time[i];
            // cout << position[ET_trait<ET_SPACE>::DIM] << ",\t";
          }
          const ngfem::ScalarSpaceTimeFEEvaluator<ET_trait<ET_SPACE>::DIM> & eval (numint.lset);
          const double lsetval = eval(position);
          if (lsetval >= 0.0)
            haspos = true;
          else
            hasneg = true;
          // cout << " :: " << lsetval << ",\t";
          // cout << endl;
          
          if(haspos && hasneg)
          {
            // finish = true;
            // break;
            // cout << " IF " << endl;
            return IF;
          }
        }
        
        I[0]++;
        int sum = 0;
        for (int checkdim = 0; checkdim < ET_trait<ET_SPACE>::DIM; ++checkdim)
        {
          sum = 0;
          for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
            sum += I[i];

          if (sum >= np1ds + 1)
          {
            if ( checkdim == ET_trait<ET_SPACE>::DIM - 1)
            {
              finish = true;
              break;
            }
            else
            {
              I[checkdim] = 0;
              I[checkdim+1]++;
            }
          }
          else
            break;
        }
        // if (!finish)
        // {
        //   cout << " index = ";
        //   for (int i = 0; i < ET_trait<ET_SPACE>::DIM; ++i)
        //     cout << I[i] << ", \t";
        //   cout << endl;
        // }
      }
      if (haspos)
        return POS;
      else
        return NEG;
    }
    break;
    default:  // for the compiler
      throw Exception(" this CheckIfCut::ELEMENT_TYPE is not treated... yet ");
      break;
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE MakeQuadRule(const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {

    const int D = ET_trait<ET_SPACE>::DIM;
    const int SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM;

    DOMAIN_TYPE dt_self = CheckIfCut(numint);

    if (dt_self == IF)
    {
      bool refine_time = numint.ref_level_time > 0;
      bool refine_space = numint.ref_level_space > 0;

      if (refine_time && refine_space)
      {
        if (numint.ref_level_time >= numint.ref_level_space)
        {
          refine_time = true;
          refine_space = false;
        }
        else
        {
          refine_time = false;
          refine_space = true;
        }
      }
      
      if (refine_time)
      {
        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_upper (numint, 0, 1);
        numint_upper.SetVerticesSpace(numint.verts_space);
        numint_upper.SetVerticesTimeFromUpperHalf(numint.verts_time);
        // cout << " call upper " << endl;
        DOMAIN_TYPE dt_upper = MakeQuadRule(numint_upper);

        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_lower (numint, 0, 1);
        numint_lower.SetVerticesSpace(numint.verts_space);
        numint_lower.SetVerticesTimeFromLowerHalf(numint.verts_time);
        // cout << " call lower " << endl;
        DOMAIN_TYPE dt_lower = MakeQuadRule(numint_lower);

        // cout << " dt_upper: " << dt_upper << endl;
        // cout << " dt_lower: " << dt_lower << endl;

        // if (dt_upper == IF || dt_lower == IF)
        //   return IF;
        // else
        //   if (dt_upper == dt_lower)
        //   {
        //     throw Exception(" sollte nicht passieren ");
        //     return dt_upper;
        //   }
        //   else
        //   {
        //     throw Exception(" tut tut - sollte auch nicht passieren");
        //     return dt_upper;
        //   }
      }
      else if (refine_space)
      {
        if ( D == 2)
        {
          Array< Vec<3> > baryc(6);
          Array< INT<3> > trigs(4);
          
          baryc[0] = Vec<3>(0.0,0.0,1.0);
          baryc[1] = Vec<3>(0.5,0.0,0.5);
          baryc[2] = Vec<3>(1.0,0.0,0.0);
          baryc[3] = Vec<3>(0.0,0.5,0.5);
          baryc[4] = Vec<3>(0.5,0.5,0.0);
          baryc[5] = Vec<3>(0.0,1.0,0.0);

          trigs[0] = INT<3>(0,1,3);
          trigs[1] = INT<3>(1,2,4);
          trigs[2] = INT<3>(1,3,4);
          trigs[3] = INT<3>(3,4,5);
          
          for (int i = 0; i < 4; ++i)
          {
            NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_i (numint, 1, 0);
            numint_i.SetVerticesTime(numint.verts_time);
            
            Array< Vec<D> > newverts(3);

            for (int j = 0; j < 3; ++j) //vertices
            {
              newverts[j] = Vec<D>(0.0);
              for (int d = 0; d < 3; ++d) 
                newverts[j] += baryc[trigs[i][j]][d] * numint.verts_space[d];
            }

            // cout << " old verts \n";
            // for (int j = 0; j < 3; ++j) //vertices
            //   cout << " vert(" << j << "):" << numint.verts_space[j] << endl;
            // cout << endl;
            // cout << " new verts \n";
            // for (int j = 0; j < 3; ++j) //vertices
            //   cout << " new vert(" << j << "):" << newverts[j] << endl;
            // cout << endl;

            numint_i.SetVerticesSpace(newverts);
            // getchar();
            // cout << " call s(" << i << ")" << endl;
            DOMAIN_TYPE dt_i = MakeQuadRule(numint_i);
            // cout << " dt_s(" << i << "): " << dt_i << endl;

          }
          
        }
        else
          throw Exception(" tut tut - 3D not yet implemented");
      }
      // cout << "TODO: Call Decomposition into uncut simplices" << endl;
      return IF;
    }
    else
    {
      double trafofac = 1.0; 

      if (D==2)
      {
        Vec<D> a = numint.verts_space[1] - numint.verts_space[0];
        Vec<D> b = numint.verts_space[2] - numint.verts_space[0];
        trafofac = abs(a(0) * b(1) - a(1) * b(0));
      }
      else
      {
        Vec<D> a = numint.verts_space[1] - numint.verts_space[0];
        Vec<D> b = numint.verts_space[2] - numint.verts_space[0];
        Vec<D> c = numint.verts_space[3] - numint.verts_space[0];
        trafofac = abs(Determinant(a,b,c));
      }

      // cout << " trafofac = " << trafofac << endl;
      const double dt = numint.verts_time[numint.verts_time.Size()-1] - numint.verts_time[0];
      // cout << " dt = " << dt << endl;
      const double t0 = numint.verts_time[0];

      const IntegrationRule & ir_space = SelectIntegrationRule (ET_TIME, numint.int_order_time);
      const IntegrationRule & ir_time = SelectIntegrationRule (ET_SPACE, numint.int_order_space);
      

      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        double current = t0 + ir_time[l](0) * dt;
        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          Vec<ET_trait<ET_SPACE>::DIM> point(0.0);
          point = numint.verts_space[0];
          for (int m = 0; m < D ;++m)
            point += ir_space[k](m) * numint.verts_space[m];

          const double weight = ir_time[l].Weight() * dt * ir_space[k].Weight() * trafofac;
          Vec<ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM> ipoint(0.0);
          
          for (int m = 0; m < D ;++m)
            ipoint(m) = point(m);
          if (D != SD)
            ipoint(SD-1) = current;
          
          numint.compquadrule.GetRule(dt_self).points.Append(ipoint);
          numint.compquadrule.GetRule(dt_self).weights.Append(weight);

        }
      }
      // fill integration rule already ? 
      return dt_self;
    }
  }


  template <int D>
  class Simplex
  {
  protected:
  public:
    bool cut;
    Array< const Vec<D> * > p;
    Simplex(const Array< const Vec<D> * > & a_p): p(a_p)
    {
      ;
    }
  };

  template <int D>
  inline ostream & operator<< (ostream & ost, const Simplex<D> & s)
  {
    for (int i = 0; i < D+1; i++)
      ost << i << ":" << *(s.p[i]) << "\t";
    ost << endl;
    return ost;
  }

  // Decompose the geometry K = T x I with T \in {trig,tet} and I \in {segm, point} into simplices of corresponding dimensions
  // D is the dimension of the spatial object, SD is the dimension of the resulting object T
  template<int D, int SD>
  void DecomposeIntoSimplices(ELEMENT_TYPE et_space, ELEMENT_TYPE et_time, 
                              Array<Simplex<SD> *>& ret, PointContainer<SD> & pc, LocalHeap & lh)
  {
    ret.SetSize(0);
    cout << " et_space : " << et_space << endl;
    cout << " et_time : " << et_space << endl;
    if (et_time == ET_SEGM)
    {
      switch (et_space)
	  {
	  case ET_TRIG:
	  case ET_TET:
      {
        const POINT3D * verts = ElementTopology::GetVertices(et_space);
        Array< const Vec<SD> * > p(2*(SD));
        for (int i = 0; i < SD; ++i)
        {
          Vec<SD> newpoint;
          Vec<SD> newpoint2;
          for (int d = 0; d < D; ++d)
          {
            newpoint[d] = verts[i][d];
            newpoint2[d] = verts[i][d];
          }
          newpoint[D] = 0.0;
          newpoint2[D] = 1.0;

          p[i] = pc(newpoint);
          p[i+SD] = pc(newpoint2);
        }

        Array< const Vec<SD> * > tet(SD+1);
        for (int i = 0; i < SD; ++i)
        {
          for (int j = 0; j < SD+1; ++j)
            tet[j] = p[i+j];
          ret.Append(new (lh) Simplex<SD> (tet));
        }

        cout << " report \n";
        for (int i = 0; i < ret.Size(); ++i)
          cout << " simplex " << i << ": " << *ret[i] << endl;
      }
      break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }
    }
    else if (et_time == ET_POINT)
    {
      switch (et_space)
	  {
	  case ET_TRIG:
	  case ET_TET:
      {
        const POINT3D * verts = ElementTopology::GetVertices(et_space);
        Array< const Vec<SD> * > tet(D+1);
        for (int i = 0; i < D+1; ++i)
        {
          Vec<SD> newpoint;
          for (int d = 0; d < D; ++d)
          {
            newpoint[d] = verts[i][d];
          }
          tet[i] = pc(newpoint);
        }

        ret.Append(new (lh) Simplex<SD> (tet));

        cout << " report \n";
        for (int i = 0; i < ret.Size(); ++i)
          cout << " simplex " << i << ": " << *ret[i] << endl;
      }
      break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }
    }
    else
      throw Exception(" this ELEMENT_TYPE et_time is not treated... yet ");
  }

  template <int D, int SD>
  void EvaluateLevelset (const ScalarSpaceTimeFEEvaluator<D> & lset, Array<Simplex<SD> *>& ret)
  {
    for (int i = 0; i < ret.Size(); ++i)
    {
      Simplex<SD> & simp = *ret[i];
      for (int j = 0; j < SD+1; ++j)
      {

      }
    }
  }


  template<int D>
  void Decompose(ELEMENT_TYPE et, const ScalarSpaceTimeFEEvaluator<D> & stfeeval)
  {

    
	switch (et)
	  {
	  case ET_TRIG:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
	  case ET_QUAD:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
	  case ET_TET:
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
	    break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }


  }
}

#endif
