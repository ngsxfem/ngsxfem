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


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  class NumericalIntegrationStrategy
  {
  public:

    enum { D = ET_trait<ET_SPACE>::DIM };
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM };

    const ScalarSpaceTimeFEEvaluator<D> & lset;
    PointContainer<SD> & pc;

    Array< Vec<D> * > verts_space;
    Array< double > verts_time;

    int ref_level_space = 0;
    int ref_level_time = 0;
    int int_order_space = 0;
    int int_order_time = 0;

    int GetIntegrationOrderMax() const{ return max(int_order_space, int_order_time); }

    
    NumericalIntegrationStrategy(const NumericalIntegrationStrategy & a)
      : lset(a.lset), pc(a.pc), 
        int_order_space(a.int_order_space), int_order_time(a.int_order_time),
        ref_level_space(a.ref_level_space), ref_level_time(a.ref_level_time)
    {
    }

    NumericalIntegrationStrategy(const ScalarSpaceTimeFEEvaluator<D> & a_lset, PointContainer<SD> & a_pc, 
                                 int a_int_order_space = 2, int a_int_order_time = 2, 
                                 int a_ref_level_space = 0, int a_ref_level_time = 0 )
      : lset(a_lset), pc(a_pc), 
        ref_level_space(a_ref_level_space), ref_level_time(a_ref_level_time),
        int_order_space(a_int_order_space), int_order_time(a_int_order_time)
    {
      ;
    }
    
    void SetVerticesSpace(const Array<Vec<D> *> verts)
    {
      verts_space.SetSize(verts.Size());
      for (int i = 0; i < verts.Size(); ++i)
        verts_space[i] = verts[i];
    }

    void SetVerticesTime(const Array<double> verts)
    {
      verts_time.SetSize(verts.Size());
      for (int i = 0; i < verts.Size(); ++i)
        verts_time[i] = verts[i];
    }
  };

  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;

  enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};

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
    if (ET_TIME == ET_SEGM)
    {
      ht = 1.0 / np1dt;
      time.SetSize(np1dt+1);
      for (int i = 0; i < np1dt + 1; ++i)
      {
        time[i] = i * ht;
      }
    }
    else
    {
      
    }
    
    const POINT3D * verts = ElementTopology::GetVertices(ET_SPACE);

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
            position[d] = verts[0][d]; 

          for (int j = 0; j < ET_trait<ET_SPACE>::DIM; ++j)
          {
            for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
            {
               position[d] += I[j] * dx_scalar * (verts[j+1][d]-verts[0][d]);
            }
          }
          // for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
          //   cout << position[d] << ",\t";

          if (ET_TIME == ET_SEGM)
          {
            position[ET_trait<ET_SPACE>::DIM] = time[i];
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
