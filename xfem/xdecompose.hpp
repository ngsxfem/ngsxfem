#ifndef FILE_XDECOMPOSE_HPP
#define FILE_XDECOMPOSE_HPP

#include "../common/spacetimefe.hpp"   // for ScalarSpaceTimeFiniteElement

#include <set>
#include <vector>

using namespace ngfem;

namespace xintegration
{

  template<int SD> class PointContainer;

  /// domain types: two domains: POS/NEG and the diving interface IF
  enum DOMAIN_TYPE { POS = 0, NEG = 1, IF = 2};

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

    DOMAIN_TYPE CheckIfCut(const ScalarSpaceTimeFEEvaluator<D-1> & lset) const
    {
      bool haspos = false;
      bool hasneg = false;

      for (int i = 0; i < D+1; ++i)
      {
        const double lsetval = lset(*(p[i]));
        if (lsetval >= 0.0)
          haspos = true;
        else
          hasneg = true;
      }

      if(haspos && hasneg)
      {
        return IF;
      }
      else
        if (haspos)
          return POS;
        else
          return NEG;
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
  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void DecomposeIntoSimplices(Array<Vec<ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM> > & verts,
                              Array<Simplex<ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM> *>& ret, 
                              PointContainer<ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM> & pc, 
                              LocalHeap & lh)
  {
    enum { D = ET_trait<ET_SPACE>::DIM };
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM};

    ret.SetSize(0);
    // cout << " et_space : " << ET_SPACE << endl;
    // cout << " et_time : " << ET_TIME << endl;
    if (ET_TIME == ET_SEGM)
    {
      switch (ET_SPACE)
	  {
	  case ET_TRIG:
	  case ET_TET:
      {
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

        // cout << " report \n";
        // for (int i = 0; i < ret.Size(); ++i)
        //   cout << " simplex " << i << ": " << *ret[i] << endl;
      }
      break;
      default:  // for the compiler
        throw Exception(" this ELEMENT_TYPE is not treated... yet ");
        break;
	  }
    }
    else if (ET_TIME == ET_POINT)
    {
      switch (ET_SPACE)
	  {
	  case ET_TRIG:
	  case ET_TET:
      {
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

        // cout << " report \n";
        // for (int i = 0; i < ret.Size(); ++i)
        //   cout << " simplex " << i << ": " << *ret[i] << endl;
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
