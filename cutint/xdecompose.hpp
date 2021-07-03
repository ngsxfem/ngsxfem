#ifndef FILE_XDECOMPOSE_HPP
#define FILE_XDECOMPOSE_HPP

#include "../cutint/fieldeval.hpp"
#include "../utils/ngsxstd.hpp"

#include <set>
#include <vector>

using namespace ngfem;

namespace xintegration
{

  template<int SD> class PointContainer;

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

    Simplex(Simplex<D> & a_s): p(a_s.p)
    {
      ;
    }

    DOMAIN_TYPE CheckIfCut(const ScalarFieldEvaluator & lset) const
    {
      static Timer timer("Simplex::CheckifCut"); RegionTimer reg (timer);

      bool haspos = false;
      bool hasneg = false;

      double sumneg = 0.0;
      double sumpos = 0.0;
      double sumposneg = 0.0;

      for (int i = 0; i < D+1; ++i)
      {
        const double lsetval = lset(*(p[i]));
        if (lsetval >= 0.0)
        {
          sumpos += lsetval;
          haspos = true;
        }
        else
        {
          sumneg -= lsetval;
          hasneg = true;
        }
      }

      sumposneg = sumpos+sumneg;
      sumpos/=sumposneg;
      sumneg/=sumposneg;

      if (sumpos < 1e-14)
        haspos = false;

      if (sumneg < 1e-14)
        hasneg = false;

      if(haspos && hasneg)
      {
        return IF;
      }
      else
        if (haspos)
          return POS;
        else
          if (hasneg)
            return NEG;
          else
          {
            throw Exception(" this is not possible, is it?");
            return IF;
          }
    }
  };

  /// helper function:
  /// Calculate Determinant of a Matrix (used for Transformation Weights)
  template < int D >
  inline double Determinant (const Vec<D> & col1,
			     const Vec<D> & col2,
			     const Vec<D> & col3)
  {
    if (D==3)
      return
        col1[0] * ( col2[1] * col3[2] - col2[2] * col3[1]) +
        col1[1] * ( col2[2] * col3[0] - col2[0] * col3[2]) +
        col1[2] * ( col2[0] * col3[1] - col2[1] * col3[0]);
    else
      return
        col1(0) * col2(1) - col1(1) * col2(0);
  }

  template<int D, int SD>
  double Measure(const Array< const Vec<SD> *> & s);
  template<>
  double Measure<1,1>(const Array< const Vec<1> *> & s);
  template<>
  double Measure<1,2>(const Array< const Vec<2> *> & s);
  template<>
  double Measure<2,2>(const Array< const Vec<2> *> & s);
  template<>
  double Measure<2,3>(const Array< const Vec<3> *> & s);
  template<>
  double Measure<3,3>(const Array< const Vec<3> *> & s);


  template <int D>
  inline ostream & operator<< (ostream & ost, const Simplex<D> & s)
  {
    for (int i = 0; i < D+1; i++)
      ost << i << ":" << *(s.p[i]) << "\t";
    ost << endl;
    return ost;
  }

  // Decompose the geometry K = T x I with T \in {trig,tet} and I \in {segm, point} into simplices of corresponding dimensions
  template <int SD>
  void DecomposePrismIntoSimplices(Array<const Vec<SD> *> & verts,
                                    Array<Simplex<SD> *>& ret, 
                                    PointContainer<SD> & pc, 
                                    LocalHeap & lh)
  {
    static int timer = NgProfiler::CreateTimer ("DecomposePrismIntoSimplices"); NgProfiler::RegionTimer reg (timer);

    ret.SetSize(SD);
    Array< const Vec<SD> * > tet(SD+1);
    for (int i = 0; i < SD; ++i)
    {
      for (int j = 0; j < SD+1; ++j)
        tet[j] = verts[i+j];
      ret[i] = new Simplex<SD> (tet);
    }
  }

}

#endif
