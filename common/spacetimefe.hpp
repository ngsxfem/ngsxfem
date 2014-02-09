#ifndef FILE_SPACETIMEFE_HPP
#define FILE_SPACETIMEFE_HPP

/// from ngsolve
#include <comp.hpp>
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

using namespace ngsolve;

namespace ngfem
{

typedef std::pair<double,double> TimeInterval;

/**
     
 */
  class SpaceTimeFiniteElement : public FiniteElement
  {
  protected:
    const FiniteElement & base_space;
    const FiniteElement & base_time;
    int order_space;
    int order_time;
    int ndof_space;
    int ndof_time;
  public:
    SpaceTimeFiniteElement(const FiniteElement & base_space,
                           const FiniteElement & base_time);
    virtual ~SpaceTimeFiniteElement();

    /// the name
    virtual string ClassName(void) const;

    int GetNDofSpace() const { return ndof_space;}
    int GetNDofTime() const { return ndof_time;}
    int OrderSpace() const { return order_space;}
    int OrderTime() const { return order_time;}
    virtual int GetNDof() const { return ndof;}

  };

  template <int D>  
  class ScalarSpaceTimeFiniteElement : public SpaceTimeFiniteElement
  {
  protected:
    const ScalarFiniteElement<D> & scalar_space;
    const DGFiniteElement<1> & scalar_time;
  public:
    ScalarSpaceTimeFiniteElement(const ScalarFiniteElement<D> & base_space,
                                 const DGFiniteElement<1> & base_time);
    virtual ~ScalarSpaceTimeFiniteElement();

    virtual int GetNDof() const { return ndof;}

    /// the name
    virtual string ClassName(void) const;

    virtual bool IsDGFiniteElement() const;
    /// mass diag
    virtual void GetDiagMassMatrix(FlatVector<> diagmass,
                                   LocalHeap & lh) const;

    /// compute shape
    virtual void CalcShapeTime (double time,
                                FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpace (const IntegrationPoint & ip,
                                 FlatVector<> shape) const;
    /// compute shape
    virtual void CalcShapeSpaceTime (const IntegrationPoint & ip, double time,
                                     FlatVector<> shape, LocalHeap & lh) const;

    /// compute dshape, vector: ndof 
    virtual void CalcDtShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatVector<> dshape, LocalHeap & lh) const;

    virtual ELEMENT_TYPE ElementType() const{ return base_space.ElementType(); }

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDxShapeSpaceTime (const IntegrationPoint & ip, double time,
                                       FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDxShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time,
                                             FlatMatrixFixWidth<D> dshape, LocalHeap & lh) const;
  };
    
    // /// compute dshape, vector: ndof 
    // virtual void CalcMappedDtShapeSpaceTime (const SpecificIntegrationPoint<D,D> & sip, double time
    //                                          FlatVector<> dshape) const;

    // /// compute dshape, matrix: ndof x (spacedim spacedim)
    // // virtual void CalcDDShape (const IntegrationPoint & ip, 
    // //                           FlatMatrix<> ddshape) const;

  template <int D>
  class ScalarSpaceTimeFEEvaluator
  {
  protected:
    const ScalarSpaceTimeFiniteElement<D> & fe;
    FlatVector<> linvec;
    IntegrationPoint ip;
    FlatVector<> shape;
    LocalHeap & lh;
  public:
    ScalarSpaceTimeFEEvaluator(const ScalarSpaceTimeFiniteElement<D> & a_fe, FlatVector<> a_linvec, LocalHeap & a_lh)
      : fe(a_fe), linvec(a_linvec), lh(a_lh), shape(linvec.Size(), lh)
    {
    }
      
    double operator()(const Vec<D+1>& point);
  };



} // end of namespace



  /// Hat hier nichts zu suchen
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

    double scaling_space = 1.0;
    double scaling_time = 1.0;
    
    Vec<D> offset_space = Vec<D>(0.);
    double offset_time = 0.;

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
  };

  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;

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
