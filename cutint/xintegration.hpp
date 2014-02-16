#ifndef FILE_XINTEGRATION_HPP
#define FILE_XINTEGRATION_HPP

#include "../common/spacetimefe.hpp"   // for ScalarSpaceTimeFiniteElement
#include "xdecompose.hpp"

#include <set>
#include <vector>

using namespace ngfem;

namespace xintegration
{
  /// struct which defines the relation a < b for Point4DCL 
  /// (in order to use std::set-features)
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


  /// Container set constitutes a collection of Vec<D> 
  /// main feature: the operator()(const PointXDCL & p)
  /// The points in the container are owned and later 
  /// released by PointContainer
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
    PointContainer();
    
    /// Access operator to points
    /// Either point is already in the Container, 
    ///   then return pointer to that point
    /// or point is not in the Container yet,
    ///   then add point to container and return pointer to new Vec<D>
    /// The return value (pointer ) points to a Vec<D> which is owned
    /// and later released by PointContainer
    const Vec<SD>* operator()(const Vec<SD> & p);

    void Report(std::ostream & out) const;

    ~PointContainer(){};
  };

  /// outstream which add the identifier for the domain types
  ostream & operator<<(ostream & s, DOMAIN_TYPE dt);

  /// simple class that constitutes a quadrature
  /// (here we don't use) IntegrationPoints or IntegrationRule
  /// because of a potential 4th dimension
  template < int SD >
  struct QuadratureRule
  {
    /// the quadrature points
    Array < Vec<SD> > points;
    /// the quadrature weights
    Array < double > weights;
    /// return number of integration points 
    int Size() { return points.Size(); }
  };

  /// simple class that constitutes a quadrature
  /// (here we don't use) IntegrationPoints or IntegrationRule
  /// because of a potential 4th dimension
  template < int SD >
  struct QuadratureRuleCoDim1
  {
    /// the quadrature points
    Array < Vec<SD> > points;
    /// the quadrature weights
    Array < double > weights;
    /// the quadrature weights
    Array < Vec<SD> > normals;
    /// return number of integration points 
    int Size() { return points.Size(); }
  };

  /// class that constitutes the components of a composite 
  /// quadrature rule 
  template < int SD >
  struct CompositeQuadratureRule
  {
    /// Quadrature rule for positive domain
    QuadratureRule<SD> quadrule_pos;
    /// Quadrature rule for negative domain
    QuadratureRule<SD> quadrule_neg;
    /// Quadrature rule for interface
    QuadratureRuleCoDim1<SD> quadrule_if;
    /// Interface for accessing the rules via DOMAIN_TYPE
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
      default:
        throw Exception(" DOMAIN_TYPE not known ");
        return quadrule_neg;
      }
    }

    QuadratureRuleCoDim1<SD> & GetInterfaceRule()
    {
      return quadrule_if;
    }

  };

  class XLocalGeometryInformation
  {
      // empty 
  public:
    XLocalGeometryInformation() {;}
    ~XLocalGeometryInformation() {;}
    virtual double EvaluateLsetAtPoint( const IntegrationPoint & ip, double time = 0) const;
    virtual DOMAIN_TYPE MakeQuadRule() const ;
    // static XLocalGeometryInformation * Create(ELEMENT_TYPE ET_SPACE,
    //                                           ELEMENT_TYPE ET_TIME,
    //                                           const ScalarFEEvaluator<D> & a_lset, 
    //                                           CompositeQuadratureRule<SD> & a_compquadrule,
    //                                           LocalHeap & a_lh,
    //                                           int a_int_order_space, int a_int_order_time, 
    //                                           int a_ref_level_space, int a_ref_level_time);
  };

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  class NumericalIntegrationStrategy : public XLocalGeometryInformation
  {
  public:

    /// Space dimension 
    enum { D = ET_trait<ET_SPACE>::DIM };
    /// Space-time dimension (if not space time SD==D)
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM };

    /// Levelset function through the evaluator
    // const ScalarFEEvaluator<D> & lset;
    const ScalarFieldEvaluator & lset;

    virtual double EvaluateLsetAtPoint( const IntegrationPoint & ip, double time = 0)
    {
        Vec<SD> p;
        for (int i = 0; i < D; ++i)
          p[i] = ip(i);
        if (ET_trait<ET_TIME>::DIM==1)
          p[SD-1] = time;
        return lset(p);
    }
      
    /// PointContainer contains all points that are used during the decomposition
    /// With this container multiple allocations are avoided
    /// Since this (adaptive) strategy only needs to store the vertices locally
    /// This Container might be dropped sooner or later
    PointContainer<SD> & pc;

    /// vertices of the spatial element
    Array< Vec<D> > verts_space;
    /// vertices of the temporal element (t0,t1,t2,t3,..,tN) 
    /// with N = 2^rn + 1 with rn = ref_level_time
    Array< double > verts_time;

    /// maximum number of refinements in space for numerical integration
    int ref_level_space = 0;
    /// maximum number of refinements in time for numerical integration
    int ref_level_time = 0;
    /// integration order in space (in case a tensor-product rule can be applied)
    int int_order_space = 0;
    /// integration order in time (in case a tensor-product rule can be applied)
    int int_order_time = 0;

    /// once a level absolute value is larger than threshold the prism is considered non-intersected
    double distance_threshold = 1e99;

    void SetDistanceThreshold( const double & a_distance_threshold ){ distance_threshold = a_distance_threshold; }

    LocalHeap & lh;
    /// 
    CompositeQuadratureRule<SD> & compquadrule;

    /// top level
    bool ownpc = false;

    /// Integration Order which is used on decomposed geometries
    /// it's the maximum of int_order_space and int_order_time
    int GetIntegrationOrderMax() const
    { 
      return max(int_order_space, int_order_time); 
    }
    
    /// constructor: basically copying from input, typically adapting ref_level 
    /// in space or time. Typically called from an adaptive strategy.
    NumericalIntegrationStrategy(const NumericalIntegrationStrategy & a, 
                                 int reduce_ref_space = 0, 
                                 int reduce_ref_time = 0);

    /// constructor: prescribing all the input except for the vertices (space and time)
    NumericalIntegrationStrategy(const ScalarFieldEvaluator & a_lset, 
                                 PointContainer<SD> & a_pc,
                                 CompositeQuadratureRule<SD> & a_compquadrule,
                                 LocalHeap & lh,
                                 int a_int_order_space = 2, 
                                 int a_int_order_time = 2, 
                                 int a_ref_level_space = 0, 
                                 int a_ref_level_time = 0 );

    /// constructor: prescribing all the input except for the vertices (space and time)
    NumericalIntegrationStrategy(const ScalarFieldEvaluator & a_lset, 
                                 CompositeQuadratureRule<SD> & a_compquadrule,
                                 LocalHeap & lh,
                                 int a_int_order_space = 2, 
                                 int a_int_order_time = 2, 
                                 int a_ref_level_space = 0, 
                                 int a_ref_level_time = 0 );

    /// Set Vertices according to input
    void SetVerticesSpace(const Array<Vec<D> > & verts);

    /// Set Vertices according reference geometry of ET_SPACE
    void SetVerticesSpace();

    /// Set Vertices according to input
    void SetVerticesTime(const Array<double> & verts);

    /// Set Vertices according int_order_time
    void SetVerticesTime();

    /// Set Vertices according to latest half of input
    void SetVerticesTimeFromUpperHalf(const Array< double >& verts_t);

    /// Set Vertices according to first half of input
    void SetVerticesTimeFromLowerHalf(const Array< double >& verts_t);

    /// Check prism for cut
    DOMAIN_TYPE CheckIfCut() const;

    /// Call adaptive strategy to generate quadrature rule
    /// adaptive strategy to generate composite quadrature rule on tensor product geometry
    /// ...
    virtual DOMAIN_TYPE MakeQuadRule() const;
    
  };


  /// Different Decomposition for different dimensions 

  // use a standard simplex rule to fill uncut simplices
  template <int D>
  void FillSimplexWithRule (const Simplex<D> & s, QuadratureRule<D> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<1> (const Simplex<1> & s, QuadratureRule<1> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<2> (const Simplex<2> & s, QuadratureRule<2> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<3> (const Simplex<3> & s, QuadratureRule<3> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<4> (const Simplex<4> & s, QuadratureRule<4> & quaddom, int intorder);

  template <int D>
  void FillSimplexCoDim1WithRule (const Array< const Vec<D> *> & s, const Vec<D> & pospoint,
                                  QuadratureRuleCoDim1<D> & quaddom, int intorder);
  template <>
  void FillSimplexCoDim1WithRule<2> (const Array< const Vec<2> *> & s, const Vec<2> & pospoint,
                                     QuadratureRuleCoDim1<2> & quaddom, int intorder);
  template <>
  void FillSimplexCoDim1WithRule<3> (const Array< const Vec<3> *> & s, const Vec<3> & pospoint,
                                     QuadratureRuleCoDim1<3> & quaddom, int intorder);
  template <>
  void FillSimplexCoDim1WithRule<4> (const Array< const Vec<4> *> & s, const Vec<4> & pospoint,
                                     QuadratureRuleCoDim1<4> & quaddom, int intorder);

  namespace DecompositionRules
  {
    template<int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    struct CutSimplex
    {
        static void MakeQuad(const Simplex <D> & s, 
                             const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint);
    };

    //partial specialization
    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    struct CutSimplex<1,ET_SPACE,ET_TIME>
    {
        static void MakeQuad(const Simplex <1> & s, 
                             const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint);
    };

    //partial specialization
    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    struct CutSimplex<2,ET_SPACE,ET_TIME>
    {
        static void MakeQuad(const Simplex <2> & s, 
                             const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint);
    };

    //partial specialization
    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    struct CutSimplex<3,ET_SPACE,ET_TIME>
    {
        static void MakeQuad(const Simplex <3> & s, 
                             const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint);
    };

  }

  template <int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void MakeQuadRuleOnCutSimplex(const Simplex <D> & s, 
                                const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint);


}

#endif
