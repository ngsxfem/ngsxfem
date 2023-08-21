#ifndef FILE_XINTEGRATION_HPP
#define FILE_XINTEGRATION_HPP

// #include "../spacetime/spacetimefe.hpp"   // for ScalarSpaceTimeFiniteElement
#include "xdecompose.hpp"
#include "lsetintdomain.hpp"

#include <set>
#include <vector>

using namespace ngfem;
using ngfem::ELEMENT_TYPE;
namespace xintegration
{
  //TODO: switch to FlatArray (call be reference or remove tuple for that)
  tuple<const IntegrationRule *, Array<double>> CreateCutIntegrationRule(const LevelsetIntegrationDomain & lsetintdom,
                                                                         const ElementTransformation & trafo,
                                                                         LocalHeap & lh);
  template<class SCAL> 
  /**
   * Creates a new FlatArray of SIMD values from a given FlatArray of scalar values.
   *
   * @param ns_arr The FlatArray of scalar values to convert to SIMD values.
   * @param lh The LocalHeap to allocate memory from for the new FlatArray.
   * @return A new FlatArray of SIMD values.
   */
  FlatArray<SIMD<SCAL>> CreateSIMD_FlatArray(FlatArray<SCAL> ns_arr, LocalHeap & lh);

  /// OLD STYLE (to be removed on the long run, hopefully)
  /// struct which defines the relation a < b for Point4DCL 
  tuple<const IntegrationRule *, Array<double> > CreateCutIntegrationRule(shared_ptr<CoefficientFunction> cflset,
                                                   shared_ptr<GridFunction> gflset,
                                                   const ElementTransformation & trafo,
                                                   DOMAIN_TYPE dt,
                                                   int intorder,
                                                   int time_intorder,
                                                   LocalHeap & lh,
                                                   int subdivlvl = 0,
                                                   SWAP_DIMENSIONS_POLICY quad_dir_policy = FIND_OPTIMAL);

  
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
  //ostream & operator<<(ostream & s, DOMAIN_TYPE dt);

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
    int Size() const { return points.Size(); }
  };

  /// simple class that constitutes a quadrature rule
  /// with Flat memory layout
  template < int SD >
  class FlatQuadratureRule
  {
  public:
    /// the quadrature points
    FlatMatrixFixWidth<SD> points;
    /// the quadrature weights
    FlatVector<double> weights;
    /// return number of integration points 
    int Size() const { return points.Height(); }
    FlatQuadratureRule( const QuadratureRule<SD> & orig, LocalHeap & lh)
      : points(orig.Size(),lh), weights(orig.Size(),lh)
    {
      for (int i = 0; i < orig.Size(); ++i)
      {
        points.Row(i) = orig.points[i];
        weights(i) = orig.weights[i];
      }
    }
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
    /// the quadrature normal vectors
    Array < Vec<SD> > normals;
    /// return number of integration points 
    int Size() const { return points.Size(); }
  };

  /// simple class that constitutes a quadrature rule
  /// with Flat memory layout
  template < int SD >
  class FlatQuadratureRuleCoDim1
  {
  public:
    /// the quadrature points
    FlatMatrixFixWidth<SD> points;
    /// the quadrature weights
    FlatVector<double> weights;
    /// the quadrature normal vectors
    FlatMatrixFixWidth<SD> normals;
    /// return number of integration points 
    int Size() const { return points.Height(); }
    FlatQuadratureRuleCoDim1( const QuadratureRuleCoDim1<SD> & orig, LocalHeap & lh)
      : points(orig.Size(),lh), weights(orig.Size(),lh), normals(orig.Size(),lh)
    {
      for (int i = 0; i < orig.Size(); ++i)
      {
        points.Row(i) = orig.points[i];
        weights(i) = orig.weights[i];
        normals.Row(i) = orig.normals[i];
      }
    }
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

    const QuadratureRule<SD> & GetRule(DOMAIN_TYPE dt) const
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

    const QuadratureRuleCoDim1<SD> & GetInterfaceRule() const
    {
      return quadrule_if;
    }

  };

  /// class that constitutes the components of a composite 
  /// quadrature rule 
  template < int SD >
  class FlatCompositeQuadratureRule
  {
  public:
    /// Quadrature rule for positive domain
    FlatQuadratureRule<SD> quadrule_pos;
    /// Quadrature rule for negative domain
    FlatQuadratureRule<SD> quadrule_neg;
    /// Quadrature rule for interface
    FlatQuadratureRuleCoDim1<SD> quadrule_if;

    FlatCompositeQuadratureRule( const CompositeQuadratureRule<SD> & orig, LocalHeap &lh)
      : quadrule_pos(orig.quadrule_pos,lh),
        quadrule_neg(orig.quadrule_neg,lh),
        quadrule_if(orig.quadrule_if,lh)
    {
      ;
    }

    /// Interface for accessing the rules via DOMAIN_TYPE
    FlatQuadratureRule<SD> & GetRule(DOMAIN_TYPE dt)
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

    const FlatQuadratureRule<SD> & GetRule(DOMAIN_TYPE dt) const
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

    FlatQuadratureRuleCoDim1<SD> & GetInterfaceRule()
    {
      return quadrule_if;
    }

    const FlatQuadratureRuleCoDim1<SD> & GetInterfaceRule() const
    {
      return quadrule_if;
    }

  };


  template class FlatQuadratureRule<2>;
  template class FlatQuadratureRule<3>;
  template class FlatQuadratureRule<4>;
  template class FlatQuadratureRuleCoDim1<2>;
  template class FlatQuadratureRuleCoDim1<3>;
  template class FlatQuadratureRuleCoDim1<4>;
  template class FlatCompositeQuadratureRule<2>;
  template class FlatCompositeQuadratureRule<3>;
  template class FlatCompositeQuadratureRule<4>;


  class XLocalGeometryInformation
  {
  protected:
    XLocalGeometryInformation * pasttracegeom = 0;
    XLocalGeometryInformation * futuretracegeom = 0;
    mutable bool quaded = false;
    // empty 
  public:

    /// Levelset function through the evaluator
    // const ScalarFEEvaluator<D> & lset;
    const ScalarFieldEvaluator * lset;

    XLocalGeometryInformation(const ScalarFieldEvaluator * a_lset): lset(a_lset) {;}
    virtual ~XLocalGeometryInformation() {;}
    virtual double EvaluateLsetAtPoint( const IntegrationPoint & ip, double time = 0) const;
    virtual DOMAIN_TYPE MakeQuadRule() const ;

    virtual int Dimension()const { return -1; }

    virtual const CompositeQuadratureRule<1> * GetRule1() const { return NULL; }
    virtual const CompositeQuadratureRule<2> * GetRule2() const { return NULL; }
    virtual const CompositeQuadratureRule<3> * GetRule3() const { return NULL; }
    virtual const CompositeQuadratureRule<4> * GetRule4() const { return NULL; }

    template <int SD>
    const CompositeQuadratureRule<SD> * GetCompositeRule() const
    {
      switch(SD)
      {
      case 1:
        return reinterpret_cast<const CompositeQuadratureRule<SD>*>(GetRule1()); break;
      case 2:
        return reinterpret_cast<const CompositeQuadratureRule<SD>*>(GetRule2()); break;
      case 3:
        return reinterpret_cast<const CompositeQuadratureRule<SD>*>(GetRule3()); break;
      case 4:
        return reinterpret_cast<const CompositeQuadratureRule<SD>*>(GetRule4()); break;
      default:
        return NULL;
        break;
      }
    }

    // template <int SD>
    static shared_ptr<XLocalGeometryInformation> Create(ELEMENT_TYPE ET_SPACE,
                                                        ELEMENT_TYPE ET_TIME,
                                                        const ScalarFieldEvaluator & a_lset, 
                                                        CompositeQuadratureRule<1> & a_compquadrule,
                                                        LocalHeap & a_lh,
                                                        int a_int_order_space = 1, int a_int_order_time = 1, 
                                                        int a_ref_level_space = 0, int a_ref_level_time = 0);
    static shared_ptr<XLocalGeometryInformation> Create(ELEMENT_TYPE ET_SPACE,
                                                        ELEMENT_TYPE ET_TIME,
                                                        const ScalarFieldEvaluator & a_lset, 
                                                        CompositeQuadratureRule<2> & a_compquadrule,
                                                        LocalHeap & a_lh,
                                                        int a_int_order_space = 1, int a_int_order_time = 1, 
                                                        int a_ref_level_space = 0, int a_ref_level_time = 0);
    static shared_ptr<XLocalGeometryInformation>  Create(ELEMENT_TYPE ET_SPACE,
                                                         ELEMENT_TYPE ET_TIME,
                                                         const ScalarFieldEvaluator & a_lset, 
                                                         CompositeQuadratureRule<3> & a_compquadrule,
                                                         LocalHeap & a_lh,
                                                         int a_int_order_space = 1, int a_int_order_time = 1, 
                                                         int a_ref_level_space = 0, int a_ref_level_time = 0);
    static shared_ptr<XLocalGeometryInformation> Create(ELEMENT_TYPE ET_SPACE,
                                                        ELEMENT_TYPE ET_TIME,
                                                        const ScalarFieldEvaluator & a_lset, 
                                                        CompositeQuadratureRule<4> & a_compquadrule,
                                                        LocalHeap & a_lh,
                                                        int a_int_order_space = 1, int a_int_order_time = 1, 
                                                        int a_ref_level_space = 0, int a_ref_level_time = 0);
    static shared_ptr<XLocalGeometryInformation>  Create(ELEMENT_TYPE ET_SPACE,
                                                         ELEMENT_TYPE ET_TIME,
                                                         const ScalarFieldEvaluator & a_lset, 
                                                         CompositeQuadratureRule<1> * a_compquadrule1,
                                                         CompositeQuadratureRule<2> * a_compquadrule2,
                                                         CompositeQuadratureRule<3> * a_compquadrule3,
                                                         CompositeQuadratureRule<4> * a_compquadrule4,
                                                         LocalHeap & a_lh,
                                                         int a_int_order_space = 1, int a_int_order_time = 1, 
                                                         int a_ref_level_space = 0, int a_ref_level_time = 0);
    bool IsDecomposed(){ return quaded; }
      
    void SetPastTrace(XLocalGeometryInformation * xlocal)
    {
      pasttracegeom = xlocal;
    }

    void SetFutureTrace(XLocalGeometryInformation * xlocal)
    {
      futuretracegeom = xlocal;
    }

    XLocalGeometryInformation * GetPastTrace() const 
    {
      if (pasttracegeom == NULL)
        throw Exception("XLocalGeometryInformation::GetPastTrace() called but pasttracegeom == NULL");
      if (! pasttracegeom->IsDecomposed())
        pasttracegeom->MakeQuadRule();
      return pasttracegeom;
    }

    XLocalGeometryInformation * GetFutureTrace() const
    {
      if (futuretracegeom == NULL)
        throw Exception("XLocalGeometryInformation::GetFutureTrace() called but futuretracegeom == NULL");
      if (! futuretracegeom->IsDecomposed())
        futuretracegeom->MakeQuadRule();
      return futuretracegeom;
    }
    
    virtual void SetDistanceThreshold( double a_distance_threshold )
    {  
      std::cout << IM(3) << " base class is doing nothing " << std::endl;
    }


    virtual void SetSimplexArrays(Array<Simplex<2>*> & simplex_array_neg,
                                  Array<Simplex<2>*> & simplex_array_pos)
    { cout << IM(3) << " baseclass: doing nothing" << endl;} 
    virtual void SetSimplexArrays(Array<Simplex<3>*> & simplex_array_neg,
                                  Array<Simplex<3>*> & simplex_array_pos)
    { cout << IM(3) << " baseclass: doing nothing" << endl;} 
    virtual void SetSimplexArrays(Array<Simplex<4>*> & simplex_array_neg,
                                  Array<Simplex<4>*> & simplex_array_pos)
    { cout << IM(3) << " baseclass: doing nothing" << endl;} 
    virtual void ClearArrays(){ cout << IM(3) << " baseclass: doing nothing" << endl;}


  };

  class FlatXLocalGeometryInformation
  {
  protected:
    FlatXLocalGeometryInformation * pasttracegeom = 0;
    FlatXLocalGeometryInformation * futuretracegeom = 0;
  public:
    double kappa[2]; //Sum of weights of the POS/NEG part divided by total sum of weights

    const ScalarFieldEvaluator * lset;

    FlatCompositeQuadratureRule<1> * compquadrule1 = 0;
    FlatCompositeQuadratureRule<2> * compquadrule2 = 0;
    FlatCompositeQuadratureRule<3> * compquadrule3 = 0;
    FlatCompositeQuadratureRule<4> * compquadrule4 = 0;

    int Dimension = -1;

    bool empty;

    FlatXLocalGeometryInformation()
      : lset(NULL), Dimension(-1), empty(true)
    {
      ;
    }

    FlatXLocalGeometryInformation(const XLocalGeometryInformation & xgeom, LocalHeap & lh) 
      : lset(xgeom.lset), Dimension(xgeom.Dimension()), empty(false)
        // :
    {
      kappa[NEG] = kappa[POS] = 0.0;
      switch(Dimension)
      {
      case 1:
        compquadrule1 = new (lh) FlatCompositeQuadratureRule<1>(*xgeom.GetCompositeRule<1>(),lh);
        for (DOMAIN_TYPE dt : {POS, NEG})
        {
          FlatQuadratureRule<1> qr (compquadrule1->GetRule(dt));
          for (int i = 0; i < qr.Size(); ++i)
            kappa[dt] += qr.weights(i);
        }
        break;
      case 2:
        compquadrule2 = new (lh) FlatCompositeQuadratureRule<2>(*xgeom.GetCompositeRule<2>(),lh);
        for (DOMAIN_TYPE dt : {POS, NEG})
        {
          FlatQuadratureRule<2> qr (compquadrule2->GetRule(dt));
          for (int i = 0; i < qr.Size(); ++i)
            kappa[dt] += qr.weights(i);
        }
        break;
      case 3:
        compquadrule3 = new (lh) FlatCompositeQuadratureRule<3>(*xgeom.GetCompositeRule<3>(),lh);
        for (DOMAIN_TYPE dt : {POS, NEG})
        {
          FlatQuadratureRule<3> qr (compquadrule3->GetRule(dt));
          for (int i = 0; i < qr.Size(); ++i)
            kappa[dt] += qr.weights(i);
        }
        break;
      case 4:
        compquadrule4 = new (lh) FlatCompositeQuadratureRule<4>(*xgeom.GetCompositeRule<4>(),lh);
        for (DOMAIN_TYPE dt : {POS, NEG})
        {
          FlatQuadratureRule<4> qr (compquadrule4->GetRule(dt));
          for (int i = 0; i < qr.Size(); ++i)
            kappa[dt] += qr.weights(i);
        }
        break;
      default:
        throw Exception("Dimension not in {1,2,3,4}");
        break;
      }
      
      const double sum = kappa[NEG] + kappa[POS];
      kappa[NEG] /= sum;
      kappa[POS] /= sum;
      ;
    }
    virtual ~FlatXLocalGeometryInformation() {;}

    template <int SD>
    const FlatCompositeQuadratureRule<SD> & GetCompositeRule() const
    {
      switch(SD)
      {
      case 1:
        return reinterpret_cast<FlatCompositeQuadratureRule<SD>&>(*compquadrule1); break;
      case 2:
        return reinterpret_cast<FlatCompositeQuadratureRule<SD>&>(*compquadrule2); break;
      case 3:
        return reinterpret_cast<FlatCompositeQuadratureRule<SD>&>(*compquadrule3); break;
      case 4:
        return reinterpret_cast<FlatCompositeQuadratureRule<SD>&>(*compquadrule4); break;
      default:
        throw Exception(" SD not in {1,2,3,4} ");
        return reinterpret_cast<FlatCompositeQuadratureRule<SD>&>(*compquadrule1); break;
        break;
      }
    }

    template <int D, int SD>
    double EvaluateLsetAtPoint( const IntegrationPoint & ip, double time = 0) const
    {
        Vec<SD> p;
        for (int i = 0; i < D; ++i)
          p[i] = ip(i);
        if (SD>D)
          p[SD-1] = time;
        return (*lset)(p);
    }

    void SetPastTrace(FlatXLocalGeometryInformation * xlocal)
    {
      pasttracegeom = xlocal;
    }

    void SetFutureTrace(FlatXLocalGeometryInformation * xlocal)
    {
      futuretracegeom = xlocal;
    }

    FlatXLocalGeometryInformation * GetPastTrace() const 
    {
      if (pasttracegeom == NULL)
        throw Exception("FlatXLocalGeometryInformation::GetPastTrace() called but pasttracegeom == NULL");
      // if (! pasttracegeom->IsDecomposed())
      //     pasttracegeom->MakeQuadRule();
      return pasttracegeom;
    }

    FlatXLocalGeometryInformation * GetFutureTrace() const
    {
      if (futuretracegeom == NULL)
        throw Exception("FlatXLocalGeometryInformation::GetFutureTrace() called but futuretracegeom == NULL");
      // if (! futuretracegeom->IsDecomposed())
      //     futuretracegeom->MakeQuadRule();
      return futuretracegeom;
    }
  };


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  class NumericalIntegrationStrategy : public XLocalGeometryInformation
  {
  public:

    /// Space dimension 
    enum { D = ET_trait<ET_SPACE>::DIM };
    /// Space-time dimension (if not space time SD==D)
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM };

    virtual int Dimension() const{ return SD; }

    virtual double EvaluateLsetAtPoint( const IntegrationPoint & ip, double time = 0) const
    {
        Vec<SD> p;
        for (int i = 0; i < D; ++i)
          p[i] = ip(i);
        if (ET_trait<ET_TIME>::DIM==1)
          p[SD-1] = time;
        return (*lset)(p);
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

    Array< Simplex<SD> *> * simplex_array_neg = NULL;
    Array< Simplex<SD> *> * simplex_array_pos = NULL;

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

    virtual void SetDistanceThreshold( double a_distance_threshold ){ distance_threshold = a_distance_threshold; }

    virtual void SetSimplexArrays(Array<Simplex<2> *> & a_simplex_array_neg, 
                                 Array<Simplex<2> *> & a_simplex_array_pos)
    { 
      if (2==SD)
      {
        simplex_array_neg = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_neg);
        simplex_array_pos = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_pos);
      }
      else
        throw Exception("Dimensions do not match 1337!");
    } 

    virtual void SetSimplexArrays(Array<Simplex<3> *> & a_simplex_array_neg, 
                                 Array<Simplex<3> *> & a_simplex_array_pos)
    { 
      throw Exception("3D not yet checked...");
      if (3==SD)
      {
        simplex_array_neg = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_neg);
        simplex_array_pos = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_pos);
      }
      else
        throw Exception("Dimensions do not match 1337!");
    } 

    virtual void SetSimplexArrays(Array<Simplex<4> *> & a_simplex_array_neg, 
                                 Array<Simplex<4> *> & a_simplex_array_pos)
    { 
      throw Exception("4D not yet checked...");
      if (4==SD)
      {
        simplex_array_neg = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_neg);
        simplex_array_pos = reinterpret_cast<Array<Simplex<SD>*>*> (&a_simplex_array_pos);
      }
      else
        throw Exception("Dimensions do not match 1337!");
    } 

    virtual void ClearArrays()
    {
      if (simplex_array_neg != NULL)
      {
        for (int i = 0; i < simplex_array_neg->Size(); ++i)
          delete (*simplex_array_neg)[i];
        simplex_array_neg->SetSize(0);
      }
      if (simplex_array_pos != NULL)
      {
        for (int i = 0; i < simplex_array_pos->Size(); ++i)
          delete (*simplex_array_pos)[i];
        simplex_array_pos->SetSize(0);
      }
      simplex_array_neg = NULL;
      simplex_array_pos = NULL;
    }


    LocalHeap & lh;
    /// 
    CompositeQuadratureRule<SD> & compquadrule;

    virtual const CompositeQuadratureRule<1> * GetRule1() const  { 
      if (SD==1) return reinterpret_cast<CompositeQuadratureRule<1>*>(&compquadrule); 
      else return NULL; 
    }
    virtual const CompositeQuadratureRule<2> * GetRule2() const { 
      if (SD==2) return reinterpret_cast<CompositeQuadratureRule<2>*>(&compquadrule); 
      else return NULL; 
    }
    virtual const CompositeQuadratureRule<3> * GetRule3() const { 
      if (SD==3) return reinterpret_cast<CompositeQuadratureRule<3>*>(&compquadrule); 
      else return NULL; 
    }
    virtual const CompositeQuadratureRule<4> * GetRule4() const { 
      if (SD==4) return reinterpret_cast<CompositeQuadratureRule<4>*>(&compquadrule); 
      else return NULL; 
    }

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
    
    virtual ~NumericalIntegrationStrategy() 
    { 
      if (ownpc) delete &pc; 
    }

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


  const IntegrationRule * CutIntegrationRule(shared_ptr<CoefficientFunction> cf_lset,
                                             const ElementTransformation & trafo,
                                             DOMAIN_TYPE dt,
                                             int intorder,
                                             int subdivlvl,
                                             LocalHeap & lh);
  
}


#endif
