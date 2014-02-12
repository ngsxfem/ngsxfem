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
    ///   then add point to container and return pointer to new Vec<D>
    /// The return value (pointer ) points to a Vec<D> which is owned
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
    QuadratureRule<SD> quadrule_if;
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
      case IF: 
        return quadrule_if;
        break;
      default:
        throw Exception(" DOMAIN_TYPE not known ");
        return quadrule_if;
      }
    }
  };

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

  template <int D>
  void FillSimplexWithRule (const Simplex<D> & s, QuadratureRule<D> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<2> (const Simplex<2> & s, QuadratureRule<2> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<3> (const Simplex<3> & s, QuadratureRule<3> & quaddom, int intorder);
  template <>
  void FillSimplexWithRule<4> (const Simplex<4> & s, QuadratureRule<4> & quaddom, int intorder);



  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  class NumericalIntegrationStrategy
  {
  public:

    /// Space dimension 
    enum { D = ET_trait<ET_SPACE>::DIM };
    /// Space-time dimension (if not space time SD==D)
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM };
    /// Levelset function through the evaluator
    const ScalarSpaceTimeFEEvaluator<D> & lset;
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

    LocalHeap & lh;
    /// 
    CompositeQuadratureRule<SD> & compquadrule;

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
    NumericalIntegrationStrategy(const ScalarSpaceTimeFEEvaluator<D> & a_lset, 
                                 PointContainer<SD> & a_pc, 
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

  };


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE CheckIfCut(const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {

    bool haspos = false;
    bool hasneg = false;

    int np1ds = pow(2,numint.ref_level_space);
    int np1dt = pow(2,numint.ref_level_time);

    double dx_scalar = 1.0 / np1ds;
    
    Array<double> time (0);
    // double ht = 1.0;

    switch (ET_SPACE)
    {
    case ET_TRIG:
    case ET_TET:
    {
      // int sum = 0;
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
    default:
      throw Exception(" this CheckIfCut::ELEMENT_TYPE is not treated... yet ");
      break;
    }
  }


  /// adaptive strategy to generate composite quadrature rule on tensor product geometry
  /// ...
  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE MakeQuadRule(const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {
    enum { D = ET_trait<ET_SPACE>::DIM };
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM};

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
        // DOMAIN_TYPE dt_upper = 
        MakeQuadRule(numint_upper);

        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_lower (numint, 0, 1);
        numint_lower.SetVerticesSpace(numint.verts_space);
        numint_lower.SetVerticesTimeFromLowerHalf(numint.verts_time);
        // cout << " call lower " << endl;
        // DOMAIN_TYPE dt_lower = 
        MakeQuadRule(numint_lower);

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
            // DOMAIN_TYPE dt_i = 
            MakeQuadRule(numint_i);
            // cout << " dt_s(" << i << "): " << dt_i << endl;

          }
          
        }
        else
          throw Exception(" tut tut - 3D not yet implemented");
      } 
      else
      {
        Array<Simplex<SD> *> simplices;
        const int nvt = ET_TIME == ET_SEGM ? 2 : 1;
        const int nvs = numint.verts_space.Size();
        Array<Vec<SD> > verts(nvs * nvt);

        for (int K = 0; K < 2; ++K)
          for (int i = 0; i < nvs; ++i)
          {
            for (int j = 0; j < D; ++j)
              verts[i+K*nvs][j] = numint.verts_space[i][j]; 
            verts[i+K*nvs][D] = K == 0? numint.verts_time[0] : numint.verts_time[numint.verts_time.Size()-1];
          }
        
        DecomposeIntoSimplices<ET_SPACE,ET_TIME>(verts, simplices, numint.pc, numint.lh);

        const ngfem::ScalarSpaceTimeFEEvaluator<D> & eval (numint.lset);

        for (int i = 0; i < simplices.Size(); ++i)
        {
          DOMAIN_TYPE dt_simplex = simplices[i]->CheckIfCut(eval);

          if (dt_simplex == IF)
          {
            if (SD==3 && ET_TIME == ET_SEGM)
            {
              Array< Vec<SD> > cutpoints(0);
              Array< Vec<SD> > pospoints(0);
              Array< Vec<SD> > negpoints(0);
              
              const int edge[6][2] = { {0, 1},
                                       {0, 2},
                                       {0, 3},
                                       {1, 2},
                                       {1, 3},
                                       {2, 3}};

              double vvals[4];
              for (int j = 0; j < 4; ++j)
              {
                vvals[j] = numint.lset(*simplices[i]->p[j]);
                if (vvals[j] >= 0)
                  pospoints.Append(*simplices[i]->p[j]);
                else
                  negpoints.Append(*simplices[i]->p[j]);
              }
              for (int j = 0; j < 6; ++j)
              {
                const int lv = edge[j][0];
                const int rv = edge[j][0];
                const double valleft = vvals[lv];
                const double valright = vvals[rv];
                if (valleft * valright < 0)
                {
                  const double cutpos = valleft / (valleft - valright);
                  Vec<SD> p = (1-cutpos) * (*simplices[i]->p[lv]) + cutpos * (*simplices[i]->p[rv]) ;
                  cutpoints.Append(p);
                  negpoints.Append(p);
                  pospoints.Append(p);
                }
              }
              //TODO continue here...
            }
            else
              throw Exception(" only space-time in 2D so far");
          }
          else
          {
            FillSimplexWithRule<SD>(*simplices[i], 
                                    numint.compquadrule.GetRule(dt_simplex), 
                                    numint.GetIntegrationOrderMax());
          }
        }
        // cout << "TODO: Call Decomposition into uncut simplices" << endl;

      }
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
            point += ir_space[k](m) * numint.verts_space[m+1];

          const double weight = ir_time[l].Weight() * dt * ir_space[k].Weight() * trafofac;
          Vec<ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM> ipoint(0.0);
          
          for (int m = 0; m < D ;++m)
            ipoint(m) = point(m);
          if (ET_trait<ET_TIME>::DIM > 0)
            ipoint(SD-1) = current;
          
          numint.compquadrule.GetRule(dt_self).points.Append(ipoint);
          numint.compquadrule.GetRule(dt_self).weights.Append(weight);

        }
      }
      // fill integration rule already ? 
      return dt_self;
    }
  }

}

#endif
