#include "xintegration.hpp"

namespace xintegration
{
  

  template<int SD>
  PointContainer<SD>::PointContainer()
  {
#ifdef DEBUG
    k=0;
#endif
    pset.clear();
  };

  template<int SD>
  const Vec<SD>* PointContainer<SD>::operator()(const Vec<SD> & p)
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

  template<int SD>
  void PointContainer<SD>::Report(std::ostream & out) const
  {
    out << " PointContainer stored " << pset.size() << " points.\n";
#ifdef DEBUG
    out << " PointContainer rejected " << k << " points.\n";
#endif
  }

  template class PointContainer<2>;
  template class PointContainer<3>;
  template class PointContainer<4>;

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
    return s;
  }


  template<>
  void FillSimplexWithRule<4> (const Simplex<4> & s, QuadratureRule<4> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 4");
  }

  template<>
  void FillSimplexWithRule<3> (const Simplex<3> & s, QuadratureRule<3> & quaddom, int intorder)
  {
    const double trafofac = Measure<3,3>(s.p) * 6.0;
    const IntegrationRule & ir = SelectIntegrationRule (ET_TET, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 3 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s.p[0]);
      for (int m = 0; m < 3 ;++m)
        point += ir[k](m) * (*s.p[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }

  template<>
  void FillSimplexWithRule<2> (const Simplex<2> & s, QuadratureRule<2> & quaddom, int intorder)
  {
    const double trafofac = Measure<2,2>(s.p) * 2.0;
    const IntegrationRule & ir = SelectIntegrationRule (ET_TRIG, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<2> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 2 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s.p[0]);
      for (int m = 0; m < 2 ;++m)
        point += ir[k](m) * (*s.p[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }


  template<>
  void FillSimplexCoDim1WithRule<4> (const Array< const Vec<4> *> & s, QuadratureRuleCoDim1<4> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 4");
  }

  template<>
  void FillSimplexCoDim1WithRule<3> (const Array< const Vec<3> *> & s, QuadratureRuleCoDim1<3> & quaddom, int intorder)
  {
    Vec<3> a = *s[1] - *s[0];
    Vec<3> b = *s[2] - *s[0];
    Vec<3> c = Cross(a,b);
    const double trafofac = L2Norm(c);
    c /= trafofac;
    const IntegrationRule & ir = SelectIntegrationRule (ET_TRIG, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      double originweight = 1.0;
      for (int m = 0; m < 2 ;++m)
        originweight -= ir[k](m);
      point = originweight * (*s[0]);
      for (int m = 0; m < 2 ;++m)
        point += ir[k](m) * (*s[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
      quaddom.normals.Append(c);
    }
  }

  template<>
  void FillSimplexCoDim1WithRule<2> (const Array< const Vec<2> *> & s, QuadratureRuleCoDim1<2> & quaddom, int intorder)
  {
    Vec<2> a = *s[1] - *s[0];
    Vec<2> n = (-a(1),a(0));
    const double trafofac = L2Norm(a);
    n /= trafofac;
    const IntegrationRule & ir = SelectIntegrationRule (ET_SEGM, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      point = (1.0-ir[k](0)) * (*s[0]);
      point += ir[k](0) * (*s[1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
      quaddom.normals.Append(n);
    }
  }



  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: NumericalIntegrationStrategy(const NumericalIntegrationStrategy & a, int reduce_ref_space, int reduce_ref_time)
    : lset(a.lset), pc(a.pc), 
      ref_level_space(a.ref_level_space-reduce_ref_space), ref_level_time(a.ref_level_time-reduce_ref_time),
      int_order_space(a.int_order_space), int_order_time(a.int_order_time),
      lh(a.lh), compquadrule(a.compquadrule)
  {
  }


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: NumericalIntegrationStrategy(const ScalarSpaceTimeFEEvaluator<D> & a_lset, 
                                  PointContainer<SD> & a_pc, 
                                  CompositeQuadratureRule<SD> & a_compquadrule,
                                  LocalHeap & a_lh,
                                  int a_int_order_space, int a_int_order_time, 
                                  int a_ref_level_space, int a_ref_level_time)
    : lset(a_lset), pc(a_pc), 
      ref_level_space(a_ref_level_space), ref_level_time(a_ref_level_time),
      int_order_space(a_int_order_space), int_order_time(a_int_order_time),
      lh(a_lh), compquadrule(a_compquadrule)
  {
    ;
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesSpace()
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

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesTime(const Array<double> & verts)
  {
    verts_time.SetSize(verts.Size());
    for (int i = 0; i < verts.Size(); ++i)
      verts_time[i] = verts[i];
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesTime()
  {
    const int np1dt = pow(2,ref_level_time);
    const double ht = 1.0 / np1dt;
    verts_time.SetSize(np1dt+1);
    for (int i = 0; i < np1dt + 1; ++i)
    {
      verts_time[i] = i * ht;
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesTimeFromUpperHalf(const Array< double >& verts_t)
  {
    const int newsize = (verts_t.Size()+1)/2;
    const int offset = (verts_t.Size()-1)/2;
    verts_time.SetSize(newsize);
    for (int i = 0; i < newsize; ++i)
    {
      verts_time[i] = verts_t[offset+i];
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesTimeFromLowerHalf(const Array< double >& verts_t)
  {
    const int newsize = (verts_t.Size()+1)/2;
    verts_time.SetSize(newsize);
    for (int i = 0; i < newsize; ++i)
    {
      verts_time[i] = verts_t[i];
    }
  }

  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: SetVerticesSpace(const Array<Vec<D> > & verts)
  {
    verts_space.SetSize(verts.Size());
    for (int i = 0; i < verts.Size(); ++i)
      verts_space[i] = verts[i];
  }

  // Check prism for cut
  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: CheckIfCut() const
  {

    bool haspos = false;
    bool hasneg = false;

    int np1ds = pow(2,ref_level_space);
    int np1dt = pow(2,ref_level_time);

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
            position[d] = verts_space[0][d];

          for (int j = 0; j < ET_trait<ET_SPACE>::DIM; ++j)
          {
            for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
            {
               position[d] += I[j] * dx_scalar * (verts_space[j+1][d] - verts_space[0][d]);
            }
          }
          // for (int d = 0; d < ET_trait<ET_SPACE>::DIM; ++d)
          //   cout << position[d] << ",\t";

          if (ET_TIME == ET_SEGM)
          {
            position[ET_trait<ET_SPACE>::DIM] = verts_time[i];
            // cout << position[ET_trait<ET_SPACE>::DIM] << ",\t";
          }
          const ngfem::ScalarSpaceTimeFEEvaluator<ET_trait<ET_SPACE>::DIM> & eval (lset);
          const double lsetval = eval(position);

          if (lsetval > distance_threshold)
            return POS;

          if (lsetval < -distance_threshold)
            return NEG;

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


  template <ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  DOMAIN_TYPE NumericalIntegrationStrategy<ET_SPACE,ET_TIME> 
  :: MakeQuadRule() const
  {
    enum { D = ET_trait<ET_SPACE>::DIM }; // spatial dimension
    enum { SD = ET_trait<ET_SPACE>::DIM + ET_trait<ET_TIME>::DIM}; // total dimension (space+time)

    // check with the help of regularly distributed points if current
    // space(-time) geometry is cut (has different sign in lset-value)
    DOMAIN_TYPE dt_self = CheckIfCut();

    if (dt_self == IF)
    {
      // cout << " cut " << endl;

      bool refine_time = ET_TIME == ET_SEGM && ref_level_time > 0;
      bool refine_space = ref_level_space > 0;

      // only refinement in space or time for the recursive call
      // decide which refinement should come first here:
      if (refine_time && refine_space)
      {
        if (ref_level_time >= ref_level_space)
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
      
      // divide space-time prism into upper and lower half
      if (refine_time)
      {
        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_upper (*this, 0, 1);
        numint_upper.SetVerticesSpace(verts_space);
        numint_upper.SetVerticesTimeFromUpperHalf(verts_time);
        numint_upper.MakeQuadRule();  // recursive call!

        NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_lower (*this, 0, 1);
        numint_lower.SetVerticesSpace(verts_space);
        numint_lower.SetVerticesTimeFromLowerHalf(verts_time);
        numint_lower.MakeQuadRule();  // recursive call!
      }

      // divide space-time prism into 4/8 prism of same height
      if (refine_space)
      {
        if ( ET_SPACE == ET_TRIG)
        {
          // barycentric coordinates for new points
          static double baryc[6][3] = { { 0.0, 0.0, 1.0},
                                        { 0.5, 0.0, 0.5},
                                        { 1.0, 0.0, 0.0},
                                        { 0.0, 0.5, 0.5},
                                        { 0.5, 0.5, 0.0},
                                        { 0.0, 1.0, 0.0}};

          // new triangles as connectivity information of the vertices baryc above
          static int trigs[4][3] = { { 0, 1, 3},
                                     { 1, 2, 4},
                                     { 1, 3, 4},
                                     { 3, 4, 5}};
          
          for (int i = 0; i < 4; ++i) // triangles
          {
            NumericalIntegrationStrategy<ET_SPACE,ET_TIME> numint_i (*this, 1, 0);
            numint_i.SetVerticesTime(verts_time);
            Array< Vec<D> > newverts(3);
            for (int j = 0; j < 3; ++j) //vertices
            {
              newverts[j] = Vec<D>(0.0);
              for (int d = 0; d < 3; ++d) 
                newverts[j] += baryc[trigs[i][j]][d] * verts_space[d];
            }
            numint_i.SetVerticesSpace(newverts);
            numint_i.MakeQuadRule(); // recursive call!
          }
          
        }
        else
          throw Exception(" refine_space in 3D in NumInt::MakeQuad not yet implemented");
      }

      if (!refine_space && !refine_time) // already on finest level: deal with cut situation
      {
        // Generate list of vertices corresponding to simplex/prism
        Array<Simplex<SD> *> simplices;
        const int nvt = ET_TIME == ET_SEGM ? 2 : 1;
        const int nvs = verts_space.Size();
        Array<const Vec<SD> * > verts(nvs * nvt);
        for (int K = 0; K < nvt; ++K)
          for (int i = 0; i < nvs; ++i)
          {
            Vec<SD> newpoint;
            for (int j = 0; j < D; ++j)
            {
              newpoint[j] = verts_space[i][j]; 
              if (ET_TIME == ET_SEGM)
                newpoint[D] = K == 0 ? verts_time[0] : verts_time[verts_time.Size()-1];
            }
            verts[i+K*nvs] = pc(newpoint);
            // cout << "verts["<<i+K*nvs<<"]:" << newpoint << endl;
          }

        // spacetime: prism to simplices
        // only space: set simplix
        // in both cases we have a list of SD-dimensional simplices
        if (ET_TIME==ET_POINT)
        {
          simplices.SetSize(1);
          simplices[0] = new (lh) Simplex<SD>(verts);
        }
        else
        {
          DecomposePrismIntoSimplices<SD>(verts, simplices, pc, lh);
        }

        // cout << "simplices:\n";
        // for (int i = 0; i < simplices.Size(); ++i)
        // {
        //   cout << *simplices[i] << endl;
        // }

        const ngfem::ScalarSpaceTimeFEEvaluator<D> & eval (lset);

        for (int i = 0; i < simplices.Size(); ++i)
        {
          // Check for each simplex if it is cut. 
          // If yes call decomposition strategy for according dimension
          // If no  direction fill the composition rule accordingly
          DOMAIN_TYPE dt_simplex = simplices[i]->CheckIfCut(eval);
          if (dt_simplex == IF)
          {
            MakeQuadRuleOnCutSimplex<SD>(*simplices[i], *this);
          }
          else
          {
            FillSimplexWithRule<SD>(*simplices[i], 
                                    compquadrule.GetRule(dt_simplex), 
                                    GetIntegrationOrderMax());
          }
        }
      }
      return IF;
    }
    else // no cut
    {
      double trafofac = 1.0; 

      if (D==2)
      {
        Vec<D> a = verts_space[1] - verts_space[0];
        Vec<D> b = verts_space[2] - verts_space[0];
        trafofac = abs(a(0) * b(1) - a(1) * b(0));
      }
      else
      {
        Vec<D> a = verts_space[1] - verts_space[0];
        Vec<D> b = verts_space[2] - verts_space[0];
        Vec<D> c = verts_space[3] - verts_space[0];
        trafofac = abs(Determinant(a,b,c));
      }

      const double dt = verts_time[verts_time.Size()-1] - verts_time[0];
      const double t0 = verts_time[0];

      const IntegrationRule & ir_time = SelectIntegrationRule (ET_TIME, int_order_time);
      const IntegrationRule & ir_space = SelectIntegrationRule (ET_SPACE, int_order_space);
      

      for (int l = 0; l < ir_time.GetNIP(); l++)
      {
        double current = t0 + ir_time[l](0) * dt;
        for (int k = 0; k < ir_space.GetNIP(); k++)
        {
          Vec<D> point(0.0);
          double originweight = 1.0;
          for (int m = 0; m < D ;++m)
            originweight -= ir_space[k](m);
          point = originweight * verts_space[0];
          for (int m = 0; m < D ;++m)
            point += ir_space[k](m) * verts_space[m+1];
          const double weight = ir_time[l].Weight() * dt * ir_space[k].Weight() * trafofac;
          Vec<SD> ipoint(0.0);
          
          for (int m = 0; m < D ;++m)
            ipoint(m) = point(m);
          if (ET_trait<ET_TIME>::DIM > 0)
            ipoint(SD-1) = current;
          
          compquadrule.GetRule(dt_self).points.Append(ipoint);
          compquadrule.GetRule(dt_self).weights.Append(weight);

        }
      }
      return dt_self;
    }
  }

  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;


  namespace DecompositionRules
  {
    template<int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <D> & s, 
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    { 
      throw Exception("CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad --- no implementation for these Element Types");
    }


    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<3,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <3> & s, 
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    { 
      enum { SD = 3};

      Array< const Vec<SD> * > cutpoints(0);
      Array< const Vec<SD> * > pospoints(0);
      Array< const Vec<SD> * > negpoints(0);

      Array<int> posvidx(0);
      Array<int> negvidx(0);

      // vertex idx connected to cut idx (just in case of 4 cut positions)
      // connectivity information of cuts
      Array<int> v2cut_1(4); 
      Array<int> v2cut_2(4); 
      v2cut_1 = -1;
      v2cut_2 = -1;

      const int edge[6][2] = { {0, 1},
                               {0, 2},
                               {0, 3},
                               {1, 2},
                               {1, 3},
                               {2, 3}};
              
      double vvals[4];
      bool zero[4];
              
      for (int j = 0; j < 4; ++j)
      {
        zero[j] = false;
        vvals[j] = numint.lset(*(s.p[j]));
        if (vvals[j] > 0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
        }
        else if (vvals[j] < 0)
        {
          negpoints.Append(numint.pc(*(s.p[j])));
          negvidx.Append(j);
        }
        else // (vvals[j] == 0.0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
          zero[j] = true;
        }
      }

      // cout << " vvals = \n";
      // for (int l = 0; l < 4; ++l)
      //   cout << l << ":" << vvals[l] << endl;

      int cntcuts = 0;
      for (int j = 0; j < 6; ++j)
      {
        const int lv = edge[j][0];
        const int rv = edge[j][1];
        const double valleft = vvals[lv];
        const double valright = vvals[rv];
        bool hascut = (valleft * valright < 0);
        if (zero[lv] && valright < 0)
          hascut = true;

        if (zero[rv] && valleft < 0)
          hascut = true;

        if (hascut)
        {
          const double cutpos = valleft / (valleft - valright);
          // std::cout << " cutpos = " << cutpos << std::endl;
          Vec<SD> p = (1-cutpos) * *(s.p[lv]) + cutpos * *(s.p[rv]) ;
          cutpoints.Append(numint.pc(p));
          // collect connectivity of cut and vertices
          if (v2cut_1[lv] == -1)
            v2cut_1[lv] = cntcuts;
          else
            v2cut_2[lv] = cntcuts;
          if (v2cut_1[rv] == -1)
            v2cut_1[rv] = cntcuts;
          else
            v2cut_2[rv] = cntcuts;
          cntcuts ++;
        }
      }

      if (cutpoints.Size() == 3) // three intersections: prism + tetra
      {

        Array< const Vec<SD> *> & minorgroup ( negpoints.Size() > pospoints.Size() ? 
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_minor = negpoints.Size() > pospoints.Size() ? POS : NEG;
        Array< const Vec<SD> *> & majorgroup ( negpoints.Size() <= pospoints.Size() ? 
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_major = negpoints.Size() <= pospoints.Size() ? POS : NEG;

        Array<int> & majvidx( negvidx.Size() > posvidx.Size() ? negvidx : posvidx);

        for (int k = 0; k < 3; ++k)
          minorgroup.Append(cutpoints[k]);
        // minorgroup is a simplex of type dt_minor
        FillSimplexWithRule<SD>(minorgroup, 
                                numint.compquadrule.GetRule(dt_minor), 
                                numint.GetIntegrationOrderMax());

        Array< Simplex<SD> * > innersimplices(0);
        for (int k = 0; k < 3; ++k)
        {
          int corresponding_cut = v2cut_1[majvidx[k]];
          majorgroup.Append(cutpoints[corresponding_cut]);
        }
        DecomposePrismIntoSimplices<SD>(majorgroup, innersimplices, numint.pc, numint.lh);
        for (int l = 0; l < innersimplices.Size(); ++l)
        {
          FillSimplexWithRule<SD>(innersimplices[l]->p, 
                                  numint.compquadrule.GetRule(dt_major), 
                                  numint.GetIntegrationOrderMax());
        }
                
        // and the interface:
        FillSimplexCoDim1WithRule<SD> ( cutpoints, numint.compquadrule.GetInterfaceRule(), 
                                        numint.GetIntegrationOrderMax());


      }
      else if (cutpoints.Size() == 4) // four intersections: prism + prism
      {
        //pos domain
        {
          Array< const Vec<SD> *> posprism(0);
          posprism.Append(pospoints[0]);
          const int idxn = posvidx[0];
          const int cut1 = v2cut_1[idxn];
          const int cut2 = v2cut_2[idxn];
          posprism.Append(cutpoints[cut1]);
          posprism.Append(cutpoints[cut2]);
          posprism.Append(pospoints[1]);
          int cut3 = -1;
          for (int l = 0; l < 4; ++l)
            if (cut1 != l && cut2 != l)
              cut3 = l;
          int cut4 = 6 - cut3 - cut2 - cut1;

          // possibly switch orientation of cut3 / cut4
          if ((v2cut_1[negvidx[0]] == cut1 && v2cut_2[negvidx[0]] == cut4)
              || (v2cut_1[negvidx[0]] == cut2 && v2cut_2[negvidx[0]] == cut3))
          {
            const int cutt = cut3;
            cut3 = cut4;
            cut4 = cutt;
          }
          posprism.Append(cutpoints[cut3]);
          posprism.Append(cutpoints[cut4]);

          // cout << " cutpoints vertices: " << endl;
          // for (int l = 0; l < 4; ++l)
          //   cout << *cutpoints[l] << endl;

          // cout << " pospoints vertices: " << endl;
          // for (int l = 0; l < 2; ++l)
          //   cout << *pospoints[l] << endl;
                  
          // cout << " posprism vertices: " << endl;
          // for (int l = 0; l < 6; ++l)
          //   cout << *posprism[l] << endl;

          Array< Simplex<SD> * > innersimplices(0);
          DecomposePrismIntoSimplices<SD>(posprism, innersimplices, numint.pc, numint.lh);
          for (int l = 0; l < innersimplices.Size(); ++l)
          {
            // std::cout << " *innersimplices[l] = " << *innersimplices[l] << std::endl;
            FillSimplexWithRule<SD>(innersimplices[l]->p, 
                                    numint.compquadrule.GetRule(POS), 
                                    numint.GetIntegrationOrderMax());
          }
        }
        //neg domain
        {
          Array< const Vec<SD> *> negprism(0);
          negprism.Append(negpoints[0]);
          const int idxn = negvidx[0];
          const int cut1 = v2cut_1[idxn];
          const int cut2 = v2cut_2[idxn];
          negprism.Append(cutpoints[cut1]);
          negprism.Append(cutpoints[cut2]);
          negprism.Append(negpoints[1]);
          int cut3 = -1;
          for (int l = 0; l < 4; ++l)
            if (cut1 != l && cut2 != l)
              cut3 = l;
          int cut4 = 6 - cut3 - cut2 - cut1;

          // possibly switch orientation of cut3 / cut4
          if ((v2cut_1[posvidx[0]] == cut1 && v2cut_2[posvidx[0]] == cut4)
              || (v2cut_1[posvidx[0]] == cut2 && v2cut_2[posvidx[0]] == cut3))
          {
            const int cutt = cut3;
            cut3 = cut4;
            cut4 = cutt;
          }
          negprism.Append(cutpoints[cut3]);
          negprism.Append(cutpoints[cut4]);

          Array< Simplex<SD> * > innersimplices(0);
          DecomposePrismIntoSimplices<SD>(negprism, innersimplices, numint.pc, numint.lh);
          for (int l = 0; l < innersimplices.Size(); ++l)
          {
            FillSimplexWithRule<SD>(innersimplices[l]->p, 
                                    numint.compquadrule.GetRule(NEG), 
                                    numint.GetIntegrationOrderMax());
          }
        }
        //interface
        {
          int diag1, diag2;
          int ndiag1, ndiag2;
          if (v2cut_1[negvidx[0]] == v2cut_1[posvidx[0]] || v2cut_2[negvidx[0]] == v2cut_1[posvidx[0]])
          {
            diag1 = v2cut_1[posvidx[0]];
            ndiag1 = v2cut_2[posvidx[0]];
          }
          else
          {
            diag1 = v2cut_2[posvidx[0]];
            ndiag1 = v2cut_1[posvidx[0]];
          }

          if (v2cut_1[negvidx[1]] == v2cut_1[posvidx[1]] || v2cut_2[negvidx[1]] == v2cut_1[posvidx[1]])
          {
            diag2 = v2cut_1[posvidx[1]];
            ndiag2 = v2cut_2[posvidx[1]];
          }
          else
          {
            diag2 = v2cut_2[posvidx[1]];
            ndiag2 = v2cut_1[posvidx[1]];
          }

          Array< const Vec<SD> * > trig1(3);
          Array< const Vec<SD> * > trig2(3);

          trig1[0] = cutpoints[diag1]; 
          trig1[1] = cutpoints[ndiag1]; 
          trig1[2] = cutpoints[diag2]; 

          trig2[0] = cutpoints[diag2]; 
          trig2[1] = cutpoints[ndiag2]; 
          trig2[2] = cutpoints[diag1]; 
                  
          FillSimplexCoDim1WithRule<SD> ( trig1, numint.compquadrule.GetInterfaceRule(), 
                                          numint.GetIntegrationOrderMax());
          FillSimplexCoDim1WithRule<SD> ( trig2, numint.compquadrule.GetInterfaceRule(), 
                                          numint.GetIntegrationOrderMax());

        }

      } // end of 3 or 4 cutpoints
      else
      {
        cout << "cutpoints.Size() = " << cutpoints.Size() << endl;
        throw Exception(" did not expect this.. -3-");
      }
    }

    template<ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
    void CutSimplex<2,ET_SPACE,ET_TIME>::MakeQuad(const Simplex <2> & s, 
                                                  const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
    { 
      enum { SD = 2};

      // cout << " simplex = " << s << endl;

      Array< const Vec<SD> * > cutpoints(0);
      Array< const Vec<SD> * > pospoints(0);
      Array< const Vec<SD> * > negpoints(0);

      Array<int> posvidx(0);
      Array<int> negvidx(0);

      // vertex idx connected to cut idx (just in case of 4 cut positions)
      // connectivity information of cuts
      Array<int> v2cut_1(4); 
      Array<int> v2cut_2(4); 
      v2cut_1 = -1;
      v2cut_2 = -1;

      const int edge[3][2] = { {0, 1},
                               {0, 2},
                               {1, 2}};
              
      double vvals[3];
      bool zero[3];
              
      for (int j = 0; j < 3; ++j)
      {
        zero[j] = false;
        vvals[j] = numint.lset(*(s.p[j]));
        if (vvals[j] > 0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
        }
        else if (vvals[j] < 0)
        {
          negpoints.Append(numint.pc(*(s.p[j])));
          negvidx.Append(j);
        }
        else // (vvals[j] == 0.0)
        {
          pospoints.Append(numint.pc(*(s.p[j])));
          posvidx.Append(j);
          zero[j] = true;
        }
      }

      // cout << " Avvals = \n";
      // for (int i = 0; i < 3; ++i)
      //   cout << i << ":" << vvals[i] << endl;

      int cntcuts = 0;
      for (int j = 0; j < 3; ++j) //edges
      {
        const int lv = edge[j][0];
        const int rv = edge[j][1];
        const double valleft = vvals[lv];
        const double valright = vvals[rv];
        bool hascut = (valleft * valright < 0);
        if (zero[lv] && valright < 0)
          hascut = true;

        if (zero[rv] && valleft < 0)
          hascut = true;

        if (hascut)
        {
          const double cutpos = valleft / (valleft - valright);
          // std::cout << " cutpos = " << cutpos << std::endl;
          Vec<SD> p = (1-cutpos) * *(s.p[lv]) + cutpos * *(s.p[rv]) ;
          cutpoints.Append(numint.pc(p));
          // collect connectivity of cut and vertices
          if (v2cut_1[lv] == -1)
            v2cut_1[lv] = cntcuts;
          else
            v2cut_2[lv] = cntcuts;
          if (v2cut_1[rv] == -1)
            v2cut_1[rv] = cntcuts;
          else
            v2cut_2[rv] = cntcuts;
          cntcuts ++;
        }
      }

      // std::cout << " cutpoints[0] = " << *cutpoints[0] << std::endl;
      // std::cout << " cutpoints[1] = " << *cutpoints[1] << std::endl;

      if (cutpoints.Size() == 2) // three intersections: prism + tetra
      {

        Array< const Vec<SD> *> & minorgroup ( negpoints.Size() > pospoints.Size() ? 
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_minor = negpoints.Size() > pospoints.Size() ? POS : NEG;
        Array< const Vec<SD> *> & majorgroup ( negpoints.Size() <= pospoints.Size() ? 
                                               pospoints : negpoints);
        DOMAIN_TYPE dt_major = negpoints.Size() <= pospoints.Size() ? POS : NEG;

        Array<int> & majvidx( negvidx.Size() > posvidx.Size() ? negvidx : posvidx);

        for (int k = 0; k < 2; ++k)
          minorgroup.Append(cutpoints[k]);
        // minorgroup is a simplex of type dt_minor
        FillSimplexWithRule<SD>(minorgroup, 
                                numint.compquadrule.GetRule(dt_minor), 
                                numint.GetIntegrationOrderMax());

        Array< Simplex<SD> * > innersimplices(0);
        for (int k = 0; k < 2; ++k)
        {
          int corresponding_cut = v2cut_1[majvidx[k]];
          majorgroup.Append(cutpoints[corresponding_cut]);
        }

        DecomposePrismIntoSimplices<SD>(majorgroup, innersimplices, numint.pc, numint.lh);
        for (int l = 0; l < innersimplices.Size(); ++l)
        {
          FillSimplexWithRule<SD>(innersimplices[l]->p, 
                                  numint.compquadrule.GetRule(dt_major), 
                                  numint.GetIntegrationOrderMax());
        }
                
        // and the interface:
        FillSimplexCoDim1WithRule<SD> ( cutpoints, numint.compquadrule.GetInterfaceRule(), 
                                        numint.GetIntegrationOrderMax());

        
      }
      else
      {
        cout << "cutpoints.Size() = " << cutpoints.Size() << endl;
        throw Exception(" did not expect this.. -2-");
      }
    }
  }



  template <int D, ELEMENT_TYPE ET_SPACE, ELEMENT_TYPE ET_TIME>
  void MakeQuadRuleOnCutSimplex(const Simplex <D> & s, 
                                const NumericalIntegrationStrategy<ET_SPACE,ET_TIME> & numint)
  {
    DecompositionRules::CutSimplex<D,ET_SPACE,ET_TIME>::MakeQuad(s,numint);
  }


} // end of namespace
