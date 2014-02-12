#include "xintegration.hpp"

namespace xintegration
{


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
  void FillSimplexWithRule<3> (const Simplex<3> & s, QuadratureRule<3> & quaddom, int intorder)
  {
    Vec<3> a = *s.p[1] - *s.p[0];
    Vec<3> b = *s.p[2] - *s.p[0];
    Vec<3> c = *s.p[3] - *s.p[0];
    const double trafofac = abs(Determinant(a,b,c));
    const IntegrationRule & ir = SelectIntegrationRule (ET_TET, intorder);

    for (int k = 0; k < ir.GetNIP(); k++)
    {
      Vec<3> point(0.0);
      point = *s.p[0];
      for (int m = 0; m < 3 ;++m)
        point += ir[k](m) * (*s.p[m+1]);
      const double weight = ir[k].Weight() * trafofac;
      quaddom.points.Append(point);
      quaddom.weights.Append(weight);
    }
  }

  template<>
  void FillSimplexWithRule<4> (const Simplex<4> & s, QuadratureRule<4> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 4");
  }

  template<>
  void FillSimplexWithRule<2> (const Simplex<2> & s, QuadratureRule<2> & quaddom, int intorder)
  {
    throw Exception(" nonono - still no 2");
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


  template class NumericalIntegrationStrategy<ET_TRIG, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TET, ET_SEGM>;
  template class NumericalIntegrationStrategy<ET_TRIG, ET_POINT>;
  template class NumericalIntegrationStrategy<ET_TET, ET_POINT>;


} // end of namespace
