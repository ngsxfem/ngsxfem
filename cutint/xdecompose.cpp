#include "xdecompose.hpp"

namespace xintegration
{


  template<int SD>
  double Measure(const Array< const Vec<SD> *> & s)
  {
      cout << IM(1) << " not implemented for SD=" << SD << endl;
    throw Exception("not implemented");
    return 0;
  }

  template<>
  double Measure<1,1>(const Array< const Vec<1> *> & s)
  {
    Vec<1> a = *(s[1]) - *(s[0]);
    return L2Norm(a);
  }

  template<>
  double Measure<1,2>(const Array< const Vec<2> *> & s)
  {
    Vec<2> a = *(s[1]) - *(s[0]);
    return L2Norm(a);
  }

  template<>
  double Measure<2,2>(const Array< const Vec<2> *> & s)
  {
    Vec<2> a = *s[1] - *s[0];
    Vec<2> b = *s[2] - *s[0];
    return abs(a(0)*b(1)-b(0)*a(1)) / 2.0;
  }

  template<>
  double Measure<2,3>(const Array< const Vec<3> *> & s)
  {
    Vec<3> a = *s[1] - *s[0];
    Vec<3> b = *s[2] - *s[0];
    Vec<3> c = Cross(a,b);
    return L2Norm(c) / 2.0;
  }

  template<>
  double Measure<3,3>(const Array< const Vec<3> *> & s)
  {
    Vec<3> a = *s[1] - *s[0];
    Vec<3> b = *s[2] - *s[0];
    Vec<3> c = *s[3] - *s[0];
    return abs(Determinant(a,b,c)) / 6.0;
  }

}
