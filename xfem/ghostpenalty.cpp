#define FILE_GHOSTPENALTY_CPP
#include "ghostpenalty.hpp"
#include <diffop_impl.hpp>

namespace ngfem
{

  class CentralFDStencils
  {
    static const int max_accuracy_half = 8; // order 2,4,6,8,..,16
    static const int max_order = 10; // order 1,2,3,4,..,10
    Table<double> * stencils;
  public:
    static CentralFDStencils & Instance()
    {
      static CentralFDStencils myInstance;
      return myInstance;
    }
    static const FlatVector<> Get(int order, int accuracy)
    {
      CentralFDStencils & instance = Instance();
      int row = order == 0 ? 0 : 1+(order-1)*max_accuracy_half+(accuracy-1)/2;
      const FlatArray<double> fa ((*(instance.stencils))[row]);
      FlatVector<> coeffs(fa.Size(),&fa[0]);
      return coeffs;
    }
    static double GetOptimalEps(int order, int accuracy)
    {
      const double EPS = std::numeric_limits<double>::epsilon();
      const int stencil_width = order == 0 ? 0 : (accuracy+1)/2 + int(order-0.5)/2;
      const double eps = order == 0 ? 1.0 : std::pow((2 * stencil_width + 1)*EPS,1.0/(accuracy + order)); // + 0.25;
      static bool first0 = true;
      if (first0 && order==0)
      {

        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first0 = false;
        // getchar();
      }
      static bool first1 = true;
      if (first1 && order==1)
      {

        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first1 = false;
        // getchar();
      }
      static bool first2 = true;
      if (first2 && order==2)
      {
        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first2 = false;
        // getchar();
      }
      static bool first3 = true;
      if (first3 && order==3)
      {
        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first3 = false;
        // getchar();
      }
      static bool first4 = true;
      if (first4 && order==4)
      {
        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first4 = false;
        // getchar();
      }
      static bool first5 = true;
      if (first5 && order==5)
      {
        cout << IM(3) << " order, eps = " << order << ", " << eps << endl;
        first5 = false;
        // getchar();
      }
      return eps;

    }

  protected:
    CentralFDStencils()
    {
      const int rows = max_accuracy_half*max_order+1;

      Array<int> cnt(rows);

      cnt[0] = 1;
      for (int i = 0; i < max_order; ++i)
        for (int j = 0; j < max_accuracy_half; ++j)
        {
          const int row = i * max_accuracy_half + j + 1;
          const int order = i+1;
          // const int accuracy = 2*(j+1);
          const int stencil_width = j+1 + int(order-0.5)/2;
          cnt[row] = 2 * stencil_width +1;
        }

      stencils = new Table<double>(cnt);

      (*stencils)[0][0] = 1.0;
      for (int i = 0; i < max_order; ++i)
        for (int j = 0; j < max_accuracy_half; ++j)
        {
          const int row = i * max_accuracy_half + j + 1;
          const int order = i+1;
          // const int accuracy = 2*(j+1);
          const int stencil_width = j+1 + int(order-0.5)/2;
          const int size = 2 * stencil_width + 1;

          Vector<> stencil(size);
          {
            Matrix<> A(size);
            Matrix<> invA(size);
            Vector<> f(size);
            double inv_factorial=1.0;
            for (int k = 0; k < size; k++)
            {
              if (k>0)
                inv_factorial/=k;
              for (int l = 0; l < size; l++)
                A(k,l) = std::pow(l-stencil_width,k) * inv_factorial;
              f(k) = k==order ? 1.0 : 0;
            }
            CalcInverse(A,invA);
            stencil = invA * f;
          }

          for (int col = 0; col < size; ++col)
            (*stencils)[row][col] = stencil(col);
        }
    }

    ~CentralFDStencils()
    {
      delete stencils;
    }
  };


  template <int D, int ORDER>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDuDnkHDiv<D,ORDER>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                                 MAT & mat, LocalHeap & lh)
  {
    const int FD_ACCURACY = 4;
    const HDivFiniteElement<D> & hdivfel =
      dynamic_cast<const HDivFiniteElement<D> & > (bfel);
    const int ndof = hdivfel.GetNDof();


    Vec<D> normal = static_cast<const DimMappedIntegrationPoint<D>&>(mip).GetNV();

    // cout << "normal: " << normal << endl;
    // Vec<D> normal; normal(0) = -1.0; normal(1) = 1.0;
    // normal /= L2Norm(normal);

    Vec<D> invjac_normal = mip.GetJacobianInverse() * normal;

    const double h = D==2 ? sqrt(mip.GetJacobiDet()) : cbrt(mip.GetJacobiDet());

    FlatVector<> fdstencil (CentralFDStencils::Get(ORDER,FD_ACCURACY));
    const double eps = h * CentralFDStencils::GetOptimalEps(ORDER,FD_ACCURACY);

    const int stencilpoints = fdstencil.Size();
    const int stencilwidth = (stencilpoints-1)/2;

    FlatMatrixFixWidth<D> shape (ndof, lh);
    mat = 0.0;
    Vec<D> normal_dir_wrt_ref = invjac_normal;
    const double eps_fac = std::pow(1.0/eps,ORDER);
    for (int i = 0; i < stencilpoints; ++i)
    {
      Vec<D> vec = mip.GetPoint();
      vec += (i-stencilwidth) * eps * normal;

      IntegrationPoint ip_x0(mip.IP());
      for (int d = 0; d < D; ++d)
        ip_x0(d) += (i-stencilwidth) * eps * normal_dir_wrt_ref(d);
      MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());

      Vec<D> diff = vec - mip_x0.GetPoint();
      Vec<D> update = 0;
      int its = 0;
      while (L2Norm(diff) > 1e-8*h && its < 20)
      {
        MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());
        diff = vec - mip_x0.GetPoint();
        // cout << "mip_x0.GetPoint()" << mip_x0.GetPoint() << endl;
        // cout << "diff" << diff << endl;
        update = mip_x0.GetJacobianInverse() * diff;
        for (int d = 0; d < D; ++d)
          ip_x0(d) += update(d);
        // getchar();
        its++;
      }
      if (its >= 50)
        cerr << "its >= 50 " << endl;

      MappedIntegrationPoint<D,D> mip_now(ip_x0,mip.GetTransformation());
      hdivfel.CalcMappedShape (mip_now, shape);
      mat += eps_fac * fdstencil(i) * shape;
    }

    if (false)
    {
      const int feorder = hdivfel.Order();
      cout << "ORDER = " << ORDER << endl;
      cout << mip.GetPoint() << endl;
      // cout << "dshapedn" << endl;
      // cout << shape << endl;
      cout << "coeffs" << endl;
      cout << fdstencil << endl;
      cout << "dshapednk" << endl;
      cout << mat << endl;
      getchar();
    }
  }



  template <int D, int ORDER>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpDuDnk<D,ORDER>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                             MAT & mat, LocalHeap & lh)
  {
    const int FD_ACCURACY = 4;
    int version = 2;

    if (version == 1)
    // not higher order accurate on curved meshes (!),
    // but more stable and efficient (one derivate less to evaluate by FD)
    {
      const ScalarFiniteElement<D> & scafe =
        dynamic_cast<const ScalarFiniteElement<D> & > (bfel);
      const int ndof = scafe.GetNDof();

      Vec<D> normal = static_cast<const DimMappedIntegrationPoint<D>&>(mip).GetNV();

      // cout << "normal: " << normal << endl;
      // Vec<D> normal; normal(0) = -1.0; normal(1) = 1.0;
      // normal /= L2Norm(normal);

      Vec<D> invjac_normal = mip.GetJacobianInverse() * normal;

      FlatVector<> fdstencil (CentralFDStencils::Get(ORDER-1,FD_ACCURACY));
      const double eps = CentralFDStencils::GetOptimalEps(ORDER-1,FD_ACCURACY);
      const int stencilpoints = fdstencil.Size();
      const int stencilwidth = (stencilpoints-1)/2;
      FlatMatrix<> dshapedn  (ndof, stencilpoints, lh);
      FlatVector<> dshapednk (ndof, lh);

      const double norminvjacn = L2Norm(invjac_normal);
      Vec<D> normal_dir_wrt_ref = (1.0/norminvjacn)  * invjac_normal;


      for (int i = 0; i < stencilpoints; ++i)
      {
        IntegrationPoint ip(mip.IP());
        for (int d = 0; d < D; ++d)
          ip(d) += (i-stencilwidth) * eps * normal_dir_wrt_ref(d);
        dshapedn.Col(i) = scafe.GetDShape (ip, lh) * invjac_normal;
      }
      dshapednk = dshapedn * fdstencil;
      const double eps_fac = std::pow(norminvjacn/eps,ORDER-1);
      mat.Row(0) = eps_fac * dshapednk;

      const int feorder = scafe.Order();
      if (false && feorder > 1)
      {
        cout << "ORDER = " << ORDER << endl;
        cout << mip.GetPoint() << endl;
        cout << "dshapedn" << endl;
        cout << dshapedn << endl;
        cout << "coeffs" << endl;
        cout << fdstencil << endl;
        cout << "dshapednk" << endl;
        cout << dshapednk << endl;
        getchar();
      }
    }
    else
    {
      // VERSION 2 (fix points on phys. domain + transform points(!) back)
      // higher order accurate also on curved meshes (!),
      // but takes all derivatives by FD
      const ScalarFiniteElement<D> & scafe =
        dynamic_cast<const ScalarFiniteElement<D> & > (bfel);
      const int ndof = scafe.GetNDof();

      Vec<D> normal = static_cast<const DimMappedIntegrationPoint<D>&>(mip).GetNV();

      // cout << "normal: " << normal << endl;
      // Vec<D> normal; normal(0) = -1.0; normal(1) = 1.0;
      // normal /= L2Norm(normal);

      Vec<D> invjac_normal = mip.GetJacobianInverse() * normal;

      const double h = D==2 ? sqrt(mip.GetJacobiDet()) : cbrt(mip.GetJacobiDet());

      FlatVector<> fdstencil (CentralFDStencils::Get(ORDER,FD_ACCURACY));
      const double eps = h * CentralFDStencils::GetOptimalEps(ORDER,FD_ACCURACY);

      const int stencilpoints = fdstencil.Size();
      const int stencilwidth = (stencilpoints-1)/2;

      FlatMatrix<> shape  (ndof, stencilpoints, lh);
      FlatVector<> dshapednk (ndof, lh);

      Vec<D> normal_dir_wrt_ref = invjac_normal;

      for (int i = 0; i < stencilpoints; ++i)
      {
        Vec<D> vec = mip.GetPoint();
        vec += (i-stencilwidth) * eps * normal;

        IntegrationPoint ip_x0(mip.IP());
        for (int d = 0; d < D; ++d)
          ip_x0(d) += (i-stencilwidth) * eps * normal_dir_wrt_ref(d);
        MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());

        Vec<D> diff = vec - mip_x0.GetPoint();
        Vec<D> update = 0;
        int its = 0;
        while (L2Norm(diff) > 1e-8*h && its < 20)
        {
          MappedIntegrationPoint<D,D> mip_x0(ip_x0,mip.GetTransformation());
          diff = vec - mip_x0.GetPoint();
          // cout << "mip_x0.GetPoint()" << mip_x0.GetPoint() << endl;
          // cout << "diff" << diff << endl;
          update = mip_x0.GetJacobianInverse() * diff;
          for (int d = 0; d < D; ++d)
            ip_x0(d) += update(d);
          // getchar();
          its++;
        }
        if (its >= 50)
          cerr << "its >= 50 " << endl;
        shape.Col(i) = scafe.GetShape (ip_x0, lh);
      }

      dshapednk = shape * fdstencil;
      const double eps_fac = std::pow(1.0/eps,ORDER);
      mat.Row(0) = eps_fac * dshapednk;

      const int feorder = scafe.Order();
      if (false && feorder > 1)
      {
        cout << "ORDER = " << ORDER << endl;
        cout << mip.GetPoint() << endl;
        cout << "dshapedn" << endl;
        cout << shape << endl;
        cout << "coeffs" << endl;
        cout << fdstencil << endl;
        cout << "dshapednk" << endl;
        cout << dshapednk << endl;
        getchar();
      }
      // TODO: Add version 3:
      // VERSION 3 (fix points on phys. domain + transform points(!) back)
      // take du/dn (n the path-adjusted) normal direction and take FD from this.

    }


  }// generate matrix

  template class T_DifferentialOperator<DiffOpDuDnk<2,1>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,2>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,3>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,4>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,5>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,6>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,7>>;
  template class T_DifferentialOperator<DiffOpDuDnk<2,8>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,1>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,2>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,3>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,4>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,5>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,6>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,7>>;
  template class T_DifferentialOperator<DiffOpDuDnk<3,8>>;

  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,1>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,2>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,3>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,4>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,5>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,6>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,7>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<2,8>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,1>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,2>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,3>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,4>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,5>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,6>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,7>>;
  template class T_DifferentialOperator<DiffOpDuDnkHDiv<3,8>>;

}
