#define FILE_SHIFTEDEVALUATE_CPP
#include "shiftedevaluate.hpp"
#include <diffop_impl.hpp>


namespace ngfem
{



  void DiffOpShiftedEval ::
  CalcMatrix (const FiniteElement & bfel,
              const BaseMappedIntegrationPoint & bmip,
              SliceMatrix<double,ColMajor> mat,
              LocalHeap & lh) const
  {
    const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE> & mip =
      static_cast<const MappedIntegrationPoint<DIM_ELEMENT,DIM_SPACE>&> (bmip);

    const ScalarFiniteElement<2> & scafe =
            dynamic_cast<const ScalarFiniteElement<2> & > (bfel);
    const int ndof = scafe.GetNDof();

    FlatVector<> shape (ndof,lh);

    IntegrationPoint ip(mip.IP());
    auto elid = mip.GetTransformation().GetElementId();
    Array<int> dnums;

    back->GetFESpace()->GetDofNrs(elid,dnums);
    FlatVector<> values_back(dnums.Size()*DIM_SPACE,lh);
    FlatMatrixFixWidth<2> vector_back(dnums.Size(),&(values_back(0)));
    back->GetVector().GetIndirect(dnums,values_back);

    forth->GetFESpace()->GetDofNrs(elid,dnums);
    FlatVector<> values_forth(dnums.Size()*DIM_SPACE,lh);
    FlatMatrixFixWidth<2> vector_forth(dnums.Size(),&(values_forth(0)));
    forth->GetVector().GetIndirect(dnums,values_forth);

    FiniteElement& fe_back = back->GetFESpace()->GetFE(elid,lh);
    const ScalarFiniteElement<2> & scafe_back =
            dynamic_cast<const ScalarFiniteElement<2> & > (fe_back);
    FlatVector<> shape_back(dnums.Size(),lh);
    scafe_back.CalcShape(ip,shape_back);
    Vec<2> dvec_back = Trans(vector_back)*shape_back;

    FiniteElement& fe_forth = forth->GetFESpace()->GetFE(elid,lh);
    const ScalarFiniteElement<2> & scafe_forth =
            dynamic_cast<const ScalarFiniteElement<2> & > (fe_forth);
    FlatVector<> shape_forth(dnums.Size(),lh);
    scafe_forth.CalcShape(ip,shape_forth);
    Vec<2> dvec_forth = Trans(vector_forth)*shape_forth;

    Vec<2> z = mip.GetPoint() + dvec_forth;

    Vec<2> x = ip.Point();
    int its = 0;
    const double h = sqrt(mip.GetJacobiDet());
    Vec<2> diff = mip.GetPoint() + dvec_back - z;
    IntegrationPoint ipx(x);
    while (its==0 || (L2Norm(diff) > 1e-8*h && its < 200))
    {
      MappedIntegrationPoint<2,2> mip_x0(ipx,mip.GetTransformation());
      scafe_back.CalcShape(ipx,shape_back);
      dvec_back = Trans(vector_back)*shape_back;
      diff = z - dvec_back - mip_x0.GetPoint();
      Mat<2> dphi = mip_x0.GetJacobian();
      FlatMatrixFixWidth<2> dshape_back(dnums.Size(),lh);
      scafe_back.CalcDShape(ipx,dshape_back);
      Mat<2> dback_ref = Trans(dshape_back)*vector_back;
      Mat<2> dback = Trans(mip_x0.GetJacobianInverse())*dback_ref + dphi;
      Mat<2> invd;
      CalcInverse(dback,invd);

      // Backtracking-Armijo line search

      // scalars
      double alpha = 1;
      double beta = 0.1;
      int itsb = 0;
      // descent direction
      Vec<2> ddirection = invd * diff;
      // Jacobi matrix
      Mat<2> jacobi_old = dback;
      // function value at current Newton iterate
      Vec<2> value_old = -diff;
      // Stuff to calculate function value at
      // prospective new iterate
      IntegrationPoint ipx_new(x);
      FlatVector<> shape_back_new(dnums.Size(),lh);
      // the difference in the Armijo condition of which all components
      // need to be <= 0
      Vec<2> diffb;

      while(itsb == 0 || ( ( (diffb(0) > 1e-8*h) || (diffb(1) > 1e-8*h) ) && itsb < 20) )
      {
        if(itsb != 0)
            alpha *= 0.5;
        ipx_new.Point() = x + alpha*ddirection;
        MappedIntegrationPoint<2,2> mip_xnew(ipx_new,mip.GetTransformation());
        scafe_back.CalcShape(ipx_new,shape_back_new);
        Vec<2> dvec_back_new = Trans(vector_back)*shape_back_new;
        Vec<2> value_new = -z + mip_xnew.GetPoint() + dvec_back_new;
        diffb = value_new - value_old - alpha*beta*jacobi_old*ddirection;
        itsb ++;
        cout << "itsb = " << itsb << endl;
        cout << "diffb(0) = " << diffb(0) << endl;
        cout << "diffb(1) = " << diffb(1) << endl;
      }

      cout << "alpha = " << alpha << endl;


      its++;
      //alpha = 0.25;
      x += alpha*invd * diff;

      cout << "x = " << x << endl;
      cout << "z = " << z << endl;
      cout << "diff = " << diff << endl;
      cout << "its = " << its << endl;
      cout << "dvec_back = " << dvec_back << endl;
      cout << "dvec_forth = " << dvec_forth << endl;
      cout << "dphi = " << dphi << endl;
      cout << "dback_ref = " << dback_ref << endl;
      cout << "dback = " << dback << endl;
      //getchar();


      ipx.Point() = x;

    }



    scafe.CalcShape(ipx,shape);
    mat = 0.0;
    mat.Row(0) = shape;
  }

  void DiffOpShiftedEval ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x,
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }



}


