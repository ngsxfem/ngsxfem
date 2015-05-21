#include "traceintegrators.hpp"

namespace ngfem
{

  template<int D>
  void TraceMassIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("TraceMassIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&base_fel);

    elmat = 0.0;
    if (!xfe) return;

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    int ndof = scafe.GetNDof();
    FlatVector<> shape(ndof,lh);

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;
      const double weight = fquad.weights(i) * len;

      const double coef_val = coef->Evaluate(mip);
	    
      scafe.CalcShape (mip.IP(),shape);
      elmat += (coef_val * weight) * shape * Trans(shape);
    }
  }

  template class TraceMassIntegrator<2>;
  template class TraceMassIntegrator<3>;

  static RegisterBilinearFormIntegrator<TraceMassIntegrator<2> > inittracemass2d ("tracemass", 2, 1);
  static RegisterBilinearFormIntegrator<TraceMassIntegrator<3> > inittracemass3d ("tracemass", 3, 1);

  template<int D>
  void TraceSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
                       const ElementTransformation & eltrans,
                       FlatVector<double> elvec,
                       LocalHeap & lh) const
  {
    static Timer timer ("TraceSourceIntegrator::CalcElementVector");
    RegionTimer reg (timer);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&base_fel);

    elvec = 0.0;
    if (!xfe) return;

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    int ndof = scafe.GetNDof();
    FlatVector<> shape(ndof,lh);

    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;
      const double weight = fquad.weights(i) * len;

      const double coef_val = coef->Evaluate(mip);
      scafe.CalcShape (mip.IP(),shape);
      elvec += (coef_val * weight) * shape;
    }
  }

  template class TraceSourceIntegrator<2>;
  template class TraceSourceIntegrator<3>;

  static RegisterLinearFormIntegrator<TraceSourceIntegrator<2> > inittracesource2d ("tracesource", 2, 1);
  static RegisterLinearFormIntegrator<TraceSourceIntegrator<3> > inittracesource3d ("tracesource", 3, 1);


// ---------------------

  template <int D>
  template <typename FEL, typename MIP, typename MAT>
  void DiffOpEvalExtTrace<D>::GenerateMatrix (const FEL & bfel, const MIP & mip,
                                          MAT & mat, LocalHeap & lh)
  {
    const XFiniteElement * xfe = dynamic_cast<const XFiniteElement *> (&bfel);
    if (!xfe)
      return;
    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());
    FlatVector<> shape (scafe.GetNDof(), &mat(0,0));
    shape = scafe.GetShape(mip.IP(), lh);
  }


  template <int D>  ExtTraceIntegrator<D> :: ExtTraceIntegrator  (shared_ptr<CoefficientFunction> coeff)
    : T_BDBIntegrator<DiffOpEvalExtTrace<D>, DiagDMat<1>, CompoundFiniteElement > (DiagDMat<1> (coeff))
  { ; }

  template <int D>  ExtTraceIntegrator<D> :: ExtTraceIntegrator  (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : T_BDBIntegrator<DiffOpEvalExtTrace<D>, DiagDMat<1>, CompoundFiniteElement > (coeffs)
  { ; }

  template <int D>  ExtTraceIntegrator<D> :: ~ExtTraceIntegrator () { ; }

  template class ExtTraceIntegrator<2>;
  template class ExtTraceIntegrator<3>;

  static RegisterBilinearFormIntegrator<ExtTraceIntegrator<2> > initexttracemass0 ("exttrace", 2, 1);
  static RegisterBilinearFormIntegrator<ExtTraceIntegrator<3> > initexttracemass1 ("exttrace", 3, 1);

// ---------------------

  template<int D>
  void TraceLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("TraceLaplaceIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&base_fel);

    elmat = 0.0;
    if (!xfe) return;

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    int ndof = scafe.GetNDof();
    //FlatVector<> shape(ndof,lh);
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;
      const double weight = fquad.weights(i) * len;

      const double coef_val = coef->Evaluate(mip);

      scafe.CalcMappedDShape(mip, dshape);

      //scafe.CalcShape (mip.IP(),shape);
      elmat += (coef_val * weight) * dshape * Trans(dshape);
    }
  }

  template class TraceLaplaceIntegrator<2>;
  template class TraceLaplaceIntegrator<3>;

  static RegisterBilinearFormIntegrator<TraceLaplaceIntegrator<2> > inittracelaplace2d ("tracelaplace", 2, 1);
  static RegisterBilinearFormIntegrator<TraceLaplaceIntegrator<3> > inittracelaplace3d ("tracelaplace", 3, 1);




// ---------------------

  template<int D>
  void TraceLaplaceBeltramiIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("TraceLaplaceBeltramiIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&base_fel);

    elmat = 0.0;
    if (!xfe) return;

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    int ndof = scafe.GetNDof();
    //FlatVector<> shape(ndof,lh);
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D> proj(ndof,lh);
    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;
      const double weight = fquad.weights(i) * len;

      const double coef_val = coef->Evaluate(mip);

      scafe.CalcMappedDShape(mip, dshape);
      proj = dshape * (Id<D>() - normal * Trans(normal)) ;
      elmat += (coef_val * weight) * proj * Trans(proj);
    }
  }

  template class TraceLaplaceBeltramiIntegrator<2>;
  template class TraceLaplaceBeltramiIntegrator<3>;

  static RegisterBilinearFormIntegrator<TraceLaplaceBeltramiIntegrator<2> > inittracelaplacebeltrami2d ("tracelaplacebeltrami", 2, 1);
  static RegisterBilinearFormIntegrator<TraceLaplaceBeltramiIntegrator<3> > inittracelaplacebeltrami3d ("tracelaplacebeltrami", 3, 1);

// ---------------------

  template<int D>
  void TraceConvectionIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans,
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("TraceConvectionIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const XFiniteElement * xfe =
      dynamic_cast<const XFiniteElement *> (&base_fel);

    elmat = 0.0;
    if (!xfe) return;

    const ScalarFiniteElement<D> & scafe =
      dynamic_cast<const ScalarFiniteElement<D> & > (xfe->GetBaseFE());

    int ndof = scafe.GetNDof();
    FlatVector<> shape(ndof,lh);
    FlatMatrixFixWidth<D> dshape(ndof,lh);
    FlatMatrixFixWidth<D> proj(ndof,lh);
    const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
    const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
    const FlatQuadratureRuleCoDim1<D> & fquad(fcompr.GetInterfaceRule());

    for (int i = 0; i < fquad.Size(); ++i)
    {
      IntegrationPoint ip(&fquad.points(i,0),0.0);
      MappedIntegrationPoint<D,D> mip(ip, eltrans);
      Vec<D> convvec;
      conv->Evaluate(mip,convvec);

      Mat<D,D> Finv = mip.GetJacobianInverse();
      const double absdet = mip.GetMeasure();

      Vec<D> nref = fquad.normals.Row(i);
      Vec<D> normal = absdet * Trans(Finv) * nref ;
      double len = L2Norm(normal);
      normal /= len;
      const double weight = fquad.weights(i) * len;

      scafe.CalcShape (mip.IP(),shape);
      scafe.CalcMappedDShape(mip, dshape);
      proj = dshape * (Id<D>() - normal * Trans(normal)) ;
      elmat += weight * shape* Trans(proj * convvec);
    }
  }

  template class TraceConvectionIntegrator<2>;
  template class TraceConvectionIntegrator<3>;

  static RegisterBilinearFormIntegrator<TraceConvectionIntegrator<2> > inittraceconvection2d ("traceconvection", 2, 1);
  static RegisterBilinearFormIntegrator<TraceConvectionIntegrator<3> > inittraceconvection3d ("traceconvection", 3, 1);

}
