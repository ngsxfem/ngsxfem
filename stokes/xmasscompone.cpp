#include "xmasscompone.hpp"

namespace ngfem
{

  template<int D>
  void xMassCompOne<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("XMassCompOneInt::CalcElMat");
    RegionTimer reg (timer);

    const ScalarFiniteElement<D> * scafe;
    const XFiniteElement * xfe;
    const XDummyFE * dummfe;
    CastXScalarFiniteElements(base_fel, scafe, xfe, dummfe);

    elmat = 0.0;
    
    int ndof_x = xfe!=NULL ? xfe->GetNDof() : 0;
    int ndof = scafe->GetNDof();
    int ndof_total = ndof+ndof_x;
    
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape(ndof,&shape_total(0));
    FlatVector<> shapex(ndof_x,&shape_total(ndof));

    int p = scafe->Order();
    
    for (auto dt : {POS,NEG})
    {
      if(!xfe)
      { 
        if (dummfe->GetDomainType() != dt)
          continue;
        IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
        for (int i = 0 ; i < pir.GetNIP(); i++)
        {
          MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          shape.Range(0,ndof) = scafe->GetShape(mip.IP(), lh);
          elmat += mip.GetWeight() * coef * shape * Trans(shape);
        }
      }
      else
      {
        const FlatXLocalGeometryInformation & xgeom(xfe->GetFlatLocalGeometry());
        const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
        const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));
        for (int i = 0; i < fquad.Size(); ++i)
        {
          IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
          shape = scafe->GetShape(mip.IP(), lh);
          shapex *= 0.0;

          for (int l = 0; l < ndof_x; ++l)
            if (xfe->GetSignsOfDof()[l] != dt)
              shapex(l) = 0.0;
          elmat += mip.GetWeight() * coef * shape_total * Trans(shape_total);
        }
      }
    }
  }

  template class xMassCompOne<2>;
  template class xMassCompOne<3>;

  static RegisterBilinearFormIntegrator<xMassCompOne<2>> initxmasscompone2d ("xmasscompone", 2, 2);
  static RegisterBilinearFormIntegrator<xMassCompOne<3>> initxmasscompone3d ("xmasscompone", 3, 2);

}
