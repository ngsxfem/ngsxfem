#include "xfemIntegrators.hpp"

namespace ngfem
{

  template<int D>
  void FictXSourceIntegrator<D> ::
  CalcElementVector (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatVector<double> & elvec,
		     LocalHeap & lh) const
  {
    static Timer timer ("FictXSourceIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe[2]; xfe[0] = NULL; xfe[1] = NULL;
    const XDummyFE * dummfe[2]; dummfe[0] = NULL; dummfe[1] = NULL;

    const ScalarFiniteElement<D> * scafe = NULL;

    xfe[0] = dynamic_cast<const XFiniteElement* >(&cfel[0]);
    xfe[1] = dynamic_cast<const XFiniteElement* >(&cfel[1]);
    dummfe[0] = dynamic_cast<const XDummyFE* >(&cfel[0]);
    dummfe[1] = dynamic_cast<const XDummyFE* >(&cfel[1]);

    if (xfe[0] != NULL)
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[0]->GetBaseFE()));
    else
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[1]->GetBaseFE()));
    
    elvec = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    // if(xfe[0])
    //   std::cout << " xfe[NEG] " << std::endl;
    // if(xfe[1])
    //   std::cout << " xfe[POS] " << std::endl;

    int ndof_neg = xfe[0]!=NULL ? scafe->GetNDof() : 0;
    int ndof_pos = xfe[1]!=NULL ? scafe->GetNDof() : 0;
    int ndof_total = ndof_neg+ndof_pos;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape_neg(ndof_neg,&shape_total(0));
    FlatVector<> shape_pos(ndof_pos,&shape_total(ndof_neg));

    FlatVector<> elvec_neg(ndof_neg,&elvec(0));
    FlatVector<> elvec_pos(ndof_pos,&elvec(ndof_neg));
    
    int p = scafe->Order();
    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      // std::cout << " dt = " << dt << std::endl;
      FlatVector<> & shape_loc = dt == NEG ? shape_neg : shape_pos;
      FlatVector<> & elvec_loc = dt == NEG ? elvec_neg : elvec_pos;

      shape_total = 0.0;

      const int idxfe = dt==POS ? 1 : 0;

      if( xfe[idxfe]) 
      {
        // if (! xfe[dt]->Empty())
        {
          const FlatXLocalGeometryInformation & xgeom(xfe[idxfe]->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

          if (fquad.Size() > 0)
          {
            for (int i = 0; i < fquad.Size(); ++i)
            {
            
              IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
              MappedIntegrationPoint<D,D> mip(ip, eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
              shape_loc = scafe->GetShape(mip.IP(), lh);

              double fac = mip.GetWeight();
              elvec_loc += (fac*coef) * shape_loc;

            } // quad rule
          }
          else
          {
            IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
              shape_loc = scafe->GetShape(mip.IP(), lh);
              double fac = mip.GetWeight();
              elvec_loc += (fac*coef) * shape_loc;
            }
          }
        }
      } // if xfe
    } // loop over domain types

    // std::cout << " elvec = " << elvec << std::endl; getchar();
  }

  template class FictXSourceIntegrator<2>;
  template class FictXSourceIntegrator<3>;

  static RegisterLinearFormIntegrator<FictXSourceIntegrator<2> > initfictxsource2d ("fictxsource", 2, 2);
  static RegisterLinearFormIntegrator<FictXSourceIntegrator<3> > initfictxsource3d ("fictxsource", 3, 2);



  template<int D>
  void FictXMassIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> & elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("FictXMassIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe[2]; xfe[0] = NULL; xfe[1] = NULL;
    const XDummyFE * dummfe[2]; dummfe[0] = NULL; dummfe[1] = NULL;

    const ScalarFiniteElement<D> * scafe = NULL;

    xfe[0] = dynamic_cast<const XFiniteElement* >(&cfel[0]);
    xfe[1] = dynamic_cast<const XFiniteElement* >(&cfel[1]);
    dummfe[0] = dynamic_cast<const XDummyFE* >(&cfel[0]);
    dummfe[1] = dynamic_cast<const XDummyFE* >(&cfel[1]);

    if (xfe[0] != NULL)
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[0]->GetBaseFE()));
    else
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[1]->GetBaseFE()));
    
    elmat = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    // if(xfe[0])
    //   std::cout << " xfe[NEG] " << std::endl;
    // if(xfe[1])
    //   std::cout << " xfe[POS] " << std::endl;

    int ndof_neg = xfe[0]!=NULL ? scafe->GetNDof() : 0;
    int ndof_pos = xfe[1]!=NULL ? scafe->GetNDof() : 0;
    int ndof_total = ndof_neg+ndof_pos;
    FlatVector<> shape_total(ndof_total,lh);
    FlatVector<> shape_neg(ndof_neg,&shape_total(0));
    FlatVector<> shape_pos(ndof_pos,&shape_total(ndof_neg));

    FlatMatrix<> elmat_neg(ndof_neg,ndof_neg,lh);
    FlatMatrix<> elmat_pos(ndof_pos,ndof_pos,lh);
    
    elmat_neg = 0.0;
    elmat_pos = 0.0;

    int p = scafe->Order();
    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      // std::cout << " dt = " << dt << std::endl;
      FlatVector<> & shape_loc = dt == NEG ? shape_neg : shape_pos;
      FlatMatrix<> & elmat_loc = dt == NEG ? elmat_neg : elmat_pos;

      shape_total = 0.0;

      const int idxfe = dt==POS ? 1 : 0;

      if( xfe[idxfe]) 
      {
        // if (! xfe[dt]->Empty())
        {
          const FlatXLocalGeometryInformation & xgeom(xfe[idxfe]->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

          if (fquad.Size() > 0)
          {
            for (int i = 0; i < fquad.Size(); ++i)
            {
            
              IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
              MappedIntegrationPoint<D,D> mip(ip, eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
              shape_loc = scafe->GetShape(mip.IP(), lh);

              double fac = mip.GetWeight();
              elmat_loc += (fac*coef) * shape_loc * Trans(shape_loc);

            } // quad rule
          }
          else
          {
            IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
              shape_loc = scafe->GetShape(mip.IP(), lh);
              double fac = mip.GetWeight();
              elmat_loc += (fac*coef) * shape_loc * Trans(shape_loc);
            }
          }
        }
      } // if xfe
    } // loop over domain types

    // std::cout << " elmat_neg = " << elmat_neg << std::endl;
    // std::cout << " elmat_pos = " << elmat_pos << std::endl;
    elmat.Rows(0,ndof_neg).Cols(0,ndof_neg) = elmat_neg;
    elmat.Rows(ndof_neg,ndof_total).Cols(ndof_neg,ndof_total) = elmat_pos;
    // std::cout << " elmat = " << elmat << std::endl;
    // getchar();
  }

  template class FictXMassIntegrator<2>;
  template class FictXMassIntegrator<3>;

  static RegisterBilinearFormIntegrator<FictXMassIntegrator<2> > initfictxmass2d ("fictxmass", 2, 2);
  static RegisterBilinearFormIntegrator<FictXMassIntegrator<3> > initfictxmass3d ("fictxmass", 3, 2);

  template<int D>
  void FictXLaplaceIntegrator<D> ::
  CalcElementMatrix (const FiniteElement & base_fel,
                     const ElementTransformation & eltrans, 
                     FlatMatrix<double> & elmat,
                     LocalHeap & lh) const
  {
    static Timer timer ("FictXLaplaceIntegrator::CalcElementMatrix");
    RegionTimer reg (timer);

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const XFiniteElement * xfe[2]; xfe[0] = NULL; xfe[1] = NULL;
    const XDummyFE * dummfe[2]; dummfe[0] = NULL; dummfe[1] = NULL;

    const ScalarFiniteElement<D> * scafe = NULL;

    xfe[0] = dynamic_cast<const XFiniteElement* >(&cfel[0]);
    xfe[1] = dynamic_cast<const XFiniteElement* >(&cfel[1]);
    dummfe[0] = dynamic_cast<const XDummyFE* >(&cfel[0]);
    dummfe[1] = dynamic_cast<const XDummyFE* >(&cfel[1]);

    if (xfe[0] != NULL)
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[0]->GetBaseFE()));
    else
      scafe = &(dynamic_cast<const ScalarFiniteElement<D>& >(xfe[1]->GetBaseFE()));
    
    elmat = 0.0;

    if (!xfe && !dummfe) 
      throw Exception(" not containing X-elements?");

    // if(xfe[0])
    //   std::cout << " xfe[NEG] " << std::endl;
    // if(xfe[1])
    //   std::cout << " xfe[POS] " << std::endl;

    int ndof_neg = xfe[0]!=NULL ? scafe->GetNDof() : 0;
    int ndof_pos = xfe[1]!=NULL ? scafe->GetNDof() : 0;
    int ndof_total = ndof_neg+ndof_pos;
    FlatMatrixFixWidth<D> dshape_total(ndof_total,lh);
    FlatMatrixFixWidth<D> dshape_neg(ndof_neg,&dshape_total(0,0));
    FlatMatrixFixWidth<D> dshape_pos(ndof_pos,&dshape_total(ndof_neg,0));

    dshape_total = 0.0;

    FlatMatrix<> elmat_neg(ndof_neg,ndof_neg,lh);
    FlatMatrix<> elmat_pos(ndof_pos,ndof_pos,lh);
    
    elmat_neg = 0.0;
    elmat_pos = 0.0;

    int p = scafe->Order();
    DOMAIN_TYPE dt = POS;
    for (dt=POS; dt<IF; dt=(DOMAIN_TYPE)((int)dt+1))
    {
      // std::cout << " dt = " << dt << std::endl;
      FlatMatrixFixWidth<D> & dshape_loc = dt == NEG ? dshape_neg : dshape_pos;
      FlatMatrix<> & elmat_loc = dt == NEG ? elmat_neg : elmat_pos;

      const int idxfe = dt==POS ? 1 : 0;

      if( xfe[idxfe]) 
      {
        // if (! xfe[dt]->Empty())
        {
          const FlatXLocalGeometryInformation & xgeom(xfe[idxfe]->GetFlatLocalGeometry());
          const FlatCompositeQuadratureRule<D> & fcompr(xgeom.GetCompositeRule<D>());
          const FlatQuadratureRule<D> & fquad(fcompr.GetRule(dt));

          if (fquad.Size() > 0)
          {
            for (int i = 0; i < fquad.Size(); ++i)
            {
            
              IntegrationPoint ip(&fquad.points(i,0),fquad.weights(i));
              MappedIntegrationPoint<D,D> mip(ip, eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
          
              scafe->CalcMappedDShape(mip, dshape_loc);
              double fac = mip.GetWeight();
              elmat_loc += (fac*coef) * dshape_loc * Trans(dshape_loc);

            } // quad rule
          }
          else
          {
            IntegrationRule pir = SelectIntegrationRule (eltrans.GetElementType(), 2*p);
            for (int i = 0 ; i < pir.GetNIP(); i++)
            {
              MappedIntegrationPoint<D,D> mip(pir[i], eltrans);
              double coef = dt == POS ? coef_pos->Evaluate(mip) : coef_neg->Evaluate(mip);
              scafe->CalcMappedDShape(mip, dshape_loc);
              double fac = mip.GetWeight();
              elmat_loc += (fac*coef) * dshape_loc * Trans(dshape_loc);
            }
          }
        }
      } // if xfe
    } // loop over domain types

    // std::cout << " elmat_neg = " << elmat_neg << std::endl;
    // std::cout << " elmat_pos = " << elmat_pos << std::endl;
    elmat.Rows(0,ndof_neg).Cols(0,ndof_neg) = elmat_neg;
    elmat.Rows(ndof_neg,ndof_total).Cols(ndof_neg,ndof_total) = elmat_pos;
    // std::cout << " elmat = " << elmat << std::endl;
    // getchar();
  }

  template class FictXLaplaceIntegrator<2>;
  template class FictXLaplaceIntegrator<3>;

  static RegisterBilinearFormIntegrator<FictXLaplaceIntegrator<2> > initfictxlaplace2d ("fictxlaplace", 2, 2);
  static RegisterBilinearFormIntegrator<FictXLaplaceIntegrator<3> > initfictxlaplace3d ("fictxlaplace", 3, 2);


}

/// coefficientfunction statt function-pointer
