#include <variant>
#include "../cutint/cutintegral.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"

shared_ptr<BilinearFormIntegrator> CutIntegral :: MakeBilinearFormIntegrator()
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                    {
                      if (dynamic_cast<ProxyFunction*> (&cf))
                        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                          has_other = true;
                    });
  if (has_other && (cdx->element_vb != BND) && !cdx->skeleton)
    throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

  shared_ptr<BilinearFormIntegrator> bfi;
  if (!has_other && !dx.skeleton)
    bfi  = make_shared<SymbolicCutBilinearFormIntegrator> (*(cdx->lsetintdom), cf, dx.vb, dx.element_vb);
  else
  {
    if (cdx->lsetintdom->GetTimeIntegrationOrder() >= 0)
      throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for time_order >= 0..");
    if (cdx->vb == BND)
      throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for boundaries..");
    bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (*(cdx->lsetintdom), cf);
  }

  if (cdx->definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*cdx->definedon); definedon_bitarray)
        bfi->SetDefinedOn(*definedon_bitarray);
      /*
        // can't do that withouyt mesh
      if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
        {
          Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
          bfi->SetDefinedOn(reg.Mask());
        }
      */
    }
  bfi->SetDeformation(cdx->deformation);               
  bfi->SetBonusIntegrationOrder(cdx->bonus_intorder);
  if(cdx->definedonelements)
    bfi->SetDefinedOnElements(cdx->definedonelements);
  // for (auto both : dx.userdefined_intrules)
  //   bfi->SetIntegrationRule(both.first, *both.second);

  return bfi;
}



shared_ptr<LinearFormIntegrator> CutIntegral :: MakeLinearFormIntegrator()
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                    {
                      if (dynamic_cast<ProxyFunction*> (&cf))
                        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                          has_other = true;
                    });
  if (has_other && (cdx->element_vb != BND) && !cdx->skeleton)
    throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

  shared_ptr<LinearFormIntegrator> lfi;
  lfi  = make_shared<SymbolicCutLinearFormIntegrator> (*(cdx->lsetintdom), cf, dx.vb);
  if (cdx->definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*cdx->definedon); definedon_bitarray)
        lfi->SetDefinedOn(*definedon_bitarray);
      /*
        // can't do that withouyt mesh
      if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
        {
          Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
          lfi->SetDefinedOn(reg.Mask());
        }
      */
    }
  lfi->SetDeformation(cdx->deformation);               
  lfi->SetBonusIntegrationOrder(cdx->bonus_intorder);
  if(cdx->definedonelements)
    lfi->SetDefinedOnElements(cdx->definedonelements);
  // for (auto both : dx.userdefined_intrules)
  //   lfi->SetIntegrationRule(both.first, *both.second);

  return lfi;
}


template <typename TSCAL>
TSCAL CutIntegral :: T_CutIntegrate (const ngcomp::MeshAccess & ma,
                                  FlatVector<TSCAL> element_wise)
{
  LocalHeap glh(10000000, "lh-T_CutIntegrate");
  bool space_time = cdx->lsetintdom->GetTimeIntegrationOrder() >= 0;
  if (cdx->element_vb == BND)
    throw Exception("CutIntegrate can only deal with VOL a.t.m..");
  TSCAL sum = 0.0;

  BitArray defon;

  if (cdx->definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*cdx->definedon))
        defon = *definedon_bitarray;
      if (auto definedon_string = get_if<string> (&*cdx->definedon))
        {
          shared_ptr<MeshAccess> spma(const_cast<MeshAccess*>(&ma), NOOP_Deleter);
          Region reg(spma, dx.vb, *definedon_string);
          defon = reg.Mask();
        }
    }
  
  int DIM = ma.GetDimension();
  int cfdim = cf->Dimension();
  if(cfdim != 1)
    throw Exception("only implemented for 1 dimensional coefficientfunctions");

  ma.IterateElements(VOL, glh, [&] (Ngs_Element el, LocalHeap & lh)
  {
    if (defon.Size() && !defon.Test(el.GetIndex()))
      return;
    if (cdx->definedonelements && !cdx->definedonelements->Test(el.Nr()))
      return;

    auto & trafo1 = ma.GetTrafo (el, lh);
    auto & trafo = trafo1.AddDeformation(this->cdx->deformation.get(), lh);
    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(*cdx->lsetintdom,trafo,lh);

    if (ir != nullptr)
    {
      BaseMappedIntegrationRule & mir = trafo(*ir, lh);
      FlatMatrix<TSCAL> val(mir.Size(), 1, lh);

      cf -> Evaluate (mir, val);

      TSCAL lsum(0.0);
      for (int i = 0; i < mir.Size(); i++)
          lsum += mir[i].GetMeasure()*wei_arr[i]*val(i,0);

      if (element_wise.Size())
        element_wise(el.Nr()) += lsum;
      
      AtomicAdd(sum,lsum);
    }
  });
  return sum;
}


double CutIntegral::Integrate (const ngcomp::MeshAccess & ma,
                            FlatVector<double> element_wise)
{ 
  return T_CutIntegrate(ma, element_wise);
}

Complex CutIntegral::Integrate (const ngcomp::MeshAccess & ma,
                              FlatVector<Complex> element_wise)
{ return T_CutIntegrate(ma, element_wise);}


template double CutIntegral :: T_CutIntegrate<double> (const ngcomp::MeshAccess & ma,
                                                  FlatVector<double> element_wise);
template Complex CutIntegral :: T_CutIntegrate<Complex> (const ngcomp::MeshAccess & ma,
                                                    FlatVector<Complex> element_wise);
