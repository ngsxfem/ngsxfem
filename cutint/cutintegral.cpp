#include <variant>
#include "../cutint/cutintegral.hpp"
#include "../xfem/symboliccutbfi.hpp"
#include "../xfem/symboliccutlfi.hpp"

CutIntegral :: CutIntegral (shared_ptr<CoefficientFunction> _cf, shared_ptr<CutDifferentialSymbol> _dx)
  : Integral(_cf, *_dx), lsetintdom(_dx->lsetintdom) { ; }

shared_ptr<BilinearFormIntegrator> CutIntegral :: MakeBilinearFormIntegrator() const
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                    {
                      if (dynamic_cast<ProxyFunction*> (&cf))
                        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                          has_other = true;
                    });
  if (has_other && (dx.element_vb != BND) && !dx.skeleton)
    throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

  shared_ptr<BilinearFormIntegrator> bfi;
  if (!has_other && !dx.skeleton)
    bfi  = make_shared<SymbolicCutBilinearFormIntegrator> (*lsetintdom, cf, dx.vb, dx.element_vb);
  else
  {
    //if (lsetintdom->GetTimeIntegrationOrder() >= 0)
    //  throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for time_order >= 0..");
    if (dx.vb == BND)
      throw Exception("Symbolic cuts on facets and boundary not yet (implemented/tested) for boundaries..");
    bfi = make_shared<SymbolicCutFacetBilinearFormIntegrator> (*lsetintdom, cf);
  }

  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
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
  bfi->SetDeformation(dx.deformation);               
  bfi->SetBonusIntegrationOrder(dx.bonus_intorder);
  if(dx.definedonelements)
    bfi->SetDefinedOnElements(dx.definedonelements);
  // for (auto both : dx.userdefined_intrules)
  //   bfi->SetIntegrationRule(both.first, *both.second);

  return bfi;
}



shared_ptr<LinearFormIntegrator> CutIntegral :: MakeLinearFormIntegrator() const
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                    {
                      if (dynamic_cast<ProxyFunction*> (&cf))
                        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                          has_other = true;
                    });
  if (has_other && (dx.element_vb != BND) && !dx.skeleton)
    throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

  shared_ptr<LinearFormIntegrator> lfi;
  lfi  = make_shared<SymbolicCutLinearFormIntegrator> (*lsetintdom, cf, dx.vb);
  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
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
  lfi->SetDeformation(dx.deformation);               
  lfi->SetBonusIntegrationOrder(dx.bonus_intorder);
  if(dx.definedonelements)
    lfi->SetDefinedOnElements(dx.definedonelements);
  // for (auto both : dx.userdefined_intrules)
  //   lfi->SetIntegrationRule(both.first, *both.second);

  return lfi;
}


template <typename TSCAL>
TSCAL CutIntegral :: T_CutIntegrate (const ngcomp::MeshAccess & ma,
                                  FlatVector<TSCAL> element_wise)
{
  static Timer timer("CutIntegral::T_CutIntegrate");
  RegionTimer reg (timer);
  LocalHeap glh(1000000000, "lh-T_CutIntegrate");
  // bool space_time = lsetintdom->GetTimeIntegrationOrder() >= 0;
  if (dx.element_vb == BND)
    throw Exception("CutIntegrate can only deal with VOL a.t.m..");

  BitArray defon;

  if (dx.definedon)
    {
      if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon))
        defon = *definedon_bitarray;
      if (auto definedon_string = get_if<string> (&*dx.definedon))
        {
          shared_ptr<MeshAccess> spma(const_cast<MeshAccess*>(&ma), NOOP_Deleter);
          Region reg(spma, dx.vb, *definedon_string);
          defon = reg.Mask();
        }
    }
  
  // int DIM = ma.GetDimension();
  int cfdim = cf->Dimension();
  if(cfdim != 1)
    throw Exception("only implemented for 1 dimensional coefficientfunctions");

	if (globxvar.SIMD_EVAL) {
		try {
      TSCAL sum = 0.0;
		  ma.IterateElements(VOL, glh, [&] (Ngs_Element el, LocalHeap & lh)
		  {
		  	if (defon.Size() && !defon.Test(el.GetIndex()))
		  		return;
		  	if (dx.definedonelements && !dx.definedonelements->Test(el.Nr()))
		  		return;

		  	auto & trafo1 = ma.GetTrafo (el, lh);
		  	auto & trafo = trafo1.AddDeformation(this->dx.deformation.get(), lh);
		  	const IntegrationRule *ns_ir;
		  	Array<double> ns_wei_arr;
		  	tie (ns_ir, ns_wei_arr) = CreateCutIntegrationRule(*lsetintdom,trafo,lh);
		  	if (ns_ir == nullptr)
		  		return;
  
		  	SIMD_IntegrationRule simd_ir(*ns_ir, lh);
		  	FlatArray<SIMD<double>> simd_wei_arr = CreateSIMD_FlatArray(ns_wei_arr, lh);

		  	//if (&simd_ir != nullptr)
		  	{
		  		SIMD_BaseMappedIntegrationRule & simd_mir = trafo(simd_ir, lh);
		  		FlatMatrix<SIMD<TSCAL>> val(simd_mir.Size(), 1, lh);
  
		  		cf -> Evaluate (simd_mir, val);
  
		  		SIMD<TSCAL> lsum(0.0);
		  		for (int i = 0; i < simd_mir.Size(); i++)
		  				lsum += simd_mir[i].GetMeasure()*simd_wei_arr[i]*val(i,0);
  
		  		if (element_wise.Size())
		  			element_wise(el.Nr()) += HSum(lsum); // problem?
  
		  		AtomicAdd(sum, HSum(lsum));
		  	}
		  });
		  return ma.GetCommunicator().AllReduce(sum, NG_MPI_SUM);
    } catch (ExceptionNOSIMD e) {
      cout << IM(6) << e.What()
           << "switching to non-SIMD evaluation" << endl;
    }
	}
	
  TSCAL sum = 0.0;
  ma.IterateElements(VOL, glh, [&] (Ngs_Element el, LocalHeap & lh)
  {
    if (defon.Size() && !defon.Test(el.GetIndex()))
      return;
    if (dx.definedonelements && !dx.definedonelements->Test(el.Nr()))
      return;

    auto & trafo1 = ma.GetTrafo (el, lh);
    auto & trafo = trafo1.AddDeformation(this->dx.deformation.get(), lh);
    const IntegrationRule * ir;
    Array<double> wei_arr;
    tie (ir, wei_arr) = CreateCutIntegrationRule(*lsetintdom,trafo,lh);

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
  return ma.GetCommunicator().AllReduce(sum, NG_MPI_SUM);
}


double CutIntegral::Integrate (const ngcomp::MeshAccess & ma,
                            FlatVector<double> element_wise)
{ 
  return T_CutIntegrate(ma, element_wise);
}

Complex CutIntegral::Integrate (const ngcomp::MeshAccess & ma,
                              FlatVector<Complex> element_wise)
{ return T_CutIntegrate(ma, element_wise);}

FacetPatchIntegral::FacetPatchIntegral (shared_ptr<CoefficientFunction> _cf,
                                        shared_ptr<FacetPatchDifferentialSymbol> _dx)
      : Integral(_cf, *_dx), time_order(_dx->time_order), tref(_dx->tref) { ; }


shared_ptr<BilinearFormIntegrator> FacetPatchIntegral :: MakeBilinearFormIntegrator() const
{
  // check for DG terms
  bool has_other = false;
  cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                    {
                      if (dynamic_cast<ProxyFunction*> (&cf))
                        if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                          has_other = true;
                    });
  if (!has_other)
    cout << IM(3) << " no Other() used?!" << endl;

  auto bfi = make_shared<SymbolicFacetPatchBilinearFormIntegrator> (cf);
  if ((tref) && (time_order > -1))
    throw Exception("not reference time fixing for space-time integration domain");
  bfi->SetTimeIntegrationOrder(time_order);
  if (tref)
    bfi->SetReferenceTime(*tref);

  if (dx.definedon)
  {
    if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
      bfi->SetDefinedOn(*definedon_bitarray);
  }
  bfi->SetDeformation(dx.deformation);               
  bfi->SetBonusIntegrationOrder(dx.bonus_intorder);
  if(dx.definedonelements)
    bfi->SetDefinedOnElements(dx.definedonelements);
  return bfi;
}



template double CutIntegral :: T_CutIntegrate<double> (const ngcomp::MeshAccess & ma,
                                                  FlatVector<double> element_wise);
template Complex CutIntegral :: T_CutIntegrate<Complex> (const ngcomp::MeshAccess & ma,
                                                    FlatVector<Complex> element_wise);
