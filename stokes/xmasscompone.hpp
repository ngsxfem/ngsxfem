#ifndef FILE_XMASSCOMPONE_HPP
#define FILE_XMASSCOMPONE_HPP

/// from ngsolve
#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array

#include "../cutint/xintegration.hpp"
#include "../xfem/xfiniteelement.hpp"
//#include "../spacetime/spacetimeintegrators.hpp"
//#include "../utils/stcoeff.hpp"

namespace ngfem
{
  template<int D>
  void CastXScalarFiniteElements (const FiniteElement & base_fel,
                                  const ScalarFiniteElement<D> * & scafe,
                                  const XFiniteElement * & xfe,
                                  const XDummyFE * & dummfe); 
  template <int D>
  class xMassCompOne : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_neg;
    shared_ptr<CoefficientFunction> coef_pos;
  public:
    xMassCompOne (const Array<shared_ptr<CoefficientFunction>> & coeffs)
      : coef_neg(coeffs[0]),coef_pos(coeffs[1]) { ; } 

    virtual ~xMassCompOne(){ ; };

    virtual string Name () const { return "Mass on first component"; }

    virtual int DimElement () const { return D; }
    virtual int DimSpace () const { return D; }
    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }


    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;

    virtual int GetDimension () const { return 1; }


    virtual void
    CalcFlux (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & bmip,
              FlatVector<double> elx, 
              FlatVector<double> flux,
              bool applyd,
              LocalHeap & lh) const
    {
      //diffop -> Apply (fel, bmip, elx, flux, lh);

      cout<<"dynamic casts!!!!"<<endl;
      const CompoundFiniteElement & cfel = dynamic_cast<const CompoundFiniteElement&> (fel);                                
      cout<<"compund check"<<endl;
      
      const ScalarFiniteElement<D> & scafe = dynamic_cast<const ScalarFiniteElement<D> &> (cfel[0]);
      
      cout<<"apply"<<endl;

      //const CompoundFiniteElement & cfela = 
      //  dynamic_cast<const CompoundFiniteElement&> (fel);
      //const CompoundFiniteElement & cfel = 
      //  dynamic_cast<const CompoundFiniteElement&> (cfela[0]);
      //const ScalarFiniteElement<D> & fel_u = 
      //  dynamic_cast<const ScalarFiniteElement<D>&> (cfel[0]);

      cout<<"elx = "<<elx<<endl;
      cout<<"flux = "<<flux<<endl;


      DiffOpId<D>::Apply(scafe,bmip,elx,flux,lh);

      cout<<"endapply"<<endl;
      // FlatVec<DMATOP::DIM_DMAT,double> hflux(&flux(0));

      //if (applyd)
      //dmatop.Apply1 (fel, bmip, hflux, lh);
    }

  };

}

#endif

