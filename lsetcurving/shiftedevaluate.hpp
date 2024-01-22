#ifndef FILE_SHIFTEDEVALUATE_HPP
#define FILE_SHIFTEDEVALUATE_HPP

#include <fem.hpp>   // for ScalarFiniteElement
#include <ngstd.hpp> // for Array
#include <comp.hpp>

using namespace ngcomp;
namespace ngfem
{

  template <int SpaceD>
  class DiffOpShiftedEval : public DifferentialOperator
  {
private:
      shared_ptr<DifferentialOperator> evaluator;

  public:

    shared_ptr<GridFunction> back;
    shared_ptr<GridFunction> forth;

    //enum { DIM = 1 };          // 1 copies of the spaces
    //enum { DIM_SPACE = SpaceD };    // D-dim space
    //enum { DIM_ELEMENT = SpaceD };  // D-dim elements (in contrast to boundary elements)
    //enum { DIM_DMAT = 1 };     // 1-matrix
    //enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)

    DiffOpShiftedEval(shared_ptr<GridFunction> aback,shared_ptr<GridFunction> aforth, shared_ptr<DifferentialOperator> a_evaluator)
      : DifferentialOperator(a_evaluator->Dim(), a_evaluator->BlockDim(), VOL, a_evaluator->DiffOrder()),
        evaluator(a_evaluator), back(aback), forth(aforth)
    {
      SetDimensions(Array<int> ( { a_evaluator->Dim()} ));
    }
    /*
    virtual int Dim() const { return DIM_DMAT; }
    virtual bool Boundary() const { return int(DIM_SPACE) > int(DIM_ELEMENT); }
    virtual int DiffOrder() const { return DIFFORDER; }
    */
    virtual string Name() const { return "shifted_evaluation"; }

    virtual bool operator== (const DifferentialOperator & diffop2) const
    { return typeid(*this) == typeid(diffop2); }


    virtual void
    CalcMatrix (const FiniteElement & bfel,
        const BaseMappedIntegrationPoint & bmip,
        BareSliceMatrix<double,ColMajor> mat,
        LocalHeap & lh) const;


    virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<double> x, 
           FlatVector<double> flux,
           LocalHeap & lh) const;
    
    virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x,
                LocalHeap & lh) const;

  };



}
#endif
