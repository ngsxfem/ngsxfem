/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
// #include "geom.hpp"
// #include "cuttriang.hpp"
// #include "../fem/xfiniteelement.hpp"

using namespace ngsolve;

typedef std::pair<double,double> TimeInterval;

// /// A space-time integration point 
// template <int DIMS, int DIMR, typename SCAL = double>
// class SpaceTimeMappedIntegrationPoint : public MappedIntegrationPoint<DIMS,DIMR,SCAL>
// {
// protected:
//     /// reference time 
//     double tref; 
//     /// time interval
//     TimeInterval  ti;
// public:
//     ///
//     SpaceTimeMappedIntegrationPoint (const IntegrationPoint & aip,
//                                      const ElementTransformation & aeltrans,
//                                      const double & atref,
//                                      const TimeInterval & ati)
//         : MappedIntegrationPoint<DIMS,DIMR,SCAL>(aip,aeltrans), tref(atref), ti(ati)  { ; }
//     ///
//     double GetTime() const { return (1.0-tref) * ti.first + tref * ti.second; }
// };

template <ELEMENT_TYPE ET, int DIM>
ScalarFiniteElement<DIM> * GetSpaceFE (const Ngs_Element & ngel, int order, LocalHeap & lh)
{
    L2HighOrderFE<ET> * hofe =  new (lh) L2HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.vertices);
    hofe -> SetOrder (order);
    hofe -> ComputeNDof();
    
    return hofe;
}

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
template<int D>
class NumProcTestXFEM : public NumProc
{
protected:

    // CoefficientFunction levelset function
    CoefficientFunction *coef_lset;
    
    // is space time?
    bool isspacetime;

    double * time;

    int order_space;
    int order_time;

    TimeInterval ti;
public:
    
    NumProcTestXFEM (PDE & apde, const Flags & flags)
        : NumProc (apde)
    {
        cout << " \n\nNumProcTestXFEM - constructor start \n\n " << endl;
        string lsetstr = flags.GetStringFlag("levelset","none");
        isspacetime = flags.GetDefineFlag("spacetime");
        order_time = (int) flags.GetNumFlag("order_time",1);
        order_space = (int) flags.GetNumFlag("order_space",1);

        Array<double> ti1;
        ti1 = flags.GetNumListFlag ("timeinterval");
        if (ti1.Size() >= 2)
        {
            ti.first = ti1[0]; ti.second = ti1[1];
        }
        else if (isspacetime)
            throw Exception("please describe time interval");


        string timestring = flags.GetStringFlag("time","t");
        time = &(pde.GetVariable(timestring,1));
        if (isspacetime)
            cout << " is space time " << endl;
        else 
            cout << " is space only " << endl;

        coef_lset = pde.GetCoefficientFunction(lsetstr,0);
        if (coef_lset == NULL)
            throw Exception("NumProcTestXFEM: levelset not given");
        cout << " \n\nNumProcTestXFEM - constructor end \n\n " << endl;
    }
  
    virtual string GetClassName () const
    {
        return "NumProcTestXFEM";
    }
  

    virtual void Do (LocalHeap & lh)
    {
        HeapReset hr(lh);

        ELEMENT_TYPE et_time;
        DGFiniteElement<1> * fel_time;
        if(isspacetime)
        {
            et_time = ET_SEGM;
            fel_time = new (lh) L2HighOrderFE<ET_SEGM>  (order_time);
        }
        else
        {
            et_time = ET_POINT;
            fel_time = NULL;
        }

        for (int elnr = 0; elnr < ma.GetNE(); ++elnr)
        {
            HeapReset hr(lh);
            Ngs_Element ngel = ma.GetElement(elnr);

            ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), lh);
            ELEMENT_TYPE et_space = eltrans.GetElementType();

            // project to a local respresentation of the level set function
            
            DGFiniteElement<D> * fel_space;
            switch (et_space)
            {
            case ET_TRIG: 
                fel_space = dynamic_cast<DGFiniteElement<D> * >(GetSpaceFE<ET_TRIG,2>(ngel, order_space, lh)); 
                break;
            case ET_QUAD: 
                fel_space = dynamic_cast<DGFiniteElement<D> * >(GetSpaceFE<ET_QUAD,2>(ngel, order_space, lh)); 
                break;
            case ET_TET:
                fel_space = dynamic_cast<DGFiniteElement<D> * >(GetSpaceFE<ET_TET,3>(ngel, order_space, lh)); 
                break;
            default:
                cout << "eltype is " << et_space << endl;
                throw Exception ("illegal element for  GetSpaceFE");
            }

            int ndof_space = fel_space->GetNDof();
            int ndof_time = isspacetime ? fel_time->GetNDof() : 1;
            int ndof_total = ndof_space * ndof_time;

            // cout << " fel_space ndofs = " << ndof_space << endl;
            // cout << " fel_time ndofs = " << ndof_time << endl;
            // cout << " total ndofs = " << ndof_total << endl;

            const IntegrationRule & ir_space = SelectIntegrationRule (et_space, 2*order_space);
            const IntegrationRule & ir_time = SelectIntegrationRule (et_time, 2*order_time);
            // const IntegrationRule & ir1d = SelectIntegrationRule (etfacet, GetIntegrationOrderFacets(order));

            FlatVector<> massdiagspace(ndof_space,lh);
            fel_space->GetDiagMassMatrix(massdiagspace);

            FlatVector<> massdiagtime(ndof_time,lh);
            if (isspacetime)
                fel_time->GetDiagMassMatrix(massdiagtime);
            else
                massdiagtime = 1.0;


            FlatVector<> massdiag(ndof_total,lh);
            
            for (int m = 0; m < ndof_time; m++)
                for (int n = 0; n < ndof_space; n++)
                    massdiag(ndof_space*m+n) = (massdiagtime(m) * massdiagspace(n));


            FlatVector<> rhs(ndof_total,lh);
            rhs = 0.0;
            FlatVector<> shape_space(ndof_space,lh);
            FlatVector<> shape_time(ndof_time,lh);

            for (int l = 0; l < ir_time.GetNIP(); l++)
            {
                double current = ti.first +  ir_time[l](0) * (ti.second - ti.first);
                *time = current;
                if (isspacetime)
                    fel_time->CalcShape(ir_time[l],shape_time);
                else 
                    shape_time = 1.0;
                for (int k = 0; k < ir_space.GetNIP(); k++)
                {
                    double val = 0.0;
                    MappedIntegrationPoint<D,D> mip (ir_space[k], eltrans);
                    val = coef_lset->Evaluate(mip);
                    fel_space->CalcShape(ir_space[k],shape_space);

                    const double fac = ir_space[k].Weight() * (ti.second-ti.first);

                    for (int m = 0; m < ndof_time; m++)
                        for (int n = 0; n < ndof_space; n++)
                            rhs(ndof_space*m+n) += shape_time(m) * shape_space(n) * val * fac;
                }
            }

            FlatVector<> lset_loc(ndof_total,lh);
            
            for (int m = 0; m < ndof_total; m++)
                lset_loc(m) = rhs(m) / massdiag(m);
            
            // make CompositeRule based :
            // - on et_space and et_time:
            //    * 2D   stationary: ET_TRIG + ET_POINT
            //    * 2D instationary: ET_TRIG + ET_SEGM
            //    * 3D   stationary: ET_TET  + ET_POINT
            //    * 3D instationary: ET_TET  + ET_SEGM
            // - levelset function as STFiniteElement:
            //    * TP - FE of fel_space and fel_time 
            //    * => Eval-Object

        }
    }
};

static RegisterNumProc<NumProcTestXFEM<2> > npinittestxfem2d("testxfem");
static RegisterNumProc<NumProcTestXFEM<3> > npinittestxfem3d("testxfem");
