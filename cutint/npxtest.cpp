/*********************************************************************/
/* File:   npxtest.cpp                                               */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   2. Feb. 2014                                              */
/*********************************************************************/


/*
 */

#include <solve.hpp>
#include "xintegration.hpp"
#include "../common/spacetimefespace.hpp"

using namespace ngsolve;
using namespace xintegration;

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

// template <ELEMENT_TYPE ET, int DIM>
// ScalarFiniteElement<DIM> * GetSpaceFE (const Ngs_Element & ngel, int order, LocalHeap & lh)
// {
//     L2HighOrderFE<ET> * hofe =  new (lh) L2HighOrderFE<ET> ();
    
//     hofe -> SetVertexNumbers (ngel.vertices);
//     hofe -> SetOrder (order);
//     hofe -> ComputeNDof();
    
//     return hofe;
// }

/* ---------------------------------------- 
   numproc
   ---------------------------------------- */
template<int D>
class NumProcTestXFEM : public NumProc
{
protected:
    bool isspacetime;
    
    int order_space;
    int order_time;

    int num_int_ref_space;
    int num_int_ref_time;

    bool output = false;
    bool bound = false;

    double vmax = 1.0;
    double dt = 1.0;

    GridFunction * gf_lset;
public:
    
    NumProcTestXFEM (PDE & apde, const Flags & flags)
        : NumProc (apde)
    {
        cout << " \n\nNumProcTestXFEM - constructor start \n\n " << endl;

		gf_lset = pde.GetGridFunction (flags.GetStringFlag ("levelset", "lset"));
        isspacetime = flags.GetDefineFlag("spacetime");
        output = flags.GetDefineFlag("output");
        bound = flags.GetDefineFlag("bound");
        num_int_ref_space = (int) flags.GetNumFlag("num_int_ref_space",0);
        num_int_ref_time = (int) flags.GetNumFlag("num_int_ref_time",0);

        order_space = (int) flags.GetNumFlag("int_order_space",-1);
        order_time = (int) flags.GetNumFlag("int_order_time",-1);

        dt = flags.GetNumFlag("dt",1.0);
        vmax = flags.GetNumFlag("vmax",1.0);

        const FESpace & fes = gf_lset->GetFESpace();
        if (isspacetime)
        {
            const SpaceTimeFESpace & fes_st = dynamic_cast< const SpaceTimeFESpace & > (gf_lset->GetFESpace());
            if (order_space == -1)
                order_space = 2*fes_st.OrderSpace();
            if (order_time == -1)
                order_time = 2*fes_st.OrderTime();
        }
        else
        {
            order_time = 0;
            if (order_space == -1)
                order_space = 2 * fes.GetOrder();
        }

        cout << " \n\nNumProcTestXFEM - constructor end \n\n " << endl;
    }
  
    ~NumProcTestXFEM()
    {
    }

    virtual string GetClassName () const
    {
        return "NumProcTestXFEM";
    }
  

    virtual void Do (LocalHeap & clh)
    {
        static Timer timerdo ("NumProcTestXFEM::Do");
        RegionTimer reg (timerdo);

        HeapReset hr(clh);

        Array<int> els_of_dt(3);
        els_of_dt = 0.0;

        Array<double> meas_of_dt(3);
        meas_of_dt = 0.0;

        Array<int> sels_of_dt(3);
        sels_of_dt = 0.0;

        Array<double> meas_of_dt_bnd(3);
        meas_of_dt_bnd = 0.0;

        ofstream outneg_st("negpoints_st.out");
        ofstream outif_st("ifpoints_st.out");
        ofstream outneg_s("negpoints_s.out");
        ofstream outif_s("ifpoints_s.out");

#pragma omp parallel
        {
            LocalHeap lh(clh.Split());
#pragma omp for schedule(static)
            for (int elnr = 0; elnr < ma.GetNE(); ++elnr)
            {
                // Ngs_Element ngel = ma.GetElement(elnr);
                HeapReset hr(lh);

                ElementTransformation & eltrans = ma.GetTrafo (ElementId(VOL,elnr), lh);
                ELEMENT_TYPE et_space = eltrans.GetElementType();
                ELEMENT_TYPE et_time = isspacetime ? ET_SEGM : ET_POINT;
            
                IntegrationPoint ip(0.0);
                MappedIntegrationPoint<D,D> mip(ip,eltrans);
                const double absdet = mip.GetJacobiDet();

                const FESpace & fes = gf_lset->GetFESpace();
                const FiniteElement & fel = fes.GetFE(elnr, lh);
                Array<int> dnums;
                fes.GetDofNrs(elnr,dnums);
            
                FlatVector<> linvec(dnums.Size(),lh);
                gf_lset->GetVector().GetIndirect(dnums,linvec);

                if( et_space == ET_TRIG && et_time == ET_SEGM)
                {
                    const ScalarSpaceTimeFiniteElement<2> &  fel_st 
                        = dynamic_cast<const ScalarSpaceTimeFiniteElement<2> & >(fel);
                    ScalarFEEvaluator<2> lset_eval(fel_st, linvec, lh);
                    PointContainer<3> pc;
                    CompositeQuadratureRule<3> compositerule;
                    NumericalIntegrationStrategy<ET_TRIG,ET_SEGM> numint(lset_eval, pc,
                                                                         compositerule, 
                                                                         lh,
                                                                         order_space, order_time,
                                                                         num_int_ref_space, 
                                                                         num_int_ref_time);
                    numint.SetVerticesSpace();
                    numint.SetVerticesTime();
                    numint.SetDistanceThreshold(2.0*sqrt(absdet)+vmax/dt);
                    {
                        static Timer timer ("npxtest - MakeQuadRule - total");
                        RegionTimer reg (timer);
                        els_of_dt[numint.MakeQuadRule()]++;
                    }
#pragma omp critical(measpos)
                    for (int i = 0; i < compositerule.quadrule_pos.Size(); ++i)
                        meas_of_dt[POS] += absdet * compositerule.quadrule_pos.weights[i];
#pragma omp critical(measneg)
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        meas_of_dt[NEG] += absdet * compositerule.quadrule_neg.weights[i];
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                    {
                        static Timer timer ("npxtest - Transform n");
                        RegionTimer reg (timer);
                        Vec<3> st_point = compositerule.quadrule_if.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1));
                        MappedIntegrationPoint<D,D> mip(ip,eltrans);
                        Mat<2,2> Finv = mip.GetJacobianInverse();
                        Vec<3> nref = compositerule.quadrule_if.normals[i];
                        Vec<2> nref_sp (nref(0),nref(1));
                        Vec<2> n_sp = absdet * Trans(Finv) * nref_sp ;
                        double n_t = nref(2) * absdet;
                        Vec<3> n (n_sp(0),n_sp(1),n_t);
                        double fac = L2Norm(n);
#pragma omp critical(measif)
                        meas_of_dt[IF] +=  compositerule.quadrule_if.weights[i] * fac;
                    }
                    if (output)
                    {
                        for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        {
                            static Timer timer ("npxtest - output neg");
                            RegionTimer reg (timer);
                            Vec<3> st_point = compositerule.quadrule_neg.points[i];
                            IntegrationPoint ip(st_point(0),st_point(1));
                            MappedIntegrationPoint<D,D> mip(ip,eltrans);
                            Vec<2> mapped_point = mip.GetPoint();
                            outneg_st << mapped_point(0) << "\t" << mapped_point(1) << "\t" << st_point(2) << endl;
                        }
                        outneg_st << endl;
                        outneg_st << endl;
                        for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                        {
                            static Timer timer ("npxtest - output if");
                            RegionTimer reg (timer);
                            Vec<3> st_point = compositerule.quadrule_if.points[i];
                            IntegrationPoint ip(st_point(0),st_point(1));
                            MappedIntegrationPoint<D,D> mip(ip,eltrans);
                            Vec<2> mapped_point = mip.GetPoint();
                            outif_st << mapped_point(0) << "\t" << mapped_point(1) << "\t" << st_point(2) << endl;
                        }
                        outif_st << endl;
                        outif_st << endl;
                    }
                }

                if( et_space == ET_TRIG && et_time == ET_POINT)
                {
                    const ScalarSpaceTimeFiniteElement<2> *  scal_st_fel
                        = dynamic_cast<const ScalarSpaceTimeFiniteElement<2> * >(&fel);

                    ScalarFEEvaluator<2> lset_eval(fel, linvec, lh);
                    if (scal_st_fel != NULL)
                        lset_eval.FixTime(0.0);

                    PointContainer<2> pc;
                    CompositeQuadratureRule<2> compositerule;
                    NumericalIntegrationStrategy<ET_TRIG,ET_POINT> numint(lset_eval, pc,
                                                                          compositerule, 
                                                                          lh,
                                                                          order_space, order_time,
                                                                          num_int_ref_space, 
                                                                          num_int_ref_time);
                    numint.SetVerticesSpace();
                    numint.SetVerticesTime();
                    // numint.SetDistanceThreshold(0.125);

                    {
                        static Timer timer ("MakeQuadRule - total");
                        RegionTimer reg (timer);
                        els_of_dt[numint.MakeQuadRule()]++;
                    }
                    for (int i = 0; i < compositerule.quadrule_pos.Size(); ++i)
                        meas_of_dt[POS] += absdet * compositerule.quadrule_pos.weights[i];
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        meas_of_dt[NEG] += absdet * compositerule.quadrule_neg.weights[i];
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                    {
                        Vec<2> st_point = compositerule.quadrule_if.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1));
                        MappedIntegrationPoint<D,D> mip(ip,eltrans);
                        Mat<2,2> Finv = mip.GetJacobianInverse();
                        Vec<2> nref = compositerule.quadrule_if.normals[i];
                        Vec<2> n = absdet * Trans(Finv) * nref ;
                        double fac = L2Norm(n);
                        meas_of_dt[IF] += compositerule.quadrule_if.weights[i] * fac;
                    }
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                    {
                        Vec<2> st_point = compositerule.quadrule_neg.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1));
                        MappedIntegrationPoint<2,2> mip(ip,eltrans);
                        Vec<2> mapped_point = mip.GetPoint();
                        outneg_s << mapped_point(0) << "\t" << mapped_point(1) << endl;
                    }
                    outneg_s << endl;
                    outneg_s << endl;
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                    {
                        Vec<2> st_point = compositerule.quadrule_if.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1));
                        MappedIntegrationPoint<2,2> mip(ip,eltrans);
                        Vec<2> mapped_point = mip.GetPoint();
                        outif_s << mapped_point(0) << "\t" << mapped_point(1) << endl;
                    }
                    outif_s << endl;
                    outif_s << endl;
                }
            
                if ( et_space == ET_TET && et_time == ET_SEGM)
                {
                    const ScalarSpaceTimeFiniteElement<3> &  fel_st 
                        = dynamic_cast<const ScalarSpaceTimeFiniteElement<3> & >(fel);
                    ScalarFEEvaluator<3> lset_eval(fel_st, linvec, lh);

                    PointContainer<4> pc;
                    CompositeQuadratureRule<4> compositerule;
                    NumericalIntegrationStrategy<ET_TET,ET_SEGM> numint(lset_eval, pc,
                                                                        compositerule, 
                                                                        lh,
                                                                        order_space, order_time,
                                                                        num_int_ref_space, 
                                                                        num_int_ref_time);
                    numint.SetVerticesSpace();
                    numint.SetVerticesTime();

                    els_of_dt[numint.MakeQuadRule()]++;

                    for (int i = 0; i < compositerule.quadrule_pos.Size(); ++i)
                        meas_of_dt[POS] += compositerule.quadrule_pos.weights[i];
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        meas_of_dt[NEG] += compositerule.quadrule_neg.weights[i];
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                        meas_of_dt[IF] += compositerule.quadrule_if.weights[i];

                }

                if( et_space == ET_TET && et_time == ET_POINT)
                {
                    const ScalarSpaceTimeFiniteElement<3> *  scal_st_fel
                        = dynamic_cast<const ScalarSpaceTimeFiniteElement<3> * >(&fel);

                    ScalarFEEvaluator<3> lset_eval(fel, linvec, lh);
                    if (scal_st_fel != NULL)
                        lset_eval.FixTime(0.0);

                    PointContainer<3> pc;
                    CompositeQuadratureRule<3> compositerule;
                    NumericalIntegrationStrategy<ET_TET,ET_POINT> numint(lset_eval, pc,
                                                                         compositerule, 
                                                                         lh,
                                                                         order_space, order_time,
                                                                         num_int_ref_space, 
                                                                         num_int_ref_time);
                    numint.SetVerticesSpace();
                    numint.SetVerticesTime();
                    // numint.SetDistanceThreshold(0.125);

                    {
                        static Timer timer ("MakeQuadRule - total");
                        RegionTimer reg (timer);
                        els_of_dt[numint.MakeQuadRule()]++;
                    }

                    for (int i = 0; i < compositerule.quadrule_pos.Size(); ++i)
                        meas_of_dt[POS] += absdet * compositerule.quadrule_pos.weights[i];
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        meas_of_dt[NEG] += absdet * compositerule.quadrule_neg.weights[i];
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                    {
                        Vec<3> st_point = compositerule.quadrule_if.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1),st_point(2));
                        MappedIntegrationPoint<D,D> mip(ip,eltrans);
                        Mat<3,3> Finv = mip.GetJacobianInverse();
                        Vec<3> nref = compositerule.quadrule_if.normals[i];
                        Vec<3> n = absdet * Trans(Finv) * nref ;
                        double fac = L2Norm(n);
                        meas_of_dt[IF] += compositerule.quadrule_if.weights[i] * fac;
                    }
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                    {
                        Vec<3> st_point = compositerule.quadrule_neg.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1),st_point(2));
                        MappedIntegrationPoint<3,3> mip(ip,eltrans);
                        Vec<3> mapped_point = mip.GetPoint();
                        outneg_s << mapped_point(0) << "\t" << mapped_point(1) << "\t" << mapped_point(2) << endl;
                    }
                    outneg_s << endl;
                    outneg_s << endl;
                    for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
                    {
                        Vec<3> st_point = compositerule.quadrule_if.points[i];
                        IntegrationPoint ip(st_point(0),st_point(1),st_point(2));
                        MappedIntegrationPoint<D,D> mip(ip,eltrans);
                        Vec<3> mapped_point = mip.GetPoint();
                        outif_s << mapped_point(0) << "\t" << mapped_point(1) << "\t" << mapped_point(2) << endl;
                    }
                    outif_s << endl;
                    outif_s << endl;
                }

            
            }
        }
        cout << " pos elements : " << els_of_dt[POS] << endl;
        cout << " neg elements : " << els_of_dt[NEG] << endl;
        cout << " cut elements : " << els_of_dt[IF] << endl;
        cout << " pos measure : " << meas_of_dt[POS] << endl;
        cout << " neg measure : " << meas_of_dt[NEG] << endl;
        cout << " cut measure : " << meas_of_dt[IF] << endl;


        if (bound)
        {
            for (int selnr = 0; selnr < ma.GetNSE(); ++selnr)
            {
                LocalHeap lh(clh.Split());
                // Ngs_Element ngel = ma.GetElement(selnr);
                HeapReset hr(lh);

                ElementTransformation & seltrans = ma.GetTrafo (selnr, BND, lh);
                ELEMENT_TYPE et_space = seltrans.GetElementType();
                ELEMENT_TYPE et_time = isspacetime ? ET_SEGM : ET_POINT;
            
                IntegrationPoint ip(0.0);
                MappedIntegrationPoint<D-1,D> mip(ip,seltrans);
                const double absdet = mip.GetJacobiDet(); //??

                // std::cout << " mip.GitPoint() = " << mip.GetPoint() << std::endl;
                // IntegrationPoint ip2(1.0);
                // MappedIntegrationPoint<D-1,D> mip2(ip2,seltrans);
                // std::cout << " mip2.GitPoint() = " << mip2.GetPoint() << std::endl;


                const FESpace & fes = gf_lset->GetFESpace();
                const FiniteElement & fel = fes.GetSFE(selnr, lh);
                Array<int> dnums;
                fes.GetSDofNrs(selnr,dnums);
            
                // std::cout << " dnums = " << dnums << std::endl;
                FlatVector<> linvec(dnums.Size(),lh);
                gf_lset->GetVector().GetIndirect(dnums,linvec);

                // std::cout << " linvec = " << linvec << std::endl;

                if( et_space == ET_SEGM && et_time == ET_SEGM)
                {
                    const ScalarSpaceTimeFiniteElement<1> &  fel_st 
                        = dynamic_cast<const ScalarSpaceTimeFiniteElement<1> & >(fel);
                    ScalarFEEvaluator<1> lset_eval(fel_st, linvec, lh);
                    PointContainer<2> pc;
                    CompositeQuadratureRule<2> compositerule;
                    NumericalIntegrationStrategy<ET_SEGM,ET_SEGM> numint(lset_eval, pc,
                                                                         compositerule, 
                                                                         lh,
                                                                         order_space, order_time,
                                                                         num_int_ref_space, 
                                                                         num_int_ref_time);
                    numint.SetVerticesSpace();
                    numint.SetVerticesTime();
                    // numint.SetDistanceThreshold(2.0*absdet+vmax/dt);
                    {
                        static Timer timer ("npxtest - MakeQuadRule - total");
                        RegionTimer reg (timer);
                        sels_of_dt[numint.MakeQuadRule()]++;
                    }
#pragma omp critical(bndmeaspos)
                    for (int i = 0; i < compositerule.quadrule_pos.Size(); ++i)
                        meas_of_dt_bnd[POS] += absdet * compositerule.quadrule_pos.weights[i];
#pragma omp critical(bndmeasneg)
                    for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
                        meas_of_dt_bnd[NEG] += absdet * compositerule.quadrule_neg.weights[i];
//                 for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
//                 {
//                     static Timer timer ("npxtest - Transform n");
//                     RegionTimer reg (timer);
//                     Vec<3> st_point = compositerule.quadrule_if.points[i];
//                     IntegrationPoint ip(st_point(0),st_point(1));
//                     MappedIntegrationPoint<D,D> mip(ip,eltrans);
//                     Mat<2,2> Finv = mip.GetJacobianInverse();
//                     Vec<3> nref = compositerule.quadrule_if.normals[i];
//                     Vec<2> nref_sp (nref(0),nref(1));
//                     Vec<2> n_sp = absdet * Trans(Finv) * nref_sp ;
//                     double n_t = nref(2) * absdet;
//                     Vec<3> n (n_sp(0),n_sp(1),n_t);
//                     double fac = L2Norm(n);
// #pragma omp critical(measif)
//                     meas_of_dt[IF] +=  compositerule.quadrule_if.weights[i] * fac;
//                 }
//                 if (output)
//                 {
//                     for (int i = 0; i < compositerule.quadrule_neg.Size(); ++i)
//                     {
//                         static Timer timer ("npxtest - output neg");
//                         RegionTimer reg (timer);
//                         Vec<3> st_point = compositerule.quadrule_neg.points[i];
//                         IntegrationPoint ip(st_point(0),st_point(1));
//                         MappedIntegrationPoint<D,D> mip(ip,eltrans);
//                         Vec<2> mapped_point = mip.GetPoint();
//                         outneg_st << mapped_point(0) << "\t" << mapped_point(1) << "\t" << st_point(2) << endl;
//                     }
//                     outneg_st << endl;
//                     outneg_st << endl;
//                     for (int i = 0; i < compositerule.quadrule_if.Size(); ++i)
//                     {
//                         static Timer timer ("npxtest - output if");
//                         RegionTimer reg (timer);
//                         Vec<3> st_point = compositerule.quadrule_if.points[i];
//                         IntegrationPoint ip(st_point(0),st_point(1));
//                         MappedIntegrationPoint<D,D> mip(ip,eltrans);
//                         Vec<2> mapped_point = mip.GetPoint();
//                         outif_st << mapped_point(0) << "\t" << mapped_point(1) << "\t" << st_point(2) << endl;
//                     }
//                     outif_st << endl;
//                     outif_st << endl;
//                 }
                }

                if( et_space == ET_TRIG && et_time == ET_SEGM)
                {
                  cout << " no implementation for testing boundary of space-time 3D example " << endl;
                }

                if( et_space == ET_SEGM && et_time == ET_POINT)
                {
                  cout << " no implementation for testing boundary of (stationary) 2D example " << endl;
                }

                if( et_space == ET_TRIG && et_time == ET_POINT)
                {
                  cout << " no implementation for testing boundary of (stationary) 3D example " << endl;
                }

            }

            cout << " pos surf. elements : " << sels_of_dt[POS] << endl;
            cout << " neg surf. elements : " << sels_of_dt[NEG] << endl;
            cout << " cut surf. elements : " << sels_of_dt[IF] << endl;
            cout << " pos measure bndary : " << meas_of_dt_bnd[POS] << endl;
            cout << " neg measure bndary : " << meas_of_dt_bnd[NEG] << endl;
            cout << " cut measure bndary : " << meas_of_dt_bnd[IF] << endl;

        }
    }
};

static RegisterNumProc<NumProcTestXFEM<2> > npinittestxfem2d("testxfem");
static RegisterNumProc<NumProcTestXFEM<3> > npinittestxfem3d("testxfem3d");
