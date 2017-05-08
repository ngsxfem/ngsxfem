#include "spacetimecutrule.hpp"

namespace xintegration
{
    const IntegrationRule * SpaceTimeCutIntegrationRule(FlatVector<> cf_lset_at_element,
                                                        ELEMENT_TYPE et_space,
                                                        DOMAIN_TYPE dt,
                                                        int order_time,
                                                        int order_space,
                                                        LocalHeap & lh){
        int lset_nfreedofs = cf_lset_at_element.Size();
        int space_nfreedofs = ElementTopology::GetNVertices(et_space);
        if(lset_nfreedofs != 2*space_nfreedofs)
            throw Exception("Only P1 functions in time are allowed right now");
        if(et_space != ET_TRIG)
            throw Exception("Only Trigs as dt_space allowed right now");

        Array<Vec<2>> lset_on_space_vert(space_nfreedofs);
        for(int i=0; i<lset_nfreedofs; i++){
            lset_on_space_vert[i % space_nfreedofs][i/space_nfreedofs] = cf_lset_at_element[i];
        }
        cout << "The lset function values divided: " << lset_on_space_vert << endl;

        vector<double> cut_points{0,1};
        for(auto li : lset_on_space_vert){
            if(li[0]*li[1] < -1e-10){
                cut_points.push_back(-li[0]/(li[1] - li[0]));
            }
        }
        sort(cut_points.begin(), cut_points.end());
        cout << "The sorted cut points: " << endl;
        for(auto d: cut_points) cout << d << endl;

        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, order_time);
        auto ir = new (lh) IntegrationRule();

        for(int i=0; i<cut_points.size() -1; i++){
            double t0 = cut_points[i], t1 = cut_points[i+1];
            cout << "(t0,t1) = ( " << t0 << " , " << t1 << " )" << endl;
            for(auto ip:ir_time){
                double t = t0 + ip.Point()[0]*(t1 - t0);
                FlatVector<> cf_lset_at_t(3,lh);
                for(int j=0; j<3; j++) cf_lset_at_t[j] = lset_on_space_vert[j][0] + t*(lset_on_space_vert[j][1] - lset_on_space_vert[j][0]);
                cout << "t: " << t << endl;
                cout << "cf_lset_at_t: " << cf_lset_at_t << endl;

                DOMAIN_TYPE dt_at_t = CheckIfStraightCut(cf_lset_at_t);
                IntegrationRule quad_at_t;
                if(dt_at_t == IF){
                    CutSimplexElementGeometry geom(cf_lset_at_t, et_space, lh);
                    geom.GetIntegrationRule(order_space, dt, quad_at_t);
                }
                else if(dt_at_t == dt){
                    quad_at_t.Append(SelectIntegrationRule(ET_TRIG, order_space));
                }
                for(auto ip2 : quad_at_t) {
                    auto ip2Point = ip2.Point(); ip2Point[2] = t;
                    ir->Append(IntegrationPoint(ip2Point, ip2.Weight()*ip.Weight()*(t1-t0)));
                }
            }
        }
        return ir;
    }

    void DebugSpaceTimeCutIntegrationRule(){
        LocalHeap lh(10000);

        //vector<double> a2{-1.,-1.,1.,-1.,-1.,1.}; //Mimic the function phi(x,y,z) = 1- 2*x - 2*y
        vector<double> a2{-1.,-1., 1.,-3.,-3.,-1}; //Mimic the function phi(x,y,z) = 1- 2*x - 2*y - 2*z
        //vector<double> a2{1,1,1,-1,-2,-3}; //Maximal number of cuts in time
        FlatVector<> a(6,lh);
        for(int i=0; i<6; i++) a[i] = a2[i]; //Why isn't the Vec<6>{1,1,...} constructor working any more??
        auto ir_neg = SpaceTimeCutIntegrationRule(a, ET_TRIG, NEG, 2, 2, lh);
        auto ir_pos = SpaceTimeCutIntegrationRule(a, ET_TRIG, POS, 2, 2, lh);
        cout << "IR neg: " << *ir_neg << endl << "IR pos: " << *ir_pos << endl;
        double V_pos, V_neg;
        for(auto ip: *ir_neg) V_neg += ip.Weight();
        for(auto ip: *ir_pos) V_pos += ip.Weight();

        cout << "V_pos: " << V_pos << endl << "V_neg: " << V_neg << endl;
        cout << "V_pos + V_neg = " << V_pos + V_neg << endl;
    }
}
