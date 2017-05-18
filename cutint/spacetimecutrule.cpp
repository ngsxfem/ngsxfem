#include "spacetimecutrule.hpp"
#include "../spacetime/myElement.hpp"

namespace xintegration
{
    vector<double> root_finding(SliceVector<> li, ScalarFiniteElement<1>* fe_time, LocalHeap& lh, int subdivs=50, int bisection_iterations = 70){
        // if(li.Size() == 2){
       if(fe_time->Order() == 0)
         return {};
       else if(fe_time->Order() == 1){
            // we assume to know the basis here:
            // phi0 =  1 - x
            // phi1 =      x
            if((li[0] >= 0) != (li[1] >= 0)){
                return {-li[0]/(li[1] - li[0]) };
            }
            else return {};
        }
        else if(fe_time->Order() == 2){
           // we assume to know the basis here:
           // phi0 = (1-x)*(1-2x) = 1 - 3*x + 2*x*x,
           // phi1 =    4*x*(1-x) =     4*x - 4*x*x,
           // phi2 =     x*(2x-1) =     - x + 2*x*x.
           //  sum_i c_i phi_i(x) = c + b*x + a*x*x,
           double c = li[0], a = 2*li[0]+2*li[2]-4.*li[1], b = li[2] - a - c;
           vector<double> roots;
           const double w = b*b - 4*a*c;
           if(abs(w) < 1e-12){
               double x_ast = -b/(2*a);
               if(x_ast < 1 && x_ast > 0) roots.push_back(x_ast);
           }
           else {
               // if w<0 x_ast_i is nan (not in (0,1))
               double x_ast1 = (b + sqrt(w))/(-2.0*a);
               double x_ast2 = (b - sqrt(w))/(-2.0*a);

               if(x_ast1 < 1 && x_ast1 > 0) roots.push_back(x_ast1);
               if(x_ast2 < 1 && x_ast2 > 0) roots.push_back(x_ast2);
           }
           return roots;
       }
        else {
            static bool first = true;
            if (first)
            {
              cout << "Calling bisection ..." << endl;
              first = false;
            }
            vector<double> vals(subdivs); vector<tuple<double,double>> sign_change_intervals;
            vector<double> roots; double delta_x = 1./subdivs;
            FlatVector<> shape(li.Size(), lh);
            function<double(double)> eval = [&li, &fe_time, &shape](double xi) -> double {
                fe_time->CalcShape(IntegrationPoint(Vec<3>{xi,0,0}, 0.), shape);
                return InnerProduct(li,shape);
            };
            for(int i=0; i<subdivs+1; i++){
                double xi = delta_x*i;
                vals[i] = eval(xi);
                if(vals[i] == 0)
                  roots.push_back(xi);
                if(i >= 1) if(vals[i-1] * vals[i] < 0) sign_change_intervals.push_back(make_tuple( xi-delta_x, xi));
            }
            for(auto interval : sign_change_intervals){
                double a = get<0>(interval), b = get<1>(interval); double x_mid;
                double aval = eval(a), bval = eval(b);
                for(int j=0; j<bisection_iterations; j++){
                  
                  //secant rule stopping criteria
                  const double x_lin = a - aval*(b-a)/(bval-aval);
                  double val = eval(x_lin);
                  if (2*abs(val) < 1e-12)
                  {
                    a=b=x_lin;
                    break;
                  }

                  x_mid = 0.5*(a+b);
                  val = eval(x_mid);
                    if(val == 0) break;
                    if(val * aval < 0) {
                        b = x_mid; bval = val;
                    }
                    else if(val * bval < 0){
                        a = x_mid; aval = val;
                    }
                    else throw Exception("Strange sign structure during bisection!");
                }
                roots.push_back(0.5*(a+b));
            }
            return roots;
        }
    }

    const IntegrationRule * SpaceTimeCutIntegrationRule(FlatVector<> cf_lset_at_element,
                                                        ELEMENT_TYPE et_space,
                                                        ScalarFiniteElement<1>* fe_time,
                                                        DOMAIN_TYPE dt,
                                                        int order_time,
                                                        int order_space,
                                                        LocalHeap & lh){
        int lset_nfreedofs = cf_lset_at_element.Size();
        int space_nfreedofs = ElementTopology::GetNVertices(et_space);
        int time_nfreedofs = lset_nfreedofs / space_nfreedofs;
        FlatMatrix<> lset_st(time_nfreedofs, space_nfreedofs, &cf_lset_at_element(0,0));

        vector<double> cut_points{0,1};
        for(int i=0; i<space_nfreedofs; i++){
            auto li = lset_st.Col(i);
            auto cp = root_finding(li, fe_time, lh);
            if(cp.size() > 0) cut_points.insert(cut_points.begin(), cp.begin(), cp.end());
        }
        sort(cut_points.begin(), cut_points.end());

        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, order_time);
        auto ir = new (lh) IntegrationRule();

        for(int i=0; i<cut_points.size() -1; i++){
            double t0 = cut_points[i], t1 = cut_points[i+1];
            for(auto ip:ir_time){
                double t = t0 + ip.Point()[0]*(t1 - t0);
                FlatVector<> cf_lset_at_t(space_nfreedofs,lh);
                FlatVector<> shape(time_nfreedofs, lh);
                fe_time->CalcShape(IntegrationPoint(Vec<3>{t,0,0}, 0.), shape);
                cf_lset_at_t = Trans(lset_st)*shape;

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
        //auto time_fe = new NodalTimeFE(1);
        NodalTimeFE time_fe(1);
        auto ir_neg = SpaceTimeCutIntegrationRule(a, ET_TRIG, &time_fe, NEG, 0, 0, lh);
        auto ir_pos = SpaceTimeCutIntegrationRule(a, ET_TRIG,&time_fe, POS, 0, 0, lh);
        cout << "IR neg: " << *ir_neg << endl << "IR pos: " << *ir_pos << endl;
        double V_pos, V_neg;
        for(auto ip: *ir_neg) V_neg += ip.Weight();
        for(auto ip: *ir_pos) V_pos += ip.Weight();

        cout << "V_pos: " << V_pos << endl << "V_neg: " << V_neg << endl;
        cout << "V_pos + V_neg = " << V_pos + V_neg << endl;

        //Testing of the root_finding function
        cout << "Testing of the root finding function: " << endl;

        NodalTimeFE time_fe2(2);
        auto roots = root_finding(Vector<>{-2, 0.1, 1}, &time_fe2, lh);
        cout << "The roots: " << endl;
        for(auto d:roots) cout << d << endl;

        cout << "Testing of the problematic lset function: " << endl;
        Vector<> iamproblematic{0.00673963, -0.0999511, -0.0127287, 0.0220921, -0.0849581, 0.000203323};
        auto ir_neg2 = SpaceTimeCutIntegrationRule(iamproblematic, ET_TRIG, &time_fe, POS, 0,0,lh);
    }
}
