#include "spacetimecutrule.hpp"
#include "../spacetime/SpaceTimeFE.hpp"
#include "../utils/ngsxstd.hpp"

namespace xintegration
{

    vector<double> root_finding(SliceVector<> li, ScalarFiniteElement<1>* fe_time, LocalHeap& lh, int subdivs=50, int bisection_iterations = 70){
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
           if(abs(a) < 1e-12) {
               double x_ast = -c/b;
               if(x_ast < 1 && x_ast > 0) roots.push_back(x_ast);
           }
           else {
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
           }
           return roots;
       }
        else {
            static bool first = true;
            if (first)
            {
              cout << IM(3) << "Calling bisection for root finding ..." << endl;
              first = false;
            }
            vector<double> vals(subdivs+1); vector<tuple<double,double>> sign_change_intervals;
            vector<double> roots; double delta_x = 1./subdivs;
            FlatVector<> shape(li.Size(), lh);
            function<double(double)> eval = [&li, &fe_time, &shape](double xi) -> double {
                fe_time->CalcShape(IntegrationPoint(Vec<3>{xi,0,0}, 0.), shape);
                return InnerProduct(li,shape);
            };
            for(int i=0; i<subdivs+1; i++){
                double xi = delta_x*i;
                vals[i] = eval(xi);
                if(vals[i] == 0) vals[i] = globxvar.EPS_STCR_ROOT_SEARCH_BISECTION; //roots.push_back(xi);
                if(i >= 1) if(vals[i-1]*vals[i]<0) sign_change_intervals.push_back(make_tuple( xi-delta_x, xi));
            }

            for(auto interval : sign_change_intervals){
                double a = get<0>(interval), b = get<1>(interval); double x_mid;
                double aval = eval(a), bval = eval(b);
                int j;
                for(j=0; j<bisection_iterations; j++){
                  
                  //secant rule stopping criteria
                  const double x_lin = a - aval*(b-a)/(bval-aval);
                  double val = eval(x_lin);
                  if (2*abs(val) < globxvar.EPS_STCR_ROOT_SEARCH_BISECTION)
                  {
                    a=b=x_lin;
                    break;
                  }

                  x_mid = 0.5*(a+b);
                  val = eval(x_mid);
                    if(val == 0) break;
                    if(val * aval < 0){
                        b = x_mid; bval = val;
                    }
                    else if(val * bval < 0){
                        a = x_mid; aval = val;
                    }
                    else throw Exception("Strange sign structure during bisection!");
                }
                if(j == bisection_iterations)
                    cout << IM(2) << "WARNING: Bisection search did not converge. Residual: " << eval(0.5*(a+b)) << endl;
                roots.push_back(0.5*(a+b));
            }
            return roots;
        }
    }

    tuple<const IntegrationRule *, Array<double>> SpaceTimeCutIntegrationRule(FlatVector<> cf_lset_at_element,
                                                        const ElementTransformation &trafo,
                                                        ScalarFiniteElement<1>* fe_time,
                                                        DOMAIN_TYPE dt,
                                                        int order_time,
                                                        int order_space,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                        LocalHeap & lh){
        static Timer timer("SpaceTimeCutIntegrationRule");
        RegionTimer rt(timer);
        ELEMENT_TYPE et_space = trafo.GetElementType();
        
        int lset_nfreedofs = cf_lset_at_element.Size();
        int space_nfreedofs = ElementTopology::GetNVertices(et_space);
        int time_nfreedofs = lset_nfreedofs / space_nfreedofs;
        FlatMatrix<> lset_st(time_nfreedofs, space_nfreedofs, &cf_lset_at_element(0,0));

        vector<double> cut_points{0,1};
        if (globxvar.DO_NAIVE_TIMEINT){
            bool haspos = false;
            bool hasneg = false;
            for(auto d : cf_lset_at_element){
                if (d < 0) hasneg = true;
                if (d > 0) haspos = true;
            }

            if(globxvar.NAIVE_TIMEINT_SUBDIVS < 1) throw Exception("NAIVE_TIMEINT_SUBDIVS < 1 is not possible");
            else {
                if(hasneg && haspos){
                    for(int i=1; i<globxvar.NAIVE_TIMEINT_SUBDIVS; i++) cut_points.push_back(((double)i)/(globxvar.NAIVE_TIMEINT_SUBDIVS));
                }
            }
        }
        else {
            for(int i=0; i<space_nfreedofs; i++){
                auto li = lset_st.Col(i);
                auto cp = root_finding(li, fe_time, lh);

                if(cp.size() > 0) cut_points.insert(cut_points.begin(), cp.begin(), cp.end());
            }
        }
        sort(cut_points.begin(), cut_points.end());

        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, globxvar.DO_NAIVE_TIMEINT ? globxvar.NAIVE_TIMEINT_ORDER : order_time);
        if(order_space == -1) order_space = 5;
        const IntegrationRule & stdir = SelectIntegrationRule (et_space, order_space);
        const int MAXSIZE_PER = 5 * stdir.Size();
        const int MAXSIZE = MAXSIZE_PER * (cut_points.size()-1) * ir_time.Size();
        // MAXSIZE is an estimated maximum size for the IntegrationRule 
        // (for fixed memory allocation on the LocalHeap)
        IntegrationRule * ir = new (lh) IntegrationRule(MAXSIZE,lh); ir->SetSize0(); //local size 0, physical size MAXSIZE
        Array<double> wei_arr(MAXSIZE,lh); // will be resized at the end
        int ip_counter = 0;

        for(int i=0; i<cut_points.size() -1; i++){
            double t0 = cut_points[i], t1 = cut_points[i+1];
            for(auto ip:ir_time){
                double t = t0 + ip.Point()[0]*(t1 - t0);
                FlatVector<> cf_lset_at_t(space_nfreedofs,lh);
                FlatVector<> shape(time_nfreedofs, lh);

                fe_time->CalcShape(IntegrationPoint(Vec<3>{t,0,0}, 0.), shape);
                cf_lset_at_t = Trans(lset_st)*shape;
                for(auto &d : cf_lset_at_t) if(abs(d) < globxvar.EPS_STCR_LSET_PERTUBATION){
                    static bool first = true;
                    if (first) {
                        cout << IM(4) << "The ST eps pertubation trick has been applied" << endl;
                        first = false;
                    }
                    if(d >= 0) d = globxvar.EPS_STCR_LSET_PERTUBATION;
                    else d = -globxvar.EPS_STCR_LSET_PERTUBATION;
                }

                auto element_domain = CheckIfStraightCut(cf_lset_at_t);

                const int offset = ir->Size();
                if (element_domain == IF)
                {
                    auto spir = StraightCutIntegrationRule(cf_lset_at_t, trafo, dt, order_space, quad_dir_policy, lh, true, t);
                    ir->Append(*spir);
                    ip_counter += spir->Size();
                }
                else if (element_domain == dt)
                {
                    ir->Append(stdir);
                    ip_counter += stdir.Size();
                }

                for(int k = offset; k < ir->Size(); k++) {
                    wei_arr[k] = (*ir)[k].Weight()*ip.Weight()*(t1-t0);
                    (*ir)[k].SetWeight(t);
                    MarkAsSpaceTimeIntegrationPoint((*ir)[k]);
                }
                /*                                     
                CutSimplexElementGeometry geom(cf_lset_at_t, et_space, lh);
                CutQuadElementGeometry geom_quad(cf_lset_at_t, et_space, lh);

                if (element_domain == IF)
                {
                  if((et_space == ET_QUAD)||(et_space == ET_HEX)) geom_quad.GetIntegrationRule(order_space, dt, quad_untrafo);
                  else geom.GetIntegrationRule(order_space, dt, quad_untrafo);
                }
                else if(dt == element_domain)
                {
                  quad_untrafo.Append(SelectIntegrationRule (et_space, order_space));
                }
                for(IntegrationPoint& ip2 : quad_untrafo) {
                    if(trafo.SpaceDim() == 1) ip2.Point()[1] = t;
                    if(trafo.SpaceDim() == 2) ip2.Point()[2] = t;
                    ip2.SetWeight(ip2.Weight()*ip.Weight()*(t1-t0));
                }
                if((element_domain == IF)&&(dt == IF)){
                    auto ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
                    if(trafo.SpaceDim() == 1) TransformQuadUntrafoToIRInterface<1>(quad_untrafo, trafo, geom, ir_interface);
                    else if (trafo.SpaceDim() == 2){
                        if(et_space == ET_QUAD){
                            TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom_quad, ir_interface);
                        }
                        else TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom, ir_interface);
                    }
                    else{
                        if(et_space == ET_HEX){
                            TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom_quad, ir_interface);
                        }
                        else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom, ir_interface);
                    }
                    ir->Append(*ir_interface);
                }
                else ir->Append(quad_untrafo);
                */
            }
        }
        wei_arr.SetSize(ir->Size());
        if (ip_counter > MAXSIZE)
            throw Exception("memory allocation for integration rule was insufficient");

        if (ir->Size() == 0)
            return make_tuple(nullptr, wei_arr);
        else
            return make_tuple(ir, wei_arr);
    }

    tuple<const IntegrationRule *, Array<double>> SpaceTimeCutIntegrationRuleUntransformed(FlatVector<> cf_lset_at_element,
                                                        ELEMENT_TYPE et_space,
                                                        ScalarFiniteElement<1>* fe_time,
                                                        DOMAIN_TYPE dt,
                                                        int order_time,
                                                        int order_space,
                                                        SWAP_DIMENSIONS_POLICY quad_dir_policy,
                                                        LocalHeap & lh){
        static Timer timer("SpaceTimeCutIntegrationRule");
        RegionTimer rt(timer);

        int lset_nfreedofs = cf_lset_at_element.Size();
        int space_nfreedofs = ElementTopology::GetNVertices(et_space);
        int time_nfreedofs = lset_nfreedofs / space_nfreedofs;
        FlatMatrix<> lset_st(time_nfreedofs, space_nfreedofs, &cf_lset_at_element(0,0));

        vector<double> cut_points{0,1};
        if (globxvar.DO_NAIVE_TIMEINT){
            bool haspos = false;
            bool hasneg = false;
            for(auto d : cf_lset_at_element){
                if (d < 0) hasneg = true;
                if (d > 0) haspos = true;
            }

            if(globxvar.NAIVE_TIMEINT_SUBDIVS < 1) throw Exception("NAIVE_TIMEINT_SUBDIVS < 1 is not possible");
            else {
                if(hasneg && haspos){
                    for(int i=1; i<globxvar.NAIVE_TIMEINT_SUBDIVS; i++) cut_points.push_back(((double)i)/(globxvar.NAIVE_TIMEINT_SUBDIVS));
                }
            }
        }
        else {
            for(int i=0; i<space_nfreedofs; i++){
                auto li = lset_st.Col(i);
                auto cp = root_finding(li, fe_time, lh);

                if(cp.size() > 0) cut_points.insert(cut_points.begin(), cp.begin(), cp.end());
            }
        }
        sort(cut_points.begin(), cut_points.end());

        const IntegrationRule & ir_time = SelectIntegrationRule(ET_SEGM, globxvar.DO_NAIVE_TIMEINT ? globxvar.NAIVE_TIMEINT_ORDER : order_time);
        if(order_space == -1) order_space = 5;
        const IntegrationRule & stdir = SelectIntegrationRule (et_space, order_space);
        const int MAXSIZE_PER = 5 * stdir.Size();
        const int MAXSIZE = MAXSIZE_PER * (cut_points.size()-1) * ir_time.Size();
        // MAXSIZE is an estimated maximum size for the IntegrationRule
        // (for fixed memory allocation on the LocalHeap)
        IntegrationRule * ir = new (lh) IntegrationRule(MAXSIZE,lh); ir->SetSize0(); //local size 0, physical size MAXSIZE
        Array<double> wei_arr(MAXSIZE,lh); // will be resized at the end
        int ip_counter = 0;

        for(int i=0; i<cut_points.size() -1; i++){
            double t0 = cut_points[i], t1 = cut_points[i+1];
            for(auto ip:ir_time){
                double t = t0 + ip.Point()[0]*(t1 - t0);
                FlatVector<> cf_lset_at_t(space_nfreedofs,lh);
                FlatVector<> shape(time_nfreedofs, lh);

                fe_time->CalcShape(IntegrationPoint(Vec<3>{t,0,0}, 0.), shape);
                cf_lset_at_t = Trans(lset_st)*shape;
                for(auto &d : cf_lset_at_t) if(abs(d) < globxvar.EPS_STCR_LSET_PERTUBATION){
                    static bool first = true;
                    if (first) {
                        cout << IM(4) << "The ST eps pertubation trick has been applied" << endl;
                        first = false;
                    }
                    if(d >= 0) d = globxvar.EPS_STCR_LSET_PERTUBATION;
                    else d = -globxvar.EPS_STCR_LSET_PERTUBATION;
                }

                auto element_domain = CheckIfStraightCut(cf_lset_at_t);

                const int offset = ir->Size();
                if (element_domain == IF)
                {
                    auto spir = StraightCutIntegrationRuleUntransformed(cf_lset_at_t, et_space, dt, order_space, quad_dir_policy, lh);
                    ir->Append(*spir);
                    ip_counter += spir->Size();
                }
                else if (element_domain == dt)
                {
                    ir->Append(stdir);
                    ip_counter += stdir.Size();
                }

                for(int k = offset; k < ir->Size(); k++) {
                    wei_arr[k] = (*ir)[k].Weight()*ip.Weight()*(t1-t0);
                    (*ir)[k].SetWeight(t);
                    MarkAsSpaceTimeIntegrationPoint((*ir)[k]);
                }
                /*
                CutSimplexElementGeometry geom(cf_lset_at_t, et_space, lh);
                CutQuadElementGeometry geom_quad(cf_lset_at_t, et_space, lh);

                if (element_domain == IF)
                {
                  if((et_space == ET_QUAD)||(et_space == ET_HEX)) geom_quad.GetIntegrationRule(order_space, dt, quad_untrafo);
                  else geom.GetIntegrationRule(order_space, dt, quad_untrafo);
                }
                else if(dt == element_domain)
                {
                  quad_untrafo.Append(SelectIntegrationRule (et_space, order_space));
                }
                for(IntegrationPoint& ip2 : quad_untrafo) {
                    if(trafo.SpaceDim() == 1) ip2.Point()[1] = t;
                    if(trafo.SpaceDim() == 2) ip2.Point()[2] = t;
                    ip2.SetWeight(ip2.Weight()*ip.Weight()*(t1-t0));
                }
                if((element_domain == IF)&&(dt == IF)){
                    auto ir_interface  = new (lh) IntegrationRule(quad_untrafo.Size(),lh);
                    if(trafo.SpaceDim() == 1) TransformQuadUntrafoToIRInterface<1>(quad_untrafo, trafo, geom, ir_interface);
                    else if (trafo.SpaceDim() == 2){
                        if(et_space == ET_QUAD){
                            TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom_quad, ir_interface);
                        }
                        else TransformQuadUntrafoToIRInterface<2>(quad_untrafo, trafo, geom, ir_interface);
                    }
                    else{
                        if(et_space == ET_HEX){
                            TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom_quad, ir_interface);
                        }
                        else TransformQuadUntrafoToIRInterface<3>(quad_untrafo, trafo, geom, ir_interface);
                    }
                    ir->Append(*ir_interface);
                }
                else ir->Append(quad_untrafo);
                */
            }
        }
        wei_arr.SetSize(ir->Size());
        if (ip_counter > MAXSIZE)
            throw Exception("memory allocation for integration rule was insufficient");

        if (ir->Size() == 0)
            return make_tuple(nullptr, wei_arr);
        else
            return make_tuple(ir, wei_arr);
    }

}
