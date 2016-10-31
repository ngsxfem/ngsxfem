#include <iostream>
#include <vector>
#include <fstream>
#include <vector>
#include <memory>
#include <tuple>
#include <set>
#include <cmath>

using namespace std;

ofstream visu;
string new_picture = "\\end{tikzpicture} \n \\\\ \\hfill \\\\ \\begin{tikzpicture} [scale=5] \n";

bool operator==(const vector<double> &a, const vector<double> &b){
    double sum=0;
    for(unsigned int i=0; i<a.size(); i++){
        sum += pow(a[i]-b[i],2);
        if(sum > 1e-10) return false;
    }
    return true;
}

vector<double> operator*(const vector<double> &a, double b){
    vector<double> c(a.size());
    for(unsigned int i=0; i<a.size(); i++) c[i] = a[i]*b;
    return c;
}

vector<double> operator*(double a, const vector<double> &b){
    return b*a;
}

vector<double> operator+(const vector<double>& a, const vector<double>& b){
    vector<double> c(a.size());
    for(unsigned int i=0; i<a.size(); i++) c[i] = a[i] + b[i];
    return c;
}

template<unsigned int D>
class PolytopE{
public:
    vector<PolytopE<D-1>> childs;
    PolytopE(initializer_list<PolytopE<D-1>> a_childs) : childs(a_childs) {;}
    PolytopE() {;}

    void Partition(bool sign);

    template<typename U>
    void Execute_on_every_point(function<void(PolytopE<0>&, U&)> f, U& arg) {
        for(auto &t: childs) t.Execute_on_every_point(f, arg);
    }

    double GetLsetSum();

    void FindAllElements(set<int> &s){
        for(auto& t:childs) t.FindAllElements(s);
    }

    void FindNeighbours(int i_searched, set<int> &s){
        for(auto& t:childs) t.FindNeighbours(i_searched, s);
    }

    bool CanStay(bool sign){
        for(auto t:childs) if(!t.CanStay(sign)) return false;
        return true;
    }

    vector<PolytopE<D>> PartitionViaArtificalLevelset(double dominant_sign=0);
    vector<PolytopE<D>> DecomposeViaArtificalLevelset(const PolytopE<D-1> &element_repl);

    bool IsSimplex();

    vector<PolytopE<D>> Partition();
    vector<PolytopE<D>> Decompose(double dominant_sign=0);

    shared_ptr<vector<vector<double>>> GetPointCointainerPtr(){
        return childs[0].GetPointCointainerPtr();
    }
    shared_ptr<vector<double>> GetLsetSetPointer(){
        return childs[0].GetLsetSetPointer();
    }
};

template<unsigned int D>
ostream& operator <<(ostream& stream, const PolytopE<D>& p){
    for(auto t: p.childs) stream << t << " ";
    stream << endl;
    return stream;
}

template<>
class PolytopE<0> {
public:
    int i;
    shared_ptr<vector<vector<double>>> pointset;
    shared_ptr<vector<double>> lsetset;

    template<typename U>
    void Execute_on_every_point(function<void(PolytopE<0>&, U&)> f, U &arg){
        f(*this, arg);
    }

    vector<double> p() const{
        return (*pointset)[i];
    }
    double lset() const{
        return (*lsetset)[i];
    }
    void Setlset(double a_lset) {
        (*lsetset)[i] = a_lset;
    }

    bool CanStay(bool sign){
        if(sign){
            if (lset() > -1e-10) return true;
            else return false;
        }
        else{
            if(lset() < 1e-10) return true;
            else return false;
        }
    }
    void FindAllElements(set<int> &s){
        s.insert(i);
    }
    PolytopE<0>(int a_i, shared_ptr<vector<vector<double>>> a_pointset, shared_ptr<vector<double>> a_lsetset): i(a_i), pointset(a_pointset), lsetset(a_lsetset) {;}

    shared_ptr<vector<vector<double>>> GetPointCointainerPtr(){
        return pointset;
    }
    shared_ptr<vector<double>> GetLsetSetPointer(){
        return lsetset;
    }
};

bool operator==(const PolytopE<0>& a, const PolytopE<0>& b){
    return a.i == b.i;
}

template<>
ostream& operator <<(ostream& stream, const PolytopE<1>& p){
    stream << "\\draw[->] (" << p.childs[0].p()[0] << ", " << p.childs[0].p()[1] << ", " << p.childs[0].p()[2] << ") -- ";
    stream << " (" << p.childs[1].p()[0] << ", " << p.childs[1].p()[1] << ", " << p.childs[1].p()[2] << "); ";
    stream << endl;
    return stream;
}

template<>
void PolytopE<1>::FindNeighbours(int i_searched, set<int> &s) {
    if(i_searched == childs[0].i) s.insert(childs[1].i);
    if(i_searched == childs[1].i) s.insert(childs[0].i);
}

template<unsigned int D>
bool IsCut(const PolytopE<D> &p){
    for(auto t: p.childs) if(IsCut(t)) return true;
    return false;
}

template<>
bool IsCut(const PolytopE<1> &p){
    return p.childs[0].lset()*p.childs[1].lset() < -1e-10;
}

template<unsigned int D>
double PolytopE<D>::GetLsetSum(){
    double lsetsum = 0;
    function<void(PolytopE<0>&, double&)> f = [] (PolytopE<0>& p, double& d) {d += p.lset();};
    Execute_on_every_point(f, lsetsum);
    return lsetsum;
}

template<>
vector<PolytopE<1>> PolytopE<1>::Partition(){
    vector<PolytopE<1>> new_lines;

    if(IsCut(*this)){
        cout << "Cut found!" << endl;
        double a = childs[0].lset(); double b = childs[1].lset();
        auto p1 = childs[0].p(); auto p2 = childs[1].p();
        vector<double> p = p1+(a/(a-b)*(p2+(-1.)*p1));

        auto pointsetptr = childs[0].pointset;
        auto lsetptr = childs[0].lsetset;

        int idx = -1;
        for(unsigned int i=0; i<pointsetptr->size(); i++) {
            if(p == (*pointsetptr)[i]) idx = i;
        }
        if(idx == -1) { (*pointsetptr).push_back(p); (*lsetptr).push_back(0); idx = pointsetptr->size()-1; cout << "Inserted new point!" << endl; }

        PolytopE<0> new_point = PolytopE<0>(idx, pointsetptr, lsetptr);
        //The order in the vector for the parts is -,+
        if(a<b){
            new_lines.push_back(PolytopE<1>{childs[0], new_point});
            new_lines.push_back(PolytopE<1>{childs[1], new_point});
        }
        else {
            new_lines.push_back(PolytopE<1>{childs[1], new_point});
            new_lines.push_back(PolytopE<1>{childs[0], new_point});
        }
    }
    else {
        cout << "No cut found!" << endl;
        new_lines.push_back(*this);
    }
    cout << "Finished with cutting this line" << endl;
    return new_lines;
}

template<unsigned int D>
vector<PolytopE<D>> PolytopE<D>::Partition(){
    vector<PolytopE<D>> new_figures;
    if(IsCut(*this)){
        new_figures.resize(2);
        PolytopE<D-1> new_Codim1Fig;
        for(auto t: childs) {
            vector<PolytopE<D-1>> tp = t.Partition();
            if(tp.size() > 1){
                for(int i=0; i<2; i++) new_figures[i].childs.push_back(tp[i]);
                new_Codim1Fig.childs.push_back(tp[0].childs[tp[0].childs.size()-1]);
            }
            else {
                double s = tp[0].GetLsetSum();
                if(s<-1e-8) new_figures[0].childs.push_back(tp[0]);
                else if (s>1e-8) new_figures[1].childs.push_back(tp[0]);
                else {
                    new_figures[0].childs.push_back(tp[0]);
                    new_figures[1].childs.push_back(tp[0]);
                }
            }
        }
        new_figures[0].childs.push_back(new_Codim1Fig);
        new_figures[1].childs.push_back(new_Codim1Fig);
    }
    else {
        new_figures.push_back(*this);
    }
    return new_figures;
}

template<unsigned int D>
bool PolytopE<D>::IsSimplex(){
    set<int> elem; FindAllElements(elem);
    return elem.size() == D+1;
}

template<>
vector<PolytopE<2>> PolytopE<2>::DecomposeViaArtificalLevelset(const PolytopE<1>& element_repl){
    vector<PolytopE<2>> new_figures(2);
    for(auto t:childs){
        double s = t.GetLsetSum();
        if(s<-1e-10) new_figures[0].childs.push_back(t);
        else if(s>1e-10) new_figures[1].childs.push_back(t);
        else {
            new_figures[0].childs.push_back(t);
            new_figures[1].childs.push_back(t);
        }
    }
    new_figures[0].childs.push_back(element_repl);
    new_figures[1].childs.push_back(element_repl);
    return new_figures;
}

template<>
vector<PolytopE<3>> PolytopE<3>::DecomposeViaArtificalLevelset(const PolytopE<2>& element_repl){
    auto pointsetptr = GetPointCointainerPtr();
    auto lsetptr = GetLsetSetPointer();
    vector<PolytopE<3>> new_figures(2);
    for(auto t:childs){
        //visu << t << new_picture;
        if((!t.CanStay(true)) && !(t.CanStay(false))){
            //Maybe this could be done with better Performance using element_repl...
            set<int> s1; t.FindAllElements(s1);
            PolytopE<1> p;
            for(int i : s1) if(abs((*lsetptr)[i]) < 1e-10) p.childs.push_back(PolytopE<0>(i, pointsetptr, lsetptr));
            auto dec = t.DecomposeViaArtificalLevelset(p);
            new_figures[0].childs.push_back(dec[0]);
            new_figures[1].childs.push_back(dec[1]);
        }
        else{
            double s = t.GetLsetSum();
            if(s<-1e-8) { new_figures[0].childs.push_back(t); cout << "This plane belongs to neg side." << endl; }
            else if(s>1e-8){ new_figures[1].childs.push_back(t); cout << "This plane belongs to pos side." << endl; }
            else {
                new_figures[0].childs.push_back(t);
                new_figures[1].childs.push_back(t);
                cout << "This plane belongs to both sides." << endl;
            }
        }
    }
    new_figures[0].childs.push_back(element_repl);
    new_figures[1].childs.push_back(element_repl);
    return new_figures;
}

template<>
vector<PolytopE<2> > PolytopE<2>::PartitionViaArtificalLevelset(double dominant_sign){
    if(dominant_sign == 0) dominant_sign = GetLsetSum();
    auto pointsetptr = GetPointCointainerPtr();
    auto lsetptr = GetLsetSetPointer();

    PolytopE<1> new_element = childs[childs.size()-1];
    set<int> neighbours; int i_isolate = new_element.childs[0].i;
    (*lsetptr)[i_isolate] = dominant_sign > 0? -1 : 1;
    FindNeighbours(i_isolate, neighbours);

    PolytopE<1> new_element_repl;
    for(auto n : neighbours) {
        (*lsetptr)[n] = 0;
        new_element_repl.childs.push_back(PolytopE<0>(n, pointsetptr, lsetptr));
    }
    return DecomposeViaArtificalLevelset(new_element_repl);
}

PolytopE<2> make_Trig(PolytopE<0>& a, PolytopE<0>& b, PolytopE<0>& c){
    return {{a,b},{b,c},{a,c}};
}

template<>
vector<PolytopE<3> > PolytopE<3>::PartitionViaArtificalLevelset(double dominant_sign){
    if(dominant_sign == 0) dominant_sign = GetLsetSum();
    auto pointsetptr = GetPointCointainerPtr();
    auto lsetptr = GetLsetSetPointer();

    PolytopE<2> new_element = childs[childs.size()-1];
    int i_isolate; set<int> neighbours;
    set<int> new_element_i; new_element.FindAllElements(new_element_i);
    int least_neighbours = 100;
    for(auto i: new_element_i) {
        FindNeighbours(i, neighbours);
        cout << "Found " << neighbours.size() << " nbs of " << i << endl;
        if(neighbours.size() < least_neighbours) {i_isolate = i; least_neighbours = neighbours.size();}
        neighbours.clear();
        (*lsetptr)[i] = dominant_sign > 0? 1 : -1;
    }
    (*lsetptr)[i_isolate] = dominant_sign > 0? -1 : 1;
    cout << "Taking " << i_isolate << " as i_isolate." << endl;
    FindNeighbours(i_isolate, neighbours);

    vector<PolytopE<0>> new_element_repl;
    for(auto n : neighbours) {
        (*lsetptr)[n] = 0;
        new_element_repl.push_back(PolytopE<0>(n, pointsetptr, lsetptr));
    }
    auto element_repl = make_Trig(new_element_repl[0], new_element_repl[1], new_element_repl[2]);

    return DecomposeViaArtificalLevelset(element_repl);
}

template<unsigned int D>
vector<PolytopE<D>> PolytopE<D>::Decompose(double dominant_sign){
    if(dominant_sign == 0) dominant_sign = GetLsetSum();
    PolytopE<D> iteration_copy(*this);
    vector<PolytopE<D>> container;
    int i;
    for(i=0; i<20; i++){
        container.push_back(iteration_copy);
        if(iteration_copy.IsSimplex()) break;
        cout << "Part Neg is Not jet a Simplex" << endl;
        auto w = iteration_copy.PartitionViaArtificalLevelset(dominant_sign);
        visu << w[0] << new_picture << w[1] << new_picture;

        if(dominant_sign>1e-10) iteration_copy = w[1];
        else iteration_copy = w[0];
    }
    if(i == 19) cout << "Decomposition did not converge!" << endl;
    return container;
}

int main(int argc, char *argv[])
{
    vector<double> a{0,0,0}, b{1,0,0}, c{0,1,0}, d{0,0,1}, e{1,0,1};

    auto pointsetpointer = make_shared<vector<vector<double>>>(vector<vector<double>>{a,b,c,d,e});
    //auto lsetsetpointer = make_shared<vector<double>>(vector<double>{-1,-1,1,1,1}); //for trig
    auto lsetsetpointer = make_shared<vector<double>>(vector<double>{-1,-4,1,1,1});

    PolytopE<0> ap(0,pointsetpointer, lsetsetpointer);
    PolytopE<0> bp(1,pointsetpointer, lsetsetpointer);
    PolytopE<0> cp(2,pointsetpointer, lsetsetpointer);
    PolytopE<0> dp(3,pointsetpointer, lsetsetpointer);
    PolytopE<0> ep(4,pointsetpointer, lsetsetpointer);

    //PolytopE<2> tet = make_Trig(ap,bp,cp);
    PolytopE<3> tet{make_Trig(ap,bp,cp), make_Trig(ap,bp,dp), make_Trig(ap,cp,dp), make_Trig(bp,cp,dp)};
    cout << tet << endl;
    visu.open("visu.tex"); visu << " \\begin{tikzpicture}[scale=5]" << endl;
    visu << tet << new_picture;

    auto tet_part_pos = tet; auto tet_part_neg = tet;

    cout << "Will do the Partition now..." <<endl;
    auto v = tet.Partition();
    tet_part_neg = v[0]; tet_part_pos = v[1];

    visu << tet_part_neg << new_picture << tet_part_pos << new_picture;

    auto deco1 = tet_part_neg.Decompose();
    auto deco2 = tet_part_pos.Decompose();

    cout << "Finished with Partition. "<< endl;

    visu << " \\end{tikzpicture}" << endl;
    return 0;
}
