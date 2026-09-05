// Exhaustive constants for ranks 2--6, with parity pruning and resumable tasks.
// Generalizes research/rank4/enumerate_constants.cpp; all arithmetic is exact.
// Usage: enumerate_constants RANK OUTPUT_DIR SHARD SHARDS [ONLY_UPPER_INDEX]
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using namespace std;
namespace fs = std::filesystem;
constexpr int MAX_N=6;
using Int=__int128_t;
using Vec=array<Int,MAX_N>;
using Mat=array<Vec,MAX_N>;
int n, bound;
Mat p{},m{};
set<vector<int>> answers,positive_cache;
unsigned long long nodes[MAX_N+1]{},labelled=0,lp_checks=0;
int current_upper=-1,current_partner=-1;
int ordered_limits[MAX_N][MAX_N]{};
double last_heartbeat=0;
auto started=chrono::steady_clock::now();
double seconds(){return chrono::duration<double>(chrono::steady_clock::now()-started).count();}
Int absolute(Int x){return x<0?-x:x;}
Int gcd128(Int a,Int b){while(b){Int r=a%b;a=b;b=r;}return a;}

Int det(Mat a,int size){
    if(!size)return 1;
    if(size==1)return a[0][0];
    if(size==2)return a[0][0]*a[1][1]-a[0][1]*a[1][0];
    // Fraction-free Bareiss elimination. All intermediate divisions are exact.
    Int previous=1,sign=1;
    for(int k=0;k<size-1;k++){
        int pivot=k;while(pivot<size&&!a[pivot][k])pivot++;
        if(pivot==size)return 0;
        if(pivot!=k){swap(a[pivot],a[k]);sign=-sign;}
        Int d=a[k][k];
        for(int i=k+1;i<size;i++)for(int j=k+1;j<size;j++){
            Int value=a[i][j]*d-a[i][k]*a[k][j];
            if(value%previous)throw runtime_error("Nonexact Bareiss division");
            a[i][j]=value/previous;
        }
        for(int i=k+1;i<size;i++)a[i][k]=0;
        previous=d;
    }
    return sign*a[size-1][size-1];
}
Mat A(const Mat& b){Mat a{};for(int i=0;i<n;i++)for(int j=0;j<n;j++)a[i][j]=2*(i==j)-b[i][j];return a;}
bool principal(const Mat& b,int j){
    for(int mask=1<<j;mask<(1<<(j+1));mask++){
        Mat a{};int r=0;
        for(int i=0;i<=j;i++)if(mask>>i&1){
            int c=0;for(int k=0;k<=j;k++)if(mask>>k&1)a[r][c++]=2*(i==k)-b[i][k];r++;
        }
        if(det(a,r)<=0)return false;
    }
    return true;
}
int weight(int i,int j){return int(max(p[i][j],m[i][j]));}
void path_bound(int at,int target,int visited,int length,long long product,int& limit){
    if(at==target&&length){limit=min(limit,int(((1LL<<(length+1))-1)/product));return;}
    for(int k=0;k<n;k++)if(!(visited>>k&1)&&at!=k&&weight(at,k))
        path_bound(k,target,visited|(1<<k),length+1,product*weight(at,k),limit);
}
bool cycles_at(int j){
    for(int k=0;k<j;k++)if(weight(j,k)){
        int limit=bound;path_bound(k,j,1<<k,0,1,limit);
        if(weight(j,k)>limit)return false;
    }
    return true;
}
bool parity(int i){Int x=0;for(int j=0;j<n;j++)x+=(2*(i==j)-p[i][j])*(2*(i==j)-m[i][j]);return x%2==0;}
bool strong(){
    bool reach[MAX_N][MAX_N]{};
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)reach[i][j]=i==j||weight(i,j);
    for(int k=0;k<n;k++)for(int i=0;i<n;i++)for(int j=0;j<n;j++)reach[i][j]|=reach[i][k]&&reach[k][j];
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(!reach[i][j])return false;
    return true;
}
vector<int> key(){
    vector<int> perm(n),best;iota(perm.begin(),perm.end(),0);
    do{for(int s=0;s<2;s++){
        vector<int> v;v.reserve(2*n*n);
        for(int sign=0;sign<2;sign++){const Mat& b=sign==s?p:m;for(int i:perm)for(int j:perm)v.push_back(int(b[i][j]));}
        if(best.empty()||v<best)best=move(v);
    }}while(next_permutation(perm.begin(),perm.end()));
    return best;
}
bool positive_ordered(){
    // The pointed cone's extreme rays span it. A strict positive vector exists
    // iff every strict inequality is positive on some feasible extreme ray.
    lp_checks++;
    vector<Vec> normals;
    for(int i=0;i<n;i++){Vec v{};v[i]=1;normals.push_back(v);}
    for(auto b:{p,m})for(int j=0;j<n;j++){Vec v{};for(int i=0;i<n;i++)v[i]=2*(i==j)-b[i][j];normals.push_back(v);}
    int strict=normals.size();
    for(int i=0;i<n-1;i++){Vec v{};v[i]=1;v[i+1]=-1;normals.push_back(v);}
    unsigned seen=0,all=(1U<<strict)-1;
    vector<int> chosen;
    function<bool(int)> visit=[&](int start){
        if(int(chosen.size())<n-1){
            for(int k=start;k<=int(normals.size())-(n-1-int(chosen.size()));k++){
                chosen.push_back(k);bool found=visit(k+1);chosen.pop_back();if(found)return true;
            }
            return false;
        }
        Vec ray{};
        for(int omit=0;omit<n;omit++){
            Mat minor{};
            for(int r=0;r<n-1;r++)for(int k=0,t=0;k<n;k++)if(k!=omit)minor[r][t++]=normals[chosen[r]][k];
            ray[omit]=(omit%2?-1:1)*det(minor,n-1);
        }
        Int g=0;for(int i=0;i<n;i++)g=gcd128(g,absolute(ray[i]));if(!g)return false;
        int orientation=0;for(int i=0;i<n;i++)if(ray[i]){orientation=ray[i]>0?1:-1;break;}
        for(auto& x:ray)x=x/g*orientation;
        for(int i=0;i<n;i++)if(ray[i]<0)return false;
        for(int i=0;i<n-1;i++)if(ray[i]<ray[i+1])return false;
        unsigned mask=0;
        for(int r=0;r<strict;r++){
            Int dot=0;for(int i=0;i<n;i++)dot+=normals[r][i]*ray[i];
            if(dot<0)return false;
            if(dot>0)mask|=1U<<r;
        }
        seen|=mask;return seen==all;
    };
    return visit(0);
}
void rows(int j);
struct Linear {Vec coefficients{};Int constant=0,lower=0,upper=0;};
Int floor_div(Int a,Int b){if(b<0){a=-a;b=-b;}Int q=a/b;if(a%b<0)q--;return q;}
Int ceil_div(Int a,Int b){return -floor_div(-a,b);}
void fill_row(int j,int k,const int* limits,const Mat& adj,Int denom,const vector<Linear>& constraints){
    if(k<j){
        Int low=0,high=limits[k];
        for(const auto& c:constraints){
            Int partial=c.constant,min_rest=0,max_rest=0;
            for(int i=0;i<k;i++)partial+=c.coefficients[i]*p[j][i];
            for(int i=k+1;i<j;i++){
                Int v=c.coefficients[i]*limits[i];min_rest+=min(Int(0),v);max_rest+=max(Int(0),v);
            }
            Int a=c.coefficients[k],lower=c.lower-partial-max_rest,upper=c.upper-partial-min_rest;
            if(a>0){low=max(low,ceil_div(lower,a));high=min(high,floor_div(upper,a));}
            else if(a<0){high=min(high,floor_div(lower,a));low=max(low,ceil_div(upper,a));}
            else if(lower>0||upper<0)return;
            if(low>high)return;
        }
        for(int x=int(low);x<=int(high);x++){p[j][k]=x;fill_row(j,k+1,limits,adj,denom,constraints);}
        p[j][k]=0;return;
    }
    nodes[j]++;
    if(!(nodes[j]&262143)&&seconds()-last_heartbeat>15){
        last_heartbeat=seconds();cerr<<"HEARTBEAT upper "<<current_upper<<" partner "<<current_partner<<" depth "<<j<<" nodes "<<nodes[j]<<" seconds "<<last_heartbeat<<endl;
    }
    if(!principal(p,j))return;
    Mat ap=A(p),am=A(m);Vec rhs{};
    for(int i=0;i<j;i++){
        for(int t=0;t<n;t++)rhs[i]+=am[i][t]*ap[j][t];
        for(int t=j;t<n;t++)rhs[i]-=ap[i][t]*am[j][t];
    }
    bool ok=true;
    for(int i=0;i<j;i++){
        Int x=0;for(int t=0;t<j;t++)x-=adj[i][t]*rhs[t];
        if(x%denom||x<0||x>limits[i]*denom){ok=false;break;}
        m[j][i]=x/denom;
    }
    if(ok&&parity(j)&&cycles_at(j)&&principal(m,j))rows(j+1);
    for(int i=0;i<j;i++)m[j][i]=0;
}
void rows(int j){
    if(j==n){
        if(!strong())return;
        labelled++;auto v=key();
        if(answers.count(v))return;
        if(!positive_cache.count(v)&&!positive_ordered())return;
        answers.insert(v);positive_cache.insert(move(v));return;
    }
    int limits[MAX_N]{};
    for(int k=0;k<j;k++){limits[k]=ordered_limits[j][k];path_bound(k,j,1<<k,0,1,limits[k]);}
    Mat a=A(p),adj{};Int denom=det(a,j);if(denom<=0)throw runtime_error("Nonpositive pivot block");
    for(int i=0;i<j;i++)for(int t=0;t<j;t++){
        Mat minor{};
        for(int r=0,rr=0;r<j;r++)if(r!=t){for(int c=0,cc=0;c<j;c++)if(c!=i)minor[rr][cc++]=a[r][c];rr++;}
        adj[i][t]=((i+t)%2?-1:1)*det(minor,j-1);
    }
    // The solved M row is affine-linear in x_k=P[j][k]. Enforce its bounds
    // on partial assignments, before constructing any complete candidate.
    Mat am=A(m),beta{};Vec alpha{},base_rhs{};
    for(int i=0;i<j;i++)for(int t=j;t<n;t++)base_rhs[i]+=am[i][t]*a[j][t]-a[i][t]*am[j][t];
    vector<Linear> constraints;
    for(int i=0;i<j;i++){
        for(int t=0;t<j;t++){
            alpha[i]-=adj[i][t]*base_rhs[t];
            for(int k=0;k<j;k++)beta[i][k]+=adj[i][t]*am[t][k];
        }
        constraints.push_back({beta[i],alpha[i],0,limits[i]*denom});
    }
    Linear leading;leading.constant=(2-p[j][j])*denom;leading.lower=1;
    for(int k=0;k<j;k++)for(int t=0;t<j;t++)leading.coefficients[k]+=adj[k][t]*a[t][j];
    leading.upper=leading.constant;
    for(int k=0;k<j;k++)leading.upper+=max(Int(0),leading.coefficients[k]*limits[k]);
    constraints.push_back(leading);
    fill_row(j,0,limits,adj,denom,constraints);
}
void triangular(int j,Mat b,vector<Mat>& out){
    if(j==n){out.push_back(b);return;}
    for(int diag=0;diag<2;diag++){
        b[j][j]=diag;triangular(j+1,b,out);
        if(!diag)for(int i=0;i<j;i++){b[i][j]=1;triangular(j+1,b,out);b[i][j]=0;}
    }
}
int main(int argc,char** argv){
    if(argc<5)return 2;
    n=stoi(argv[1]);if(n<2||n>MAX_N)return 2;
    bound=(1<<n)-1;fs::path directory=argv[2];fs::create_directories(directory);
    int shard=stoi(argv[3]),shards=stoi(argv[4]),only=argc>5?stoi(argv[5]):-1;
    if(shards<1||shard<0||shard>=shards)return 2;
    vector<Mat> triangle;triangular(0,Mat{},triangle);
    for(int a=0;a<int(triangle.size());a++){
        if(a%shards!=shard||(only>=0&&a!=only))continue;
        auto path=directory/("upper-"+to_string(a)+".json");
        if(fs::exists(path)){cout<<"SKIP "<<a<<endl;continue;}
        answers.clear();double before=seconds();
        current_upper=a;
        auto labelled_before=labelled,lp_before=lp_checks;
        array<unsigned long long,MAX_N+1> nodes_before{};copy(begin(nodes),end(nodes),nodes_before.begin());
        for(int b=a;b<int(triangle.size());b++){
            current_partner=b;
            p=triangle[a];m=triangle[b];if(!parity(0))continue;
            bool ok=true;
            for(int cut=1;cut<n;cut++){
                bool crossing=false;for(int i=0;i<cut;i++)for(int j=cut;j<n;j++)crossing|=weight(i,j)>0;
                ok&=crossing;
            }
            if(ok){
                // v_0 >= ... >= v_{n-1}. Moving backwards in the index order
                // costs no factor; each known upper arrow gives v_i < 2 v_j.
                // A path from k to j with d upper arrows bounds every lower
                // weight C[j][k] by 2^(d+1)-1, often much below 2^n-1.
                int distance[MAX_N][MAX_N];
                for(int i=0;i<n;i++)for(int j=0;j<n;j++)distance[i][j]=i>=j?0:(weight(i,j)?1:100);
                for(int k=0;k<n;k++)for(int i=0;i<n;i++)for(int j=0;j<n;j++)distance[i][j]=min(distance[i][j],distance[i][k]+distance[k][j]);
                for(int i=0;i<n;i++)for(int j=0;j<i;j++){
                    if(distance[j][i]>=n)throw runtime_error("Disconnected ordered upper graph");
                    ordered_limits[i][j]=(1<<(distance[j][i]+1))-1;
                }
                rows(1);
            }
        }
        auto temporary=path;temporary+=".tmp";ofstream out(temporary);
        out<<"{\"rank\":"<<n<<",\"upper_index\":"<<a<<",\"triangular_count\":"<<triangle.size()
           <<",\"completed\":true,\"parity_pruning\":true,\"seconds\":"<<seconds()-before
           <<",\"labelled\":"<<labelled-labelled_before<<",\"positivity_checks\":"<<lp_checks-lp_before<<",\"row_nodes\":[";
        for(int j=1;j<n;j++){if(j>1)out<<",";out<<nodes[j]-nodes_before[j];}out<<"],\"keys\":[";
        bool first=true;for(auto v:answers){if(!first)out<<",";first=false;out<<"[";for(size_t i=0;i<v.size();i++){if(i)out<<",";out<<v[i];}out<<"]";}
        out<<"]}\n";out.close();fs::rename(temporary,path);
        cout<<"upper "<<a<<"/"<<triangle.size()<<" orbits "<<answers.size()<<" nodes ";
        for(int j=1;j<n;j++)cout<<nodes[j]-nodes_before[j]<<" ";
        cout<<" task_seconds "<<seconds()-before<<" total_seconds "<<seconds()<<endl;
    }
    cout<<"DONE rank "<<n<<" shard "<<shard<<"/"<<shards<<" seconds "<<seconds()<<endl;
}
