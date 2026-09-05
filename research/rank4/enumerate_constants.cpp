// Exhaustive integer enumeration after ordering a simultaneous positive vector.
// Compile: g++ -O3 -std=c++17 enumerate_constants.cpp -o bin/enumerate_constants
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using namespace std;
using Mat = array<array<long long,4>,4>;
int n=4, bound=15;
Mat p{},m{};
set<vector<int>> answers;
long long nodes[5]{}, labelled=0, positive_count=0;

long long det(Mat a,int size) {
    if(size==0) return 1;
    if(size==1) return a[0][0];
    long long d=0;
    for(int j=0;j<size;j++) {
        Mat b{};
        for(int i=1;i<size;i++) for(int k=0,t=0;k<size;k++)
            if(k!=j) b[i-1][t++]=a[i][k];
        d+=(j%2?-1:1)*a[0][j]*det(b,size-1);
    }
    return d;
}
Mat A(const Mat& b) {
    Mat a{};
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) a[i][j]=2*(i==j)-b[i][j];
    return a;
}
bool principal(const Mat& b,int j) {
    for(int mask=1<<j;mask<(1<<(j+1));mask++) {
        Mat a{}; int r=0;
        for(int i=0;i<=j;i++) if(mask>>i&1) {
            int c=0;
            for(int k=0;k<=j;k++) if(mask>>k&1) a[r][c++]=2*(i==k)-b[i][k];
            r++;
        }
        if(det(a,r)<=0) return false;
    }
    return true;
}
int max_weight(int i,int j) {return max(p[i][j],m[i][j]);}
void path_bound(int start,int at,int target,int visited,int length,long long product,int& limit) {
    if(at==target && length) {limit=min(limit,int(((1<<(length+1))-1)/product));return;}
    for(int k=0;k<n;k++) if(!(visited>>k&1) && at!=k && max_weight(at,k))
        path_bound(start,k,target,visited|(1<<k),length+1,product*max_weight(at,k),limit);
}
bool cycles_at(int j) {
    // Every newly completed cycle contains j, because later lower rows are zero.
    for(int k=0;k<j;k++) if(max_weight(j,k)) {
        int limit=bound;
        path_bound(k,k,j,1<<k,0,1,limit);
        if(max_weight(j,k)>limit) return false;
    }
    return true;
}
bool strong() {
    bool reach[4][4]{};
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) reach[i][j]=i==j || max_weight(i,j);
    for(int k=0;k<n;k++) for(int i=0;i<n;i++) for(int j=0;j<n;j++) reach[i][j]|=reach[i][k]&&reach[k][j];
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(!reach[i][j]) return false;
    return true;
}
vector<int> key() {
    vector<int> perm(n),best;
    iota(perm.begin(),perm.end(),0);
    do { for(int s=0;s<2;s++) {
        vector<int> v;
        for(int sign=0;sign<2;sign++) {
            const Mat& b=(sign==s?p:m);
            for(int i:perm) for(int j:perm) v.push_back(b[i][j]);
        }
        if(best.empty() || v<best) best=v;
    }} while(next_permutation(perm.begin(),perm.end()));
    return best;
}
bool positive_ordered() {
    // Extreme rays of {v>=0, vA+>=0, vA->=0, v_i>=v_{i+1}}.
    vector<array<long long,4>> normals;
    for(int i=0;i<n;i++) {array<long long,4> v{};v[i]=1;normals.push_back(v);}
    for(auto b:{p,m}) for(int j=0;j<n;j++) {
        array<long long,4> v{};for(int i=0;i<n;i++) v[i]=2*(i==j)-b[i][j];normals.push_back(v);
    }
    int strict=normals.size();
    for(int i=0;i<n-1;i++) {array<long long,4> v{};v[i]=1;v[i+1]=-1;normals.push_back(v);}
    array<long long,4> sum{};
    for(int a=0;a<int(normals.size());a++) for(int b=a+1;b<int(normals.size());b++)
    for(int c=(n==4?b+1:b); c<(n==4?int(normals.size()):b+1);c++) {
        array<long long,4> ray{};
        for(int omit=0;omit<n;omit++) {
            Mat minor{};int rs[3]={a,b,c};
            for(int r=0;r<n-1;r++) for(int k=0,t=0;k<n;k++) if(k!=omit) minor[r][t++]=normals[rs[r]][k];
            ray[omit]=(omit%2?-1:1)*det(minor,n-1);
        }
        long long g=0;for(int i=0;i<n;i++) g=gcd(g,abs(ray[i]));
        if(!g) continue;
        for(auto& x:ray) x/=g;
        for(int sign:{-1,1}) {
            bool ok=true;
            for(auto v:normals) {long long d=0;for(int i=0;i<n;i++) d+=v[i]*ray[i]*sign;if(d<0){ok=false;break;}}
            if(ok) for(int i=0;i<n;i++) sum[i]+=sign*ray[i];
        }
    }
    for(int r=0;r<strict;r++) {long long d=0;for(int i=0;i<n;i++) d+=normals[r][i]*sum[i];if(d<=0)return false;}
    return true;
}

void rows(int j);
void fill_row(int j,int k,const int* limits,const Mat& adj,long long denom) {
    if(k<j) {
        for(int x=0;x<=limits[k];x++) {p[j][k]=x;fill_row(j,k+1,limits,adj,denom);}
        p[j][k]=0;return;
    }
    nodes[j]++;
    if(!principal(p,j))return;
    // A+_{<j,<j} A-_{j,<j}^t = A-_{<j,*} A+_{j,*}^t - A+_{<j,>=j} A-_{j,>=j}^t.
    Mat ap=A(p),am=A(m);
    long long rhs[4]{};
    for(int i=0;i<j;i++) {
        for(int t=0;t<n;t++) rhs[i]+=am[i][t]*ap[j][t];
        for(int t=j;t<n;t++) rhs[i]-=ap[i][t]*am[j][t];
    }
    bool ok=true;
    for(int i=0;i<j;i++) {
        long long x=0;for(int t=0;t<j;t++) x-=adj[i][t]*rhs[t];
        if(x%denom || x<0 || x>limits[i]*denom) {ok=false;break;}
        m[j][i]=x/denom;
    }
    if(ok && cycles_at(j) && principal(m,j)) rows(j+1);
    for(int i=0;i<j;i++)m[j][i]=0;
}
void rows(int j) {
    if(j==n) {
        if(!strong())return;
        labelled++;
        auto v=key();
        if(answers.count(v)) return;
        if(!positive_ordered())return;
        positive_count++;answers.insert(v);return;
    }
    int limits[4]{};
    for(int k=0;k<j;k++) {limits[k]=bound;path_bound(k,k,j,1<<k,0,1,limits[k]);}
    Mat a=A(p),adj{};
    long long denom=det(a,j);
    if(denom<=0) abort();
    for(int i=0;i<j;i++)for(int t=0;t<j;t++) {
        Mat minor{};
        for(int r=0,rr=0;r<j;r++)if(r!=t) {
            for(int c=0,cc=0;c<j;c++)if(c!=i)minor[rr][cc++]=a[r][c];
            rr++;
        }
        adj[i][t]=((i+t)%2?-1:1)*det(minor,j-1);
    }
    fill_row(j,0,limits,adj,denom);
}
void triangular(int j,Mat b,vector<Mat>& out) {
    if(j==n){out.push_back(b);return;}
    for(int diag=0;diag<2;diag++) {
        b[j][j]=diag;triangular(j+1,b,out);
        if(!diag)for(int i=0;i<j;i++) {b[i][j]=1;triangular(j+1,b,out);b[i][j]=0;}
    }
}
int main(int argc,char** argv) {
    if(argc>1)n=stoi(argv[1]);
    if(n!=3 && n!=4)return 2;
    bound=(1<<n)-1;
    int shard=argc>2?stoi(argv[2]):0, shards=argc>3?stoi(argv[3]):1;
    string file=argc>4?argv[4]:"constant_keys.txt";
    vector<Mat> triangle;triangular(0,Mat{},triangle);
    int tasks=0;
    auto start=chrono::steady_clock::now();
    for(int a=0;a<int(triangle.size());a++) {
        if(a%shards!=shard)continue;
        for(int b=a;b<int(triangle.size());b++) {
            p=triangle[a];m=triangle[b];
            bool ok=true;
            for(int cut=1;cut<n;cut++) {
                bool crossing=false;
                for(int i=0;i<cut;i++)for(int j=cut;j<n;j++)crossing|=max_weight(i,j)>0;
                ok &= crossing;
            }
            if(ok)rows(1);
            tasks++;
        }
        cerr<<"upper "<<a<<"/"<<triangle.size()<<" tasks "<<tasks<<" orbits "<<answers.size()<<" nodes ";
        for(int j=1;j<n;j++)cerr<<nodes[j]<<" ";
        cerr<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
    }
    ofstream out(file);
    for(auto v:answers) {for(size_t i=0;i<v.size();i++){if(i)out<<" ";out<<v[i];}out<<"\n";}
    cout<<"rank "<<n<<" triangular "<<triangle.size()<<" labelled "<<labelled<<" orbits "<<answers.size()<<endl;
}
