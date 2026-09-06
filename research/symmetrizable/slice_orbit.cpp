// Exhaust the same canonical rotations and commuting moves as rank4_slices.py.
// No orbit-size cutoff. Compact exact skew matrices; explicit overflow failure.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
using Matrix=vector<vector<int>>;
int n,r;
vector<int> weights;
void put(string& s,int x){if(x<0||x>65535)throw runtime_error("Encoding overflow");s.push_back(char(x>>8));s.push_back(char(x&255));}
int get(const string& s,int& at){int a=(unsigned char)s[at++];return (a<<8)|(unsigned char)s[at++];}
string canonical(const Matrix& b,const vector<int>& word,const vector<int>& p){
    vector<int> order,lengths;
    for(int root:word){int k=root,length=0;do{order.push_back(k);length++;k=p[k];if(length>n)throw runtime_error("Invalid permutation");}while(k!=root);lengths.push_back(length);}
    vector<int> check=order;sort(check.begin(),check.end());for(int i=0;i<n;i++)if(int(check.size())!=n||check[i]!=i)throw runtime_error("Invalid roots");
    string s;s.reserve(2*r+n*(n-1));for(int length:lengths)put(s,length);for(int root:word)put(s,weights[root]);
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)put(s,b[order[i]][order[j]]+32768);
    return s;
}
void decode(const string& s,Matrix& b,vector<int>& word,vector<int>& p,vector<int>& lengths){
    int at=0,base=0;word.clear();lengths.clear();p.assign(n,0);b.assign(n,vector<int>(n));
    for(int i=0;i<r;i++){int length=get(s,at);lengths.push_back(length);word.push_back(base);for(int j=0;j<length;j++)p[base+j]=base+(j+1)%length;base+=length;}
    weights.clear();for(int length:lengths){int d=get(s,at);if(d<=0)throw runtime_error("Nonpositive symmetrizer");for(int k=0;k<length;k++)weights.push_back(d);}
    if(base!=n)throw runtime_error("Bad cycle lengths");
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){b[i][j]=get(s,at)-32768;if((int64_t(b[i][j])*weights[j])%weights[i])throw runtime_error("Nonintegral opposite entry");b[j][i]=-int64_t(b[i][j])*weights[j]/weights[i];}
}
Matrix mutate(const Matrix& b,int k){
    Matrix c=b;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        int64_t x=(i==k||j==k)?-int64_t(b[i][j]):int64_t(b[i][j])+int64_t(max(b[i][k],0))*max(b[k][j],0)-int64_t(max(-b[i][k],0))*max(-b[k][j],0);
        if(x<=-32768||x>=32768)throw runtime_error("Matrix weight exceeds compact exact encoding");
        c[i][j]=int(x);
    }
    return c;
}
int main(int argc,char** argv){
    if(argc!=3)return 2;
    auto started=chrono::steady_clock::now();ifstream input(argv[1]);if(!(input>>r>>n)||r<1||n<r)throw runtime_error("Invalid input");
    Matrix b(n,vector<int>(n));vector<int> word(r),p(n);
    for(auto& row:b)for(auto& x:row)input>>x;
    for(auto& x:word)input>>x;
    for(auto& x:p)input>>x;
    weights.resize(n);for(auto& d:weights)input>>d;
    for(int i=0;i<n;i++){if(weights[i]<=0||weights[i]!=weights[p[i]])throw runtime_error("Invalid symmetrizer");for(int j=0;j<n;j++)if(int64_t(b[i][j])*weights[j]!=-int64_t(b[j][i])*weights[i])throw runtime_error("Not skew symmetrizable");}
    unordered_set<string> seen;deque<const string*> pending;string minimum;
    auto insert=[&](string s){auto [it,fresh]=seen.insert(move(s));if(fresh){pending.push_back(&*it);if(minimum.empty()||*it<minimum)minimum=*it;}};
    insert(canonical(b,word,p));for(auto& row:b)for(auto& x:row)x=-x;insert(canonical(b,word,p));
    size_t processed=0;double last=0;
    while(!pending.empty()){
        const string* s=pending.front();pending.pop_front();vector<int> lengths;
        decode(*s,b,word,p,lengths);vector<int> inverse(n);for(int i=0;i<n;i++)inverse[p[i]]=i;
        vector<int> rotated(word.begin()+1,word.end());rotated.push_back(inverse[word[0]]);
        insert(canonical(mutate(b,word[0]),rotated,p));
        Matrix prefix=b;
        for(int k=0;k<r-1;k++){
            if(prefix[word[k]][word[k+1]]==0){auto other=word;swap(other[k],other[k+1]);insert(canonical(b,other,p));}
            prefix=mutate(prefix,word[k]);
        }
        processed++;
        if(!(processed%10000)){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-started).count();
            if(elapsed-last>=15){cerr<<"Slice orbit: "<<processed<<" processed, "<<pending.size()<<" pending; "<<elapsed<<" seconds"<<endl;last=elapsed;}
        }
    }
    vector<int> lengths;decode(minimum,b,word,p,lengths);
    ofstream out(argv[2]);out<<"{\"states\":"<<seen.size()<<",\"processed\":"<<processed<<",\"complete\":true,\"minimum\":[";
    bool first=true;for(int x:lengths){if(!first)out<<',';first=false;out<<x;}for(int root:word)out<<','<<weights[root];for(const auto& row:b)for(int x:row)out<<','<<x;
    out<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
}
