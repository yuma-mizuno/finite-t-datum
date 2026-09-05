// Independent exhaustive check of every integer pulse position up to the
// proved bound. Uses the A matrices directly, without the Python L pruning.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
struct Term{int delay,row,column,plus,minus;};
int main(int argc,char**argv){
    if(argc!=3)return 2;
    ifstream in(argv[1]);int n,h,a,diagonal,attachment_delay;in>>n>>h>>a>>diagonal>>attachment_delay;
    vector<vector<int>> target(2,vector<int>(n));int mass=0;
    for(auto& v:target)for(int& x:v){in>>x;mass+=x;}
    vector<Term> terms;
    for(int d=1;d<=h;d++)for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        int p,m;in>>p>>m;if(p||m)terms.push_back({d,i,j,p,m});
    }
    int bound=(mass+1)*h;unsigned long long tested=0,solutions=0;auto start=chrono::steady_clock::now();
    for(int R=1;R<=bound;R++)for(int L=diagonal?1:0;L<(diagonal?R:1);L++){
        tested++;vector<vector<int>> w(R+h+1,vector<int>(n));vector<vector<int>> counts(2,vector<int>(n));bool bad=false;
        for(int t=0;t<=R+h&&!bad;t++){
            if(t==0)w[t][a]--;
            if(diagonal&&t==L)w[t][a]++;
            if(t==R)w[t][a]--;
            for(const auto& v:terms)if(t>=v.delay){
                int old=w[t-v.delay][v.column];w[t][v.row]+=-v.minus*max(old,0)+v.plus*max(-old,0);
            }
            for(int i=0;i<n;i++){
                int x=w[t][i];int sign=x<0?1:0;counts[sign][i]+=abs(x);
                if(counts[sign][i]>target[sign][i]||(t>R-attachment_delay&&x)){bad=true;break;}
            }
        }
        if(!bad&&counts==target)solutions++;
    }
    ofstream out(argv[2]);out<<"{\"delay_bound\":"<<bound<<",\"pulse_positions_checked\":"<<tested<<",\"solutions\":"<<solutions
       <<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
