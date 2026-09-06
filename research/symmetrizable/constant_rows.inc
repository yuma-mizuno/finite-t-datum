// Integer row reduction of the linear symplectic equations in d_0,...,d_n-1.
// Each row is primitive, with positive first entry. No floating point is used.
using Equations=vector<Vec>;
Equations inherited_equations[MAX_N+1];
unordered_map<string,Vec> principal_diagonals[MAX_N];
bool weighted_principal_pruning=false;
void load_principal_diagonals(const string& path){
    ifstream input(path);if(!input)return;
    int size,value;
    while(input>>size){
        if(size<1||size>=MAX_N)throw runtime_error("Invalid weighted principal rank");
        Vec d{};for(int i=0;i<size;i++){if(!(input>>value)||value<=0)throw runtime_error("Invalid principal weight");d[i]=value;}
        string key;for(int i=0;i<2*size*size;i++){if(!(input>>value)||value<0||value>63)throw runtime_error("Invalid weighted principal entry");key.push_back(char(value));}
        if(principal_diagonals[size].count(key)&&principal_diagonals[size][key]!=d)throw runtime_error("Principal symmetrizer is not unique");
        principal_diagonals[size][key]=d;
    }
    for(int size=1;size<n;size++)for(const auto& key:principal_constants[size])if(!principal_diagonals[size].count(key))throw runtime_error("Incomplete weighted principal table");
    weighted_principal_pruning=true;
}
unordered_map<string,vector<string>> completion_tables[MAX_N];
void build_completion_tables(){
    for(int size=2;size<n;size++){
        string full(2*size*size,0);
        function<void(int)> fill=[&](int remaining){
            if(!remaining){
                string known,row;
                for(int sign=0;sign<2;sign++){
                    for(int i=0;i<size-1;i++)for(int j=0;j<size;j++)known.push_back(full[sign*size*size+i*size+j]);
                    known.push_back(full[sign*size*size+size*size-1]);
                    for(int j=0;j<size-1;j++)row.push_back(full[sign*size*size+(size-1)*size+j]);
                }
                completion_tables[size][known].push_back(row);return;
            }
            int first=remaining&-remaining,rest=remaining^first;
            for(int subset=rest;;subset=(subset-1)&rest){
                int mask=subset|first;vector<int> vertices;
                for(int i=0;i<size;i++)if(mask>>i&1)vertices.push_back(i);
                int count=vertices.size();
                for(const auto& key:principal_constants[count]){
                    int at=0;
                    for(int sign=0;sign<2;sign++)for(int i:vertices)for(int j:vertices)full[sign*size*size+i*size+j]=key[at++];
                    fill(remaining^mask);
                }
                for(int sign=0;sign<2;sign++)for(int i:vertices)for(int j:vertices)full[sign*size*size+i*size+j]=0;
                if(!subset)break;
            }
        };
        fill((1<<size)-1);
    }
}
bool normalize(Vec& v){
    Int g=0;int first=-1;
    for(int i=0;i<n;i++)if(v[i]){g=gcd128(g,absolute(v[i]));if(first<0)first=i;}
    if(first<0)return false;
    if(v[first]<0)g=-g;
    for(int i=0;i<n;i++)v[i]/=g;
    return true;
}
int pivot(const Vec& v){for(int i=0;i<n;i++)if(v[i])return i;return n;}
bool add_equation(Equations& eq,Vec v){
    for(const auto& row:eq){int k=pivot(row);if(v[k]){
        Int a=row[k],b=v[k],g=gcd128(absolute(a),absolute(b));a/=g;b/=g;
        for(int i=0;i<n;i++)v[i]=a*v[i]-b*row[i];
        normalize(v);
    }}
    if(!normalize(v))return false;
    int k=pivot(v);
    for(auto& row:eq)if(row[k]){
        Int a=v[k],b=row[k],g=gcd128(absolute(a),absolute(b));a/=g;b/=g;
        for(int i=0;i<n;i++)row[i]=a*row[i]-b*v[i];
        normalize(row);
    }
    eq.push_back(v);sort(eq.begin(),eq.end(),[](const Vec&a,const Vec&b){return pivot(a)<pivot(b);});
    return true;
}
bool possible_positive(const Equations& eq){
    if(int(eq.size())==n)return false;
    for(const auto& row:eq){bool positive=false,negative=false;
        for(int i=0;i<n;i++){positive|=row[i]>0;negative|=row[i]<0;}
        if(!positive||!negative)return false;
    }
    return true;
}
Vec unique_diagonal(const Equations& eq){
    if(int(eq.size())!=n-1)throw runtime_error("Symmetrizer is not unique");
    bool used[MAX_N]{};for(const auto& row:eq)used[pivot(row)]=true;
    int free=0;while(used[free])free++;
    Int common=1;for(const auto& row:eq){Int a=row[pivot(row)];common=common/gcd128(common,a)*a;}
    Vec d{};d[free]=common;
    for(const auto& row:eq)d[pivot(row)]=-row[free]*(common/row[pivot(row)]);
    normalize(d);return d;
}
Equations diagonal_equations(int size){
    Equations eq;
    for(int i=0;i<size;i++)for(int j=i+1;j<size;j++){
        Vec row{};
        for(int k=0;k<n;k++)row[k]=(2*(i==k)-p[i][k])*(2*(j==k)-m[j][k])-(2*(i==k)-m[i][k])*(2*(j==k)-p[j][k]);
        add_equation(eq,row);
    }
    return eq;
}
bool valid_diagonal(const Vec& d,int size){
    for(int i=0;i<n;i++)if(d[i]<=0)return false;
    for(int i=0;i<size;i++)for(int j=0;j<n;j++)for(auto& b:{p,m}){
        Int v=b[i][j]*d[j];if(v%d[i]||v/d[i]>bound)return false;
    }
    for(int i=0;i<size;i++){
        Int a=0;for(int j=0;j<n;j++)a+=(2*(i==j)-p[i][j])*(2*(i==j)-m[i][j])*d[j];
        if(a%(2*d[i]))return false;
    }
    return true;
}
vector<int> weighted_key(const Vec& d){
    vector<int> perm(n),best;iota(perm.begin(),perm.end(),0);
    do{for(int s=0;s<2;s++){
        vector<int> v;v.reserve(n+2*n*n);
        for(int i:perm){if(d[i]>2147483647)throw runtime_error("Symmetrizer key overflow");v.push_back(int(d[i]));}
        for(int sign=0;sign<2;sign++){const Mat& b=sign==s?p:m;for(int i:perm)for(int j:perm)v.push_back(int(b[i][j]));}
        if(best.empty()||v<best)best=move(v);
    }}while(next_permutation(perm.begin(),perm.end()));
    return best;
}
vector<int> ratio_choices(int i,int k);
bool partial_diagonal(int completed,const Equations& eq);
void emit_diagonals(Equations eq){
    if(!possible_positive(eq)||!partial_diagonal(n,eq))return;
    if(int(eq.size())==n-1){
        Vec d=unique_diagonal(eq);if(valid_diagonal(d,n))answers.insert(weighted_key(d));return;
    }
    // Strong connectivity ensures that some edge ratio is not yet fixed.
    // Its integral dual coefficient is in [1,2^n-1]. Each branch raises rank.
    int from=-1,to=-1;vector<int> choices;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j&&weight(i,j)){
        Vec relation{};relation[j]=weight(i,j);relation[i]=-1;
        Equations probe=eq;
        if(!add_equation(probe,relation))continue;
        // Test whether d_i and d_j are already proportional on the kernel.
        Vec ei{},ej{};ei[i]=1;ej[j]=1;auto a=eq;add_equation(a,ei);int rank=a.size();add_equation(a,ej);
        if(int(a.size())==rank)continue;
        auto values=ratio_choices(i,j);
        if(from<0||values.size()<choices.size()){from=i;to=j;choices=move(values);}
    }
    if(from>=0){for(int q:choices){auto branch=eq;Vec relation{};relation[to]=weight(from,to);relation[from]=-q;add_equation(branch,relation);emit_diagonals(move(branch));}return;}
    throw runtime_error("No independent ratio in a strongly connected graph");
}
void rows(int j);
struct Linear {Vec coefficients{};Int constant=0,lower=0,upper=0;};
Int floor_div(Int a,Int b){if(b<0){a=-a;b=-b;}Int q=a/b;if(a%b<0)q--;return q;}
Int ceil_div(Int a,Int b){return -floor_div(-a,b);}
void row_bounds(int j,int* limits,int* m_limits){
    for(int k=0;k<j;k++){
        int cap=ordered_limits[j][k];path_bound(k,j,1<<k,0,1,cap);
        auto column_cap=[&](const Mat& b){
            Int used=0;for(int i=0;i<n;i++)if(i!=j)used+=b[i][k]*(1<<(n-ordered_distance[k][i]));
            return min(cap,int(floor_div((2<<n)-1-used,1<<(n-ordered_distance[k][j]))));
        };
        limits[k]=column_cap(p);m_limits[k]=column_cap(m);
    }
}
struct CommonBudget {Vec plus{},minus{},maximum_entry{};Int maximum=0;};
vector<CommonBudget> common_budgets[MAX_N];
Int common_remaining[MAX_N][MAX_N+1][1<<MAX_N];
bool prepare_common_budgets(int j){
    auto& budgets=common_budgets[j];budgets.clear();
    // A common positive left vector makes every independent column choice
    // from A_plus and A_minus a nonsingular M-matrix. Enforce every new
    // mixed-column Schur complement, not just the two unmixed matrices.
    int neutral=0;
    for(int k=0;k<j;k++){bool same=true;for(int i=0;i<j;i++)same&=p[i][k]==m[i][k];if(same)neutral|=1<<k;}
    {bool same=true;for(int i=0;i<=j;i++)same&=p[i][j]==m[i][j];if(same)neutral|=1<<j;}
    for(int mask=0;mask<(1<<(j+1));mask++){
        if(mask&neutral)continue;
        Mat a{};for(int i=0;i<j;i++)for(int k=0;k<j;k++)a[i][k]=2*(i==k)-((mask>>k&1)?m[i][k]:p[i][k]);
        Int denominator=det(a,j);if(denominator<=0)return false;
        const Mat& top=(mask>>j&1)?m:p;CommonBudget budget;budget.maximum=(2-top[j][j])*denominator-1;
        for(int i=0;i<j;i++){
            Int coefficient=0;
            for(int t=0;t<j;t++){
                Mat minor{};for(int r=0,rr=0;r<j;r++)if(r!=t){for(int k=0,kk=0;k<j;k++)if(k!=i)minor[rr][kk++]=a[r][k];rr++;}
                coefficient+=((i+t)%2?-1:1)*det(minor,j-1)*top[t][j];
            }
            if(coefficient<0)throw runtime_error("Negative mixed M-matrix adjugate coefficient");
            if(neutral>>i&1)budget.maximum_entry[i]=coefficient;
            else ((mask>>i&1)?budget.minus:budget.plus)[i]=coefficient;
        }
        bool duplicate=false;for(const auto& previous:budgets)if(previous.plus==budget.plus&&previous.minus==budget.minus&&previous.maximum_entry==budget.maximum_entry&&previous.maximum==budget.maximum){duplicate=true;break;}
        if(!duplicate)budgets.push_back(budget);
    }
    for(int i=0;i<int(budgets.size());i++)common_remaining[j][0][i]=budgets[i].maximum;
    return true;
}
bool common_row_valid(int j){
    for(const auto& b:common_budgets[j]){
        Int used=0;for(int i=0;i<j;i++)used+=b.plus[i]*p[j][i]+b.minus[i]*m[j][i]+b.maximum_entry[i]*max(p[j][i],m[j][i]);
        if(used>b.maximum)return false;
    }
    return true;
}
unsigned row_positive[MAX_N][MAX_N+1]{},row_negative[MAX_N][MAX_N+1]{};
unsigned future_positive[MAX_N][MAX_N+1]{},future_negative[MAX_N][MAX_N+1]{};
void free_row(int j,int k,const int* limits,const int* ml){
    if(k==0){
        row_positive[j][0]=row_negative[j][0]=0;
        for(int i=0;i<j;i++)for(int t=j;t<n;t++){
            Int c=(2*(i==t)-p[i][t])*(2*(j==t)-m[j][t])-(2*(i==t)-m[i][t])*(2*(j==t)-p[j][t]);
            if(c>0)row_positive[j][0]|=1<<i;if(c<0)row_negative[j][0]|=1<<i;
        }
        future_positive[j][j]=future_negative[j][j]=0;
        for(int t=j-1;t>=0;t--){
            future_positive[j][t]=future_positive[j][t+1];future_negative[j][t]=future_negative[j][t+1];
            for(int i=0;i<j;i++){
                Int a=(2*(i==t)-m[i][t])*limits[t],b=-(2*(i==t)-p[i][t])*ml[t];
                if(a>0||b>0)future_positive[j][t]|=1<<i;
                if(a<0||b<0)future_negative[j][t]|=1<<i;
            }
        }
    }
    if((row_positive[j][k]&~row_negative[j][k]&~future_negative[j][k])||
       (row_negative[j][k]&~row_positive[j][k]&~future_positive[j][k]))return;

    if(k<j){
        int caps[2]={limits[k],ml[k]};
        for(int i=0;i<int(common_budgets[j].size());i++){
            const auto& b=common_budgets[j][i];Int remaining=common_remaining[j][k][i];
            if(remaining<0)return;
            if(b.plus[k])caps[0]=min(caps[0],int(remaining/b.plus[k]));
            if(b.minus[k])caps[1]=min(caps[1],int(remaining/b.minus[k]));
            if(b.maximum_entry[k]){int cap=int(remaining/b.maximum_entry[k]);caps[0]=min(caps[0],cap);caps[1]=min(caps[1],cap);}
        }
        for(int a=0;a<=caps[0];a++)for(int b=0;b<=caps[1];b++){
            p[j][k]=a;m[j][k]=b;
            row_positive[j][k+1]=row_positive[j][k];row_negative[j][k+1]=row_negative[j][k];
            for(int i=0;i<j;i++){
                Int c=(2*(i==k)-m[i][k])*a-(2*(i==k)-p[i][k])*b;
                if(c>0)row_positive[j][k+1]|=1<<i;if(c<0)row_negative[j][k+1]|=1<<i;
            }

            for(int i=0;i<int(common_budgets[j].size());i++)common_remaining[j][k+1][i]=common_remaining[j][k][i]-common_budgets[j][i].plus[k]*a-common_budgets[j][i].minus[k]*b-common_budgets[j][i].maximum_entry[k]*max(a,b);
            free_row(j,k+1,limits,ml);
        }
        p[j][k]=m[j][k]=0;return;
    }
    nodes[j]++;
    if(!(nodes[j]&262143)&&seconds()-last_heartbeat>15){last_heartbeat=seconds();cerr<<"HEARTBEAT upper "<<current_upper<<" partner "<<current_partner<<" depth "<<j<<" nodes "<<nodes[j]<<" seconds "<<last_heartbeat<<endl;}
    // The two exact Schur budgets above already prove both leading minors.
    rows(j+1);
}
Vec fixed_minus_coefficients[MAX_N];Int fixed_minus_constant[MAX_N];
void fixed_row(int j,int k,const int* limits,const int* ml,const Mat& adj,Int det_a,const vector<Linear>& cs,const Vec& d){
    if(k<j){
        Int low=0,high=limits[k];
        for(const auto& c:cs){
            Int partial=c.constant,lo=0,hi=0;for(int i=0;i<k;i++)partial+=c.coefficients[i]*p[j][i];
            for(int i=k+1;i<j;i++){Int v=c.coefficients[i]*limits[i];lo+=min(Int(0),v);hi+=max(Int(0),v);}
            Int a=c.coefficients[k],lower=c.lower-partial-hi,upper=c.upper-partial-lo;
            if(a>0){low=max(low,ceil_div(lower,a));high=min(high,floor_div(upper,a));}
            else if(a<0){high=min(high,floor_div(lower,a));low=max(low,ceil_div(upper,a));}
            else if(lower>0||upper<0)return;
            if(low>high)return;
        }
        for(int x=int(low);x<=int(high);x++)if((Int(x)*d[k])%d[j]==0){p[j][k]=x;fixed_row(j,k+1,limits,ml,adj,det_a,cs,d);}
        p[j][k]=0;return;
    }
    nodes[j]++;
    // The final linear constraint in cs already proves the plus minor.
    bool ok=true;
    for(int i=0;i<j;i++){
        Int value=cs[i].constant;for(int k=0;k<j;k++)value+=cs[i].coefficients[k]*p[j][k];
        Int denom=det_a*d[i];
        if(value%denom||value<0||value>ml[i]*denom){ok=false;break;}
        m[j][i]=value/denom;
    }
    Int leading=fixed_minus_constant[j];for(int i=0;i<j;i++)leading-=fixed_minus_coefficients[j][i]*m[j][i];
    if(ok&&leading>0&&common_row_valid(j)&&valid_diagonal(d,j+1))rows(j+1);
    for(int i=0;i<j;i++)m[j][i]=0;
}
void row_with_diagonal(int j,const int* limits,const int* ml,const Vec& d){
    if(!valid_diagonal(d,j))return;
    auto inherited_before=inherited_equations[j];
    for(int i=1;i<n;i++){Vec row{};row[i]=d[0];row[0]=-d[i];add_equation(inherited_equations[j],row);}
    Mat a=A(p),am=A(m),adj{};Int denom=det(a,j);
    fixed_minus_constant[j]=(2-m[j][j])*det(am,j);
    fixed_minus_coefficients[j]={};
    for(int i=0;i<j;i++)for(int t=0;t<j;t++){
        Mat minor{};for(int r=0,rr=0;r<j;r++)if(r!=t){for(int c=0,cc=0;c<j;c++)if(c!=i)minor[rr][cc++]=am[r][c];rr++;}
        fixed_minus_coefficients[j][i]+=((i+t)%2?-1:1)*det(minor,j-1)*m[t][j];
    }
    for(int i=0;i<j;i++)for(int t=0;t<j;t++){
        Mat minor{};for(int r=0,rr=0;r<j;r++)if(r!=t){for(int c=0,cc=0;c<j;c++)if(c!=i)minor[rr][cc++]=a[r][c];rr++;}
        adj[i][t]=((i+t)%2?-1:1)*det(minor,j-1);
    }
    Vec rhs{};for(int i=0;i<j;i++)for(int t=j;t<n;t++)rhs[i]+=(am[i][t]*a[j][t]-a[i][t]*am[j][t])*d[t];
    vector<Linear> cs;
    for(int i=0;i<j;i++){
        Linear c;c.upper=ml[i]*denom*d[i];
        for(int t=0;t<j;t++){
            c.constant-=adj[i][t]*rhs[t];
            for(int k=0;k<j;k++)c.coefficients[k]+=adj[i][t]*am[t][k]*d[k];
        }
        cs.push_back(c);
    }
    Linear leading;leading.constant=(2-p[j][j])*denom;leading.lower=1;
    for(int k=0;k<j;k++)for(int t=0;t<j;t++)leading.coefficients[k]+=adj[k][t]*a[t][j];
    leading.upper=leading.constant;
    for(int k=0;k<j;k++)leading.upper+=max(Int(0),leading.coefficients[k]*limits[k]);
    cs.push_back(leading);
    fixed_row(j,0,limits,ml,adj,denom,cs,d);
    inherited_equations[j]=move(inherited_before);
}
vector<int> ratio_choices(int i,int k){
    long long cycle_gcd=0;
    function<void(int,int,long long)> visit=[&](int at,int seen,long long product){
        if(at==i){cycle_gcd=gcd(cycle_gcd,product*weight(i,k));return;}
        for(int next=0;next<n;next++)if(!(seen>>next&1)&&weight(at,next))visit(next,seen|(1<<next),product*weight(at,next));
    };
    visit(k,1<<k,1);vector<int> values;
    for(int q=1;q<=bound;q++)if(!cycle_gcd||cycle_gcd%q==0)values.push_back(q);
    return values;
}
bool partial_diagonal(int completed,const Equations& eq){
    // Restrict integrality as soon as two coordinates are proportional on
    // the current kernel, before choosing ratios in other components.
    Vec functions[MAX_N]{};Int denominators[MAX_N]{};
    for(int i=0;i<n;i++){functions[i][i]=1;denominators[i]=1;}
    for(const auto& row:eq){int k=pivot(row);functions[k]={};denominators[k]=row[k];for(int i=0;i<n;i++)if(i!=k)functions[k][i]=-row[i];}
    Int numerator[MAX_N][MAX_N]{},denominator[MAX_N][MAX_N]{};
    for(int i=0;i<n;i++)for(int k=0;k<n;k++){
        int first=0;while(first<n&&!functions[i][first])first++;
        if(first==n)return false;
        bool proportional=true;
        for(int t=0;t<n;t++)if(functions[k][t]*functions[i][first]!=functions[i][t]*functions[k][first]){proportional=false;break;}
        if(!proportional)continue;
        Int a=functions[k][first]*denominators[i],b=functions[i][first]*denominators[k];
        if(a*b<=0)return false;if(b<0){a=-a;b=-b;}
        Int g=gcd128(a,b);a/=g;b/=g;numerator[i][k]=a;denominator[i][k]=b;
        if(i<completed||k>=i)for(const auto& mat:{p,m}){
            Int c=mat[i][k]*a;if(c%b||c/b>bound)return false;
        }
    }
    for(int i=0;i<completed;i++){
        Int total=0;bool fixed=true;
        for(int k=0;k<n;k++){
            Int c=(2*(i==k)-p[i][k])*(2*(i==k)-m[i][k]);if(!c)continue;
            if(!denominator[i][k]){fixed=false;break;}
            Int a=c*numerator[i][k],b=denominator[i][k];if(a%b)return false;total+=a/b;
        }
        if(fixed&&total%2)return false;
    }
    return true;
}
void branch_diagonal(int j,const int* limits,const int* ml,const Equations& eq){
    if(!possible_positive(eq)||!partial_diagonal(j,eq))return;
    if(int(eq.size())==n-1){row_with_diagonal(j,limits,ml,unique_diagonal(eq));return;}
    int from=-1,to=-1;vector<int> choices;
    for(int i=0;i<n;i++)for(int k=0;k<n;k++)if(i!=k&&(i<j||k>=i)&&weight(i,k)){
        Vec ei{},ek{};ei[i]=1;ek[k]=1;auto probe=eq;add_equation(probe,ei);int rank=probe.size();add_equation(probe,ek);
        if(int(probe.size())==rank)continue;
        auto values=ratio_choices(i,k);
        if(from<0||values.size()<choices.size()){from=i;to=k;choices=move(values);}
    }
    if(from<0){free_row(j,0,limits,ml);return;}
    for(int q:choices){auto branch=eq;Vec row{};row[to]=weight(from,to);row[from]=-q;add_equation(branch,row);branch_diagonal(j,limits,ml,branch);}
}
bool complete_principal_row(int j,const int* limits,const int* ml){
    int size=j+1;if(!hereditary_pruning||size>=n)return false;
    for(int k=size;k<n;k++){
        bool a=false,b=false;for(int i=0;i<size;i++){a|=p[i][k]!=0;b|=m[i][k]!=0;}
        if(a&&b)return false;
    }
    string key;
    for(const auto& b:{p,m}){
        for(int i=0;i<j;i++)for(int k=0;k<size;k++)key.push_back(char(b[i][k]));
        key.push_back(char(b[j][j]));
    }
    auto found=completion_tables[size].find(key);if(found==completion_tables[size].end())return true;
    for(const auto& row:found->second){
        bool okay=true;for(int k=0;k<j;k++)if(row[k]>limits[k]||row[j+k]>ml[k]){okay=false;break;}
        if(!okay)continue;
        for(int k=0;k<j;k++){p[j][k]=row[k];m[j][k]=row[j+k];}
        nodes[j]++;
        if(principal(p,j)&&principal(m,j))rows(j+1);
    }
    for(int k=0;k<j;k++)p[j][k]=m[j][k]=0;
    return true;
}
bool hereditary_completed_subsets(int size,Equations& equations){
    if(!hereditary_pruning)return true;
    int latest=1<<(size-1),previous=latest-1;
    for(int part=previous;;part=(part-1)&previous){
        int subset=part|latest;
        if(subset!=(1<<n)-1){
            bool applicable=true;
            for(int k=0;k<n;k++)if(!(subset>>k&1)){
                bool a=false,b=false;
                for(int i=0;i<size;i++)if(subset>>i&1){a|=p[i][k]!=0;b|=m[i][k]!=0;}
                if(a&&b){applicable=false;break;}
            }
            if(applicable){
                hereditary_checks++;int remaining=subset;
                while(remaining){
                    int component=remaining&-remaining;
                    for(;;){int more=component;
                        for(int i=0;i<size;i++)if(component>>i&1)for(int k=0;k<size;k++)if((remaining>>k&1)&&(weight(i,k)||weight(k,i)))more|=1<<k;
                        if(more==component)break;component=more;
                    }
                    remaining&=~component;string key;
                    for(const auto& b:{p,m})for(int i=0;i<size;i++)if(component>>i&1)for(int k=0;k<size;k++)if(component>>k&1)key.push_back(char(b[i][k]));
                    if(!principal_constants[__builtin_popcount(unsigned(component))].count(key)){hereditary_rejections++;return false;}
                    if(weighted_principal_pruning){
                        vector<int> vertices;for(int i=0;i<size;i++)if(component>>i&1)vertices.push_back(i);
                        const Vec& d=principal_diagonals[vertices.size()].at(key);
                        for(int i=1;i<int(vertices.size());i++){Vec row{};row[vertices[i]]=d[0];row[vertices[0]]=-d[i];add_equation(equations,row);}
                        if(!possible_positive(equations)){hereditary_rejections++;return false;}
                    }
                }
            }
        }
        if(!part)break;
    }
    return true;
}
void rows_after_partition(int j,Equations eq);
void rows(int j){
    auto eq=inherited_equations[j-1];
    // Only the newest row adds equations. All earlier equations, including
    // principal-block weight ratios and chosen branches, are inherited.
    for(int i=0;i<j-1;i++){
        Vec row{};bool positive=false,negative=false;
        for(int k=0;k<n;k++){
            row[k]=(2*(i==k)-p[i][k])*(2*(j-1==k)-m[j-1][k])-(2*(i==k)-m[i][k])*(2*(j-1==k)-p[j-1][k]);
            positive|=row[k]>0;negative|=row[k]<0;
        }
        if(positive!=negative)return;
        add_equation(eq,row);if(!possible_positive(eq))return;
    }
    int i=j-1;
    {
        vector<int> targets,masses;unsigned long long choices=1;bool odd=false;
        for(int k=0;k<n;k++)if(p[i][k]&&m[i][k]){
            int mass=int(p[i][k]*m[i][k]);targets.push_back(k);masses.push_back(mass);odd|=mass%2;
            choices=min(1000000ULL,choices*(mass+1));
        }
        if(targets.size()==1&&odd)return;
        if(odd&&choices<=512){
            // Cross-products have no zero exponent by disjoint support.
            // Every leading/diagonal-N product has one positive and one
            // negative exponent with the same coefficient, since 0<p<r_i.
            // Thus the cross-products have equal positive and negative
            // weighted masses. Use all possible signed counts;
            // this is a necessary relaxation, without guessing any exponent.
            Vec relation{};
            function<void(int)> branch=[&](int at){
                if(at==int(targets.size())){
                    int first=0;while(first<n&&!relation[first])first++;
                    if(first==n){rows_after_partition(j,eq);return;}
                    if(relation[first]<0)return;
                    auto next=eq;add_equation(next,relation);
                    if(possible_positive(next))rows_after_partition(j,move(next));
                    return;
                }
                int k=targets[at],mass=masses[at];
                for(int value=-mass;value<=mass;value+=2){relation[k]=value;branch(at+1);}
                relation[k]=0;
            };
            branch(0);return;
        }
    }
    rows_after_partition(j,move(eq));
}
void rows_after_partition(int j,Equations eq){
    if(!hereditary_completed_subsets(j,eq)||!possible_positive(eq))return;
    inherited_equations[j]=eq;
    if(j==n){
        if(!strong())return;
        // All independent column choices have positive leading minors by
        // the mixed Schur budgets. The mixed-column M-matrix lemma gives
        // one common strict positive left vector, without an extreme-ray LP.
        labelled++;emit_diagonals(move(eq));return;
    }
    // Completed rows must already have a path to an unfinished row.
    // Otherwise they contain a closed component that no new row can repair.
    int reachable=((1<<n)-1)^((1<<j)-1);
    for(;;){int more=reachable;for(int i=0;i<j;i++)for(int k=0;k<n;k++)if((reachable>>k&1)&&weight(i,k))more|=1<<i;if(more==reachable)break;reachable=more;}
    if(reachable!=(1<<n)-1)return;
    int limits[MAX_N]{},ml[MAX_N]{};row_bounds(j,limits,ml);
    for(int k=0;k<j;k++)if(limits[k]<0||ml[k]<0)return;
    if(complete_principal_row(j,limits,ml))return;
    if(!prepare_common_budgets(j))return;
    if(int(eq.size())==n-1){row_with_diagonal(j,limits,ml,unique_diagonal(eq));return;}
    unsigned long long diagonal_estimate=1,row_estimate=1;
    for(int k=0;k<n-1-int(eq.size());k++)diagonal_estimate*=bound;
    for(int k=0;k<j;k++)row_estimate*=ml[k]+1;
    int connected=1;for(;;){int more=connected;for(int i=0;i<n;i++)if(connected>>i&1)for(int k=0;k<n;k++)if(weight(i,k)||weight(k,i))more|=1<<k;if(more==connected)break;connected=more;}
    if(connected==(1<<n)-1&&diagonal_estimate<=row_estimate){branch_diagonal(j,limits,ml,eq);return;}
    free_row(j,0,limits,ml);
}
