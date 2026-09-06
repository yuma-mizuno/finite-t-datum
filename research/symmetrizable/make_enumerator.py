from pathlib import Path
base=Path('research/higher_rank/enumerate_constants.cpp').read_text()
start=base.index('void rows(int j);')
end=base.index('void triangular(')
base=base[:start]+Path('research/symmetrizable/constant_rows.inc').read_text()+base[end:]
start=base.index('bool principal(const Mat& b,int j){');end=base.index('int weight(',start)
base=base[:start]+'''// Earlier leading blocks are nonsingular M-matrices. For a Z-matrix,
// positivity of the next Schur complement is equivalent to the M-matrix
// property of the enlarged block, hence to all its principal minors.
bool principal(const Mat& b,int j){return det(A(b),j+1)>0;}
'''+base[end:]
base=base.replace('// Exhaustive constants for ranks 2--6, with parity pruning and resumable tasks.','// Exhaustive symmetrizable constants, with an unfixed positive diagonal.')
base=base.replace('p=triangle[a];m=triangle[b];if(!parity(0))continue;','p=triangle[a];m=triangle[b];')
base=base.replace('\\\"parity_pruning\\\":true','\\\"weighted_parity_pruning\\\":true')
base=base.replace('hereditary_pruning=true;','hereditary_pruning=true;build_completion_tables();load_principal_diagonals(string(argv[6])+".weights");')
base=base.replace('#include <unordered_set>','#include <unordered_map>\n#include <unordered_set>')
base=base.replace('int shard=stoi(argv[3]),shards=stoi(argv[4]),only=argc>5?stoi(argv[5]):-1;',
                  'int shard=stoi(argv[3]),shards=stoi(argv[4]),only=argc>5?stoi(argv[5]):-1;\n    int first_partner=argc>7?stoi(argv[7]):-1,last_partner=argc>8?stoi(argv[8]):-1;')
base=base.replace('"upper-"+to_string(a)+".json"',
                  '"upper-"+to_string(a)+(first_partner>=0?"-part-"+to_string(first_partner)+"-"+to_string(last_partner):"")+".json"')
base=base.replace('for(int b=a;b<int(triangle.size());b++){',
                  'for(int b=max(a,first_partner);b<(last_partner>=0?min(last_partner,int(triangle.size())):int(triangle.size()));b++){')
Path('research/symmetrizable/enumerate_constants.cpp').write_text(base,encoding='utf8',newline='\n')
