"""Derive exact data from the executed programs and analytic laws."""
import ast,csv,re
from fractions import Fraction as F
from itertools import permutations
from verification.branch_records import require
from src.spectra import phase_spectrum, spec_prod, overlap

def write(data,name,header,rows):
    rows=list(rows);require(rows and all(len(r)==len(header) for r in rows),name+' empty/malformed')
    with (data/(name+'.csv')).open('w',newline='') as f: w=csv.writer(f);w.writerow(header);w.writerows(rows)
def matches(text,pattern,count):
    rows=re.findall(pattern,text,re.M);require(len(rows)==count,'pattern count: '+pattern);return rows

def export(results,data,branch_max):
    text=(results/'C06.out').read_text();rows=[]
    for n,dim,h in matches(text,r'^fully charged N (\d+) dimension (\d+) rank counts (\{[^}]+\})$',4):
        n,dim=int(n),int(dim);d=ast.literal_eval(h);require(dim==3*(n-1) and sum(d.values())==2**dim and all(type(k)==type(v)==int and k%2==0 and 0<=k<=dim and v>0 for k,v in d.items()),'copy ranks');rows.extend((n,dim,k,v) for k,v in sorted(d.items()))
    require({r[0] for r in rows}==set(range(3,7)),'copy domain');write(data,'copy_ranks',['N','dimension','contraction_rank','count'],rows)
    rows=[]
    for h,n in matches(text,r'^contraction ranks (\(.*\)) references (\d+)$',4):
        pairs=ast.literal_eval(h);a=dict(pairs);require(len(a)==len(pairs) and a[0]>0 and a[0]&(a[0]-1)==0 and sum(a.values())==64,'invalid determinant rank histogram');rows.append((6-(a[0].bit_length()-1),int(n),h))
    require(sum(x[1] for x in rows)==1395 and len({x[0] for x in rows})==4,'plane domain');write(data,'determinant_planes',['nullity','planes','rank_histogram'],rows)
    rows=matches(text,r'^N (\d+) all branch identities checked (\d+)$',4);require({int(n) for n,_ in rows}==set(range(3,7)) and all(int(k)==8*2**(3*(int(n)-2)) for n,k in rows),'decoder counts');write(data,'decoder_checks',['N','assignments'],rows)
    rows=matches(text,r'^reference N (\d+) ell (\d+) eta (\d+)$',9);require({(int(n),int(e)) for n,e,_ in rows}=={(n,e) for n in range(2,5) for e in range(1,4)},'reference domain');write(data,'reference_nullities',['N','ell','nullity'],rows)
    rows=matches((results/'C07.out').read_text(),r'^labelled vertices (\d+) graphs (\d+) cycles in complete graph (\d+) active graphs (\d+)$',6);require({int(r[0]) for r in rows}==set(range(1,7)) and sum(int(r[1]) for r in rows)==33867,'graph domain');write(data,'graph_counts',['vertices','graphs','complete_graph_cycles','active'],rows)
    rows=matches((results/'C04.out').read_text(),r'^(.+) V,E,F,T,b1= (\d+) (\d+) (\d+) (\d+) (\d+) cup checks (\d+) elementary exact shifts (\d+)$',4);require(len({r[0] for r in rows})==4 and sum(int(r[-2]) for r in rows)==93 and sum(int(r[-1]) for r in rows)==2565,'mesh domain');write(data,'mesh_counts',['mesh','V','E','F','T','b1','cup_checks','exact_shift_checks'],rows)
    rows=matches((results/'C01.out').read_text(),r'^monomial (\d+) (\d+) eta (\d+) P2 (\S+) P3 (\S+)$',12);require({(int(d),int(n)) for d,n,*_ in rows}=={(d,n) for d in range(4) for n in range(1,4)} and all(0<F(r[-1])<=1 and 0<F(r[-2])<=1 for r in rows),'moment domain');write(data,'monomial_moments',['relative_type','N','nullity','P2','P3'],rows)
    write(data,'batch_yield',['N','ell_representative','yield'],((n,e,max(0,n-2) if e==0 else n-1 if e==1 else n) for e in (0,1,3) for n in range(1,13)))
    write(data,'cycle_noise',['L','u_exact','F_exact'],((L,u,((1+(1-2*u)**L)/2)**3) for L in (3,5,9,21) for u in [F(i,10000) for i in range(901)]))
    write(data,'reference_price',['b','one_source','two_sources','three_sources','minimum'],((b,1+3*b,2+b,3,min(1+3*b,2+b,3)) for b in [F(i,200) for i in range(301)]))
    C=phase_spectrum(3,lambda x:int(x==7));T=spec_prod(C,C);T3=spec_prod(T,C)
    trace=phase_spectrum(6,lambda x:int((x&7)==7)^int((x>>3)==7)^int(((x&7)^(x>>3))==7))
    det=phase_spectrum(9,lambda x:sum(all(x&(1<<(3*i+p[i])) for i in range(3)) for p in permutations(range(3)))%2)
    write(data,'fidelity_ceilings',['task','value','status'],[('A: 2 targets',F(1),'attained'),('B: 2 targets',F(branch_max),'upper bound'),('F4: 2 targets',overlap(trace,T,6),'upper bound'),('Det: 3 targets',overlap(det,T3,9),'upper bound')])
