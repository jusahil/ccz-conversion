"""Strict signed Pauli record parser and two-method certificate reconciliation."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import csv, re

def require(ok, message):
    if not ok: raise ValueError(message)

def parse(text):
    branches={}; certs={}; hist={}; total=None; maximum=None
    for line in text.splitlines():
        parts=line.split()
        if not parts: raise ValueError('unexpected empty record')
        if parts[0] in ('BR','CERT'):
            require(len(parts)==6,'invalid branch field count');tag=parts[0];a,b,e=map(int,parts[1:4]);key=(a,b,e)
            require(0<=a<128 and 0<=b<128 and (a or b) and e in (-1,1),'invalid Pauli key')
            d=branches if tag=='BR' else certs;require(key not in d,'duplicate branch key')
            v=tuple(map(int,parts[4:])) if tag=='BR' else tuple(map(F,parts[4:]));d[key]=v
            continue
        m=re.fullmatch(r'(nullity|P2|moment|spectrum bound|bound) (\S+) (?:count )?(\d+)',line)
        if m:
            tag,v,k=m.groups();key=({'P2':'moment','spectrum bound':'bound'}.get(tag,tag),F(v));require(key not in hist,'duplicate histogram row');hist[key]=int(k);require(int(k)>0,'invalid histogram mass');continue
        m=re.fullmatch(r'branches (\d+)',line)
        if m: require(total is None,'duplicate total');total=int(m[1]);continue
        m=re.fullmatch(r'max spectrum bound (\S+)(?: witness \(a,b,sign\) \d+,\d+,-?1)?',line)
        if m: require(maximum is None,'duplicate max');maximum=F(m[1]);continue
        raise ValueError('unexpected record: '+line)
    require(set(hist) and all(k[0] in {'nullity','moment','bound'} for k in hist),'empty histogram')
    return branches,certs,hist,total,maximum

def reconcile(results, data):
    a,b,c,d,e=parse((results/'C08.out').read_text()); aa,bb,cc,dd,ee=parse((results/'C09.out').read_text())
    x,y,h,_,_=parse((results/'S7-F.out').read_text())
    require((a,b,c,d,e)==(aa,bb,cc,dd,ee),'independent methods disagree');require(not x and not y and h==c,'third regression disagrees')
    keys={(i,j,s) for i in range(128) for j in range(128) if i or j for s in (-1,1)}
    require(set(a)==keys and d==len(keys),'incomplete key domain');require(set(b)=={k for k,v in a.items() if v[0]==6},'certificate domain mismatch')
    for (i,j,s),(n,den) in a.items(): require(0<=n<=6 and 0<den<256 and den+a[i,j,-s][1]==256,'invalid denominator/nullity')
    require(all(0<m<=1 and 0<=cap<1 for m,cap in b.values()),'invalid certificate entry')
    observed=Counter(('nullity',F(n)) for n,_ in a.values());observed.update(('moment',m) for m,_ in b.values());observed.update(('bound',cap) for _,cap in b.values())
    require(observed==c,'aggregate disagrees with signed keys');require(d==32766 and len(b)==31644 and e==max(v[1] for v in b.values())==F(121,128),'certificate contract failed')
    with (data/'fano_histograms.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['quantity','value','count']);w.writerows(sorted((t,str(v),n) for (t,v),n in c.items()))
    with (data/'signed_pauli_branches.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['a','b','sign','nullity','denominator','P2_if_nullity6','upper_bound_if_nullity6'])
        for k,v in sorted(a.items()):w.writerow([*k,*v,*b.get(k,('',''))])
    return {'branches':len(a),'relevant':len(b),'maximum_bound':str(e),'maximizing_keys':sum(v[1]==e for v in b.values())}
