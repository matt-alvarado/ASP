"""Parser for MCM v3.3.1 KPP files (species, equations, constants)."""
import re, collections
D='/home/claude/mcm_src/mcm331/mcm_3-3-1_windows/'
def read_eqn(path=D+'mcm_3-3-1_kpp_complete.eqn'):
    txt=open(path,errors='ignore').read()
    body=txt.split('#EQUATIONS')[1]
    rxns=[]
    for m in re.finditer(r'<(\d+)>\s*([^:]*?)=\s*([^:]*?):\s*(.*?);',body,re.S):
        n,l,r,k=m.groups()
        def sp(s):
            out=[]
            for t in s.split('+'):
                t=t.strip()
                if not t: continue
                mm=re.match(r'^([0-9.]+)\s*([A-Za-z]\S*)$',t)
                out.append((float(mm.group(1)),mm.group(2)) if mm else (1.0,t))
            return out
        rxns.append(dict(n=int(n),re=sp(l),pr=sp(r),k=' '.join(k.split())))
    return rxns
if __name__=='__main__':
    R=read_eqn(); print(len(R))
    tmpl=collections.Counter()
    for r in R:
        t=re.sub(r'(?<![A-Za-z_])[0-9]+\.?[0-9]*(?:[ED][+-]?[0-9]+)?','#',r['k'])
        tmpl[t]+=1
    print(len(tmpl))
    for t,c in tmpl.most_common(70): print(c,t)
    print(max(len(r['pr']) for r in R), collections.Counter(len(r['re']) for r in R))
    names=set(n for r in R for _,n in r['re']+r['pr']); print(len(names), max(map(len,names)), sum(len(n)>15 for n in names))
