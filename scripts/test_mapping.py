"""Benchmark PORT-EK-v2 mapping on simulated edits with known true position."""
import os, sys, random, itertools; sys.path.insert(0, os.environ.get('PORTEK_DIR') or os.getcwd())
from portek.portek_map import MappingPipeline
random.seed(7); K=15; D=3
ref="".join(random.choice("ACGT") for _ in range(3000))
mp=MappingPipeline.__new__(MappingPipeline); mp.k=K
mp.kmer_index=mp._generate_index(D,ref)
def mutate(kmer,pos,kind,n):
    s=list(kmer)
    if kind=="sub":
        for p in random.sample(range(K),n): s[p]=random.choice([c for c in "ACGT" if c!=s[p]])
    elif kind=="ins":   # query has n extra bases (insertion relative to reference)
        p=random.randint(2,K-n-3); ins="".join(random.choice("ACGT") for _ in range(n)); s=s[:p]+list(ins)+s[p:]; s=s[:K]
    elif kind=="del":   # query lacks n bases; k-mer is completed with the downstream reference bases
        p=random.randint(2,K-n-3); s=s[:p]+s[p+n:]; s=s+list(ref[pos+K:pos+K+n])
    return "".join(s)
rows=[]
N=300
for kind in ("sub","ins","del"):
    for n in (0,1,2,3):
        if kind!="sub" and n==0: continue
        ok=wrong=unm=multi=0
        for _ in range(N):
            pos=random.randint(10,len(ref)-K-10); q=mutate(ref[pos:pos+K],pos,kind,n) if n else ref[pos:pos+K]
            md=mp.map_kmer_to_index(q,D)
            hit=None
            for d in range(D+1):
                if md[d]: hit=(d,sorted(md[d])); break
            if hit is None: unm+=1
            else:
                if len(hit[1])>1: multi+=1
                if pos+1 in hit[1]: ok+=1
                else: wrong+=1
        rows.append((kind,n,ok/N,unm/N,wrong/N,multi/N))
print("edit  n  correct_pos  unmapped  mapped_wrong_pos  multi_hit")
for r in rows: print(f"{r[0]:4s} {r[1]:2d}   {r[2]:.2f}        {r[3]:.2f}      {r[4]:.2f}             {r[5]:.2f}")
# combined sub+ins
ok=0
for _ in range(N):
    pos=random.randint(10,len(ref)-K-10); q=mutate(mutate(ref[pos:pos+K],pos,"ins",1),pos,"sub",1)
    md=mp.map_kmer_to_index(q,D); ok+= any(pos+1 in md[d] for d in range(D+1))
print("1 insertion + 1 substitution (d<=3): correct", ok/N)
