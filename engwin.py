#!/usr/bin/env python3
"""engwin.py — VASP EIGENVAL 能量提取"""
import sys, os

def read_eigenval(fn):
    with open(fn) as f: lines = f.readlines()
    hdr = None
    for i, l in enumerate(lines):
        p = l.split()
        if len(p) in (2,3):
            try:
                nk, nb = int(p[0]), int(p[1])
                if 1<=nk<=500000 and 1<=nb<=50000: hdr = i; nkpts = nk; break
            except: pass
    if hdr is None: print("ERROR"); sys.exit(1)
    print(f"  nkpts={nkpts}")
    data = lines[hdr+1:]
    ik, ef, i = 0, [], 0
    while ik < nkpts and i < len(data):
        p = data[i].split(); i += 1
        if not p: continue
        if len(p) == 4:
            ik += 1; k = []
            while i < len(data):
                bp = data[i].split(); i += 1
                if not bp: continue
                if len(bp) == 4: i -= 1; break
                elif len(bp) >= 2: k.append(float(bp[1]))
            ef.append(k)
    nb = len(ef[0]) if ef else 0
    print(f"  nbands={nb} ({nkpts} k点)")
    return ef, nb

def main():
    if len(sys.argv) < 3: print("Usage: engwin.py EIGENVAL [e N | n lo hi | s]"); sys.exit(1)
    fn, op = sys.argv[1], sys.argv[2]
    if not os.path.exists(fn): print(f"ERROR: {fn}"); sys.exit(1)
    ef, nb = read_eigenval(fn)
    if op == "e":
        b = int(sys.argv[3])
        if b<1 or b>nb: print(f"band {b} 超出 [1,{nb}]"); sys.exit(1)
        e = [ek[b-1] for ek in ef if len(ek)>=b]
        print(f"band {b}: emin={min(e):.4f} emax={max(e):.4f} eV")
    elif op == "n":
        lo,hi = float(sys.argv[3]), float(sys.argv[4])
        for ik, ek in enumerate(ef):
            nc = sum(1 for e in ek if lo<=e<=hi)
            print(f"ik = {ik+1}, nbnd = {nc}")
    elif op == "s":
        print(f"{'#':>4}|{'emin':>12}|{'emax':>12}|{'ΔE':>8}")
        for b in range(min(60,nb)):
            e=[ek[b] for ek in ef if len(ek)>b]
            if e: print(f"{b+1:4d}|{min(e):12.4f}|{max(e):12.4f}|{max(e)-min(e):8.4f}")
if __name__=="__main__": main()
