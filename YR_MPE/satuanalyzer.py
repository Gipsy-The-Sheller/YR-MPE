from math import *

def CalcEn(ra, rc, rg, rt):
    En=0
    if ra!=0:
        En+=ra*log2(ra)
    if rc!=0:
        En+=rc*log2(rc)
    if rg!=0:
        En+=rg*log2(rg)
    if rt!=0:
        En+=rt*log2(rt)
    return -En

def CalcHm(alignment):
    H, N=0, 0
    for i in range(len(alignment[0])):
        na, ng, nc, nt, n=0, 0, 0, 0, 0
        for j in range(len(alignment)):
            if alignment[j][i] == 'A' or alignment[j][i] == 'a':
                na+=1
                n+=1
            elif alignment[j][i] == 'G' or alignment[j][i] =='g':
                ng+=1
                n+=1
            elif alignment[j][i] == 'T' or alignment[j][i] == 't':
                nt+=1
                n+=1
            elif alignment[j][i] == 'C' or alignment[j][i] == 'c':
                nc+=1
                n+=1
        if n!=0:
            H+=CalcEn(na/n, nc/n, ng/n, nt/n)
            N+=1
    return H/N

def CalcHfss(alignment):
    # stat pa, pc, pg, pt
    sa, sc, sg, st = 0, 0, 0, 0
    for i in alignment:
        i = i.upper()
        sa += i.count('A')
        sc += i.count('C')
        sg += i.count('G')
        st += i.count('T')
    pa = sa / (len(alignment)*len(alignment[0]))
    pc = sc / (len(alignment)*len(alignment[0]))
    pg = sg / (len(alignment)*len(alignment[0]))
    pt = st / (len(alignment)*len(alignment[0]))
    # calculate ideal Fss for random sequence set with length same of the alignment
    length = len(alignment)
    Hfss = 0
    # randomize base frequency
    for a in range(length+1):
        for c in range(length+1-a):
            for g in range(length+1-a-c):
                t = length-a-c-g
                # multinominal distribution (dimension=4)
                mdis = factorial(length)/(factorial(a)*factorial(c)*factorial(g)*factorial(t))*(pa**a)*(pc**c)*(pg**g)*(pt**t)
                entropy = CalcEn(a/length, c/length, g/length, t/length)
                Hfss += mdis*entropy
    return Hfss

def CalcHfss_jc69(alignment):
    # stat pa, pc, pg, pt
    pa, pc, pg, pt=0.25, 0.25, 0.25, 0.25
    # calculate ideal Fss for random sequence set with length same of the alignment
    length = len(alignment)
    Hfss = 0
    # randomize base frequency
    for a in range(length+1):
        for c in range(length+1-a):
            for g in range(length+1-a-c):
                t = length-a-c-g
                # multinominal distribution (dimension=4)
                mdis = factorial(length)/(factorial(a)*factorial(c)*factorial(g)*factorial(t))*(pa**a)*(pc**c)*(pg**g)*(pt**t)
                entropy = CalcEn(a/length, c/length, g/length, t/length)
                Hfss += mdis*entropy
    return Hfss