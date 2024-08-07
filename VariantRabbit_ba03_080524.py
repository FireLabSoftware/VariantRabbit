#!/usr/bin/env python
## - VariantRabbit: A fast Reference-based means to search for consistent variants in NGS Datasets
##
## - Version ba03 08-05-24
##
## - OVERVIEW OF ** VariantRabbit **
## - You have a datasets from next generation sequencing of different samples, and a reference, and want to find sample-specific variants
## - VariantRabbit is a tool providing lists of short assembled segments that are indicative of variants that meet a settable threshold
## - So your input will be a set of individual data files (fasta or fastq, can be gzipped) and (optionally a reference genome that will be used to
## -  1. Locate any perfect or imperfect match sequences amongst the assembled contigs.
## -  2. Optionally remove any perfect reference matches from the outputted list of FastAs
## - The output will be a fastA of short assemblies (mutant and reference) that can then be analyzed for their presence in different datasets

## - SYNTAX Example:
## - pypy3  DataPos='./PositiveDataFolder/*.fastq.gz'  DataPos='./NegativeDataFolder/*.fastq.gz' RefFile=myRefSeq.fa

## - Input data
DataPos1 = 'default' ## example: './Filtered/ObP_*.fastq.gz' ## Sequence Data File Name for sample "Positive Condition"
DataNeg1 = 'default'  ## example:'./Filtered/ObN_*.fastq.gz' ## Sequence Data File Name for sample "Negative Condition"

## - Optional Reference File
RefFile1 = ''  ## example:'SK36.fasta'    ## Reference File for comparisons.  

## - Parameters for k-mer list assembly
klen1 = 31               ## k-mer length for comparison and assembly
minAssemblyLength1 = 'default' ## setting this to a positive value requires that any reported assembly be at least this long
RequireFullFlank1 = True ## Relevant only if there is a reference-- the setting imposes a requirement for a flanking k-1 mer in each segment that matches the reference
                         ## Setting this to False will try to capture variants where full flanking sequences are not assembled
                         ## Setting True provides a much more stringent test for variants here, while False will provide variants in
                         ## some poorly covered or error prone regions at the expense of specificity.
              

## Parameters that control the difference in k-mer count used to identify differential incidence
upMin1 = 10              ## Minimum k-mer count in a positive sample required to assign a k-mer as possibly enriched 
downMax1 = 1             ## Maximum k-mer count in a negative sample required to assign a k-mer as possibly enriched
minFold1 = 10            ## Maximum fold difference to assign a k-mer as possibly enriched (enrichment means (Pos+Regular)>=minFold*(Neg+Regular)             
Regular1 = 1             ## Regularization constant used in calculating fold difference

## Output Control
OutFileName1 = 'default'   ## User can rely on default file name, or provide a file name here for output
ReportNonMutant1 = False   ## Setting this to true reports all differential assemblies including those that match the reference
Delimiter1 = '\r'          ## Line delimiter
ReportCadence1 = 100000    ## How often to report progress during k-mer list builds

from VSG_ModuleFP import *
from collections import Counter
from glob import glob
import gzip, os
vCommand()

if minAssemblyLength1=='default':
    minAssemblyLength1 = 2*klen1-1
def antisense(s):
    return s.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
def myOpen1(fn,mode):
    if fn.lower().endswith('.gz'):
        return gzip.open(fn, mode=mode)
    else:
        return open(fn, mode=mode)
def multiglob1(x):
    if not(type(x))==str:
        x = ','.join(x)
    l = []
    for n in x.split(','):
        n = n.strip("'").strip('"')
        l.extend(sorted(list(glob(n))))
    return l
mask1 = 4**(klen1-1)-1   ## a mask for the bit-based binarization operation (sense orientation)
ksam1 = 4**(klen1-1)     ## a second mask that allows facile calculation of an antisense k-mers binary value
BaseD1 = Counter({'G':0, 'A':1, 'T':2, 'C':3, 'N':0, 'g':0, 'a':1, 't':2, 'c':3, 'n':0, 'U':2, 'u':2})
BaseL1 = [0]*256         ## an array that allows rapid lookup of the numerical representations of G,A,T,C
BaseA1 = [0]*256         ## an array that allows rapid calculation of antisense bases
for b1 in BaseD1:
    BaseL1[ord(b1)] = BaseD1[b1]
    BaseA1[ord(b1)] = (3-BaseD1[b1])*ksam1

aBaseL1 = ['G','A','T','C']

def LotsOfKMers1(F):  ## simple demonstration of the simultaneous (kmer+anti-kmer) binarization algorithm 
    FList = multiglob1(F)
    C1 = Counter()
    for f in FList:
        if f.lower().endswith('fastq') or f.lower().endswith('fastq.gz'):
            DataCadence1 = 4
        else:
            DataCadence1 = 2
        for i,L in enumerate(myOpen1(f,'rt')):
            if i%(ReportCadence1*DataCadence1)==0:
                vLog('Building KMer List',os.path.basename(f),'at read',i//DataCadence1)
            if i%DataCadence1!=1: continue
            StartMatching1 = klen1-1
            L = L.strip()
            v0 = 0   ## v0 will hold the representation of the sense K-mer ending at the current position
            a0 = 0   ## a0 will hold the representation of the reverse complement
            for j1,c1 in enumerate(L):
                if c1 == 'N': StartMatching1 = j1+klen1
                v0 = ((v0&mask1)<<2)+BaseL1[ord(c1)]  ## a list lookup seems like the fastest operation to get this done
                a0 = (a0>>2)+BaseA1[ord(c1)]
                if j1>=StartMatching1:
                    vamin0 = min(a0,v0)  ## vamin0 is the lesser of v0 and a0, and gives a strand-independent identifier for the k-mer
                    C1[vamin0] += 1      ## this function just makes a counter from the input sequence, but many other operations are possible depending on what is done at this point
    return C1
def ValueToSeq1(v,l):
    ''' converts numerical representation back to sequence v is value, l is k-mer length'''
    return ''.join([aBaseL1[(v>>(2*i)) & 3] for i in range(l-1,-1,-1)])

def FindMe1(s,t,Circular=False):
    ''' a circle-capable finder for a target t in a sequence s.  Returns a list of positions (zero based) as tuples [position, orientation [0=sense, 1=anti]'''
    PosList = []
    if Circular:
        s += s[:len(t)-1]
    p = s.find(t)
    while p>=0:
        PosList.append((p,1))
        p = s.find(t,p+1)
    ant = antisense(t)
    p = s.find(ant)
    while p>=0:
        PosList.append((p,-1))
        p = s.find(ant,p+1)
    return PosList

def PlusMe1(s,t,Circular=True):
    ''' takes a target t putatively from a sequence s and returns a best guess of the sequence or its antisense (whichever matches best
        using a single end-longest match test)'''
    if Circular:
        s = s+s[:len(t)-1]
    for i in range(len(t),1,-1):
        if t[:i] in s:
            return t
        if antisense(t[:i]) in s:
            return antisense(t)
        if t[-i:] in s:
            return t
        if antisense(t[-i:]) in s:
            return antisense(t)
    return t

def FindMe2(s,t,mink=klen1//2,RequireFullFlank=RequireFullFlank1,klen=klen1):
    '''provides a very quick evaluation of where a probe corresponds in a sequence
    returns a list of 4-ples: start position, extrapolated end position, orientation, match-perfect [bool]
    Logic: First look for perfect matches, if these are found, return them; then look for
    longest match at one end of query and return all that match that length
    Special case if there are equally long matches at the two ends of query in which case return the combination that most closely matches the
    expected segment in size (smallest indel).  For ties, a list will be returned'''
    
    myResult=[]
    for p,o in  FindMe1(s,t): myResult.append([p,p+len(t),o,True])
    if myResult: return myResult
    maxk = len(t)-1
    if RequireFullFlank:
        mink = klen-1
        maxk = klen-1
    for i in range(maxk,mink-1,-1):
        myfindUp = FindMe1(s,t[:i])
        myfindDown = FindMe1(s,t[-i:])
        if myfindUp and myfindDown:
            minIndel = len(t)+1
            for p1,o1 in myfindUp:
                if o1==-1: p1 -= len(t)-i
                for p2,o2 in myfindDown:
                    if o2==1: p2 -= len(t)-i
                    if o1!=o2: continue
                    if p2+len(t)-p1<mink-1: continue
                    myDiff = abs(p2-p1)
                    if myDiff>minIndel: continue
                    if myDiff<minIndel:
                        minIndel=myDiff
                        myResult = []
                    myResult.append([p1,p2+len(t),o1,False])
            if myResult: return myResult
        if not(RequireFullFlank):
            for p,o in myfindUp:
                if o==-1: p -= len(t)-i
                myResult.append([p,p+len(t),o,False])
            for p,o in myfindDown:
                if o==1: p -= len(t)-i
                myResult.append([p,p+len(t),o,False])
        if myResult: return myResult
    return []

            
def myAssem1(c, klen=klen1):
    ''' a rather slow and deliberate assembler.  Input is a Counter object with keys being kmers to assemble
        and values being the observed count of each, along with a k-mer length; Output is a Counter object where the
        keys are assemblies and the values are average coverage; ambiguities (forks) in assembly are resolved in favor of the most counts
        circles will end up with terminal duplications of length klen-1'''
    myA = Counter()
    kList = [x[0] for x in c.most_common()]
    c = Counter({x:c[x] for x in kList})
    for k in kList:
        if not(k in c): continue
        mya = k
        myCount = c[k]
        del(c[k]); del(c[antisense(k)])
        myDir = 1
        while True:
            bestCount = 0
            lead1 = mya[-klen+1:]
            for b in 'GATC':
                bl = lead1+b
                bc = c[bl]
                if bc>bestCount:
                    bestCount = bc
                    bestBase = b
                    bestk = bl
            if bestCount:
                mya += bestBase
                myCount += bestCount
                del(c[bestk]); del(c[antisense(bestk)])
            elif myDir == 1:
                mya = antisense(mya)
                myDir = -1
            else:
                break
        myA[mya] = myCount/(len(mya)-klen+1)
    return myA
                    
def CounterDiff1(upC,downC,Regular=Regular1,minFold=minFold1,upMin=upMin1,downMax=downMax1):
    diffC = Counter()
    minRatio = minFold*(sum(upC.values())+Regular)/(sum(downC.values())+Regular)
    for k in upC:
        if upC[k]>=upMin and downC[k]<=downMax and (upC[k]+Regular)>=minRatio*(downC[k]+Regular):
            ks = ValueToSeq1(k,klen1)
            diffC[ks] = upC[k]
            diffC[antisense(ks)] = upC[k]
    return diffC
    
def myRef1(mySeq,myP1,myP2,myOrient):
    ''' returns an arbitrary segment of the indicated sequence'''
    if myOrient == 1:
        return mySeq[myP1:myP2]
    else:
        return antisense(mySeq[myP1:myP2])
                   
if RefFile1:
    mySeq1 = ''.join(vFastAToDict(RefFile1).values())
else:
    mySeq1 = ''
C1 = LotsOfKMers1(DataPos1)
C2 = LotsOfKMers1(DataNeg1)
vLog('Starting Difference Evaluation Pos-Neg')
myC1 = CounterDiff1(C1,C2)
vLog('Starting Difference Evaluation Neg-Pos')
myC2 = CounterDiff1(C2,C1)
vLog('Starting Assemblies Pos')
myA1 = myAssem1(myC1)
vLog('Starting Assemblies Neg')
myA2 = myAssem1(myC2)
Abbrev1 = os.path.basename(DataPos1).split('_')[0]
Abbrev2 = os.path.basename(DataNeg1).split('_')[0]
if OutFileName1=='default':
    OutFileName1 = 'DiffAssembly_'+Abbrev1+'_'+Abbrev2+'_'+vnow.replace('D','').replace('T','').replace('_','')+'.fa'
    OutFile1 = myOpen1(OutFileName1, mode='wt')
knum1 = 1
vLog('Starting FastA Output')

for a0,myA0 in zip([Abbrev1,Abbrev2],[myA1,myA2]):
    for i1,a1 in enumerate(myA0):
        if minAssemblyLength1 and len(a1)<minAssemblyLength1: continue
        if not(mySeq1):
            myCount1 = myA0[a1]
            OutFile1.write('>'+a0+'_seg_'+str(knum1)+'_length_'+str(len(a1))+'_cov_'+'%.2f' % myCount1+Delimiter1)
            OutFile1.write(a1+Delimiter1)
            knum1 += 1
        else:
            antia1 = antisense(a1)
            if not(ReportNonMutant1) and ((a1 in mySeq1) or (antia1 in mySeq1)): continue
            myCount1 = myA0[a1]
            a1 = PlusMe1(mySeq1,a1)
            myPos1 = FindMe2(mySeq1,a1)
            if not(myPos1): continue
            PStrings1 = []        
            for myP1 in myPos1:
                if myP1[2]==1:
                    PStrings1.append(str(myP1[0]+1)+'-'+str(myP1[1]))
                else:
                    PStrings1.append(str(myP1[1])+'-'+str(myP1[0]+1))
            if myPos1[0][3]:
                if ReportNonMutant1:
                    OutFile1.write('>'+a0+'_K_'+str(knum1)+'_Wt_'+','.join(PStrings1)+'_length_'+str(len(a1))+'_cov_'+'%.2f' % myCount1+Delimiter1)
                    OutFile1.write(a1+Delimiter1)
                    knum1 += 1
            else:
                OutFile1.write('>'+a0+'_K_'+str(knum1)+'_Mut_'+','.join(PStrings1)+'_length_'+str(len(a1))+'_cov_'+'%.2f' % myCount1+Delimiter1)
                OutFile1.write(a1+Delimiter1)
                knum1 += 1
                refSeqSet1 = set()
                for myP1 in myPos1:
                    Ref1 = myRef1(mySeq1,myP1[0],myP1[1],myP1[2])
                    if not(Ref1 in refSeqSet1) and not(antisense(Ref1) in refSeqSet1):
                        refSeqSet1.add(Ref1)
                        if myP1[2]==1:
                            OutFile1.write('>'+a0+'_K_'+str(knum1)+'_Ref_'+str(myP1[0]+1)+'-'+str(myP1[1])+'_length_'+str(myP1[1]-myP1[0])+Delimiter1)
                        else:
                            OutFile1.write('>'+a0+'_K_'+str(knum1)+'_Ref_'+str(myP1[1])+'-'+str(myP1[0]+1)+'_length_'+str(myP1[1]-myP1[0])+Delimiter1)
                        OutFile1.write(Ref1+Delimiter1)                
                        knum1 += 1                
OutFile1.close()
vLog('Finished FastA Output')
