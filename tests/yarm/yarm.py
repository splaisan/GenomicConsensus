#!/usr/bin/env python

# Yet Another Reference Mutator (yarm), the purpose of which is to help
# validate existing/new variant algorithms.

import sys
import pickle
import random
from os import getenv
from bisect import bisect
from libyarm import Mutator, Mutation, MutatedStrain

def create():
    '''
    Create a new set of mutation strains
        [typ=sid] [fname=strains.pkl] [fixed=0] [window=100] yarm create <fasta> [# of strains]
    '''
    fasta = sys.argv[2]
    n = int(defarg(3,1))
    strains = []
    for i in range(0,n):
        fd = open(fasta, 'r')
        header = fd.readline()[1:-1]
        strains.append(Mutator().mutate("%i_%s"%(i, header),fd))

    fname = defenv('fname', 'strains.pkl')
    pickle.dump(strains, open(fname,'w'), 2)

def apply():
    '''
    Apply strain information to a reference sequence
        [debug=0] yarm apply <strains> <fasta> [index] > <fasta>
    '''
    import cStringIO
    strains = pickle.load(open(sys.argv[2],'r'))
    inst = open(sys.argv[3],'r')
    n = int(defarg(4,0))
    stra = strains[n]
    inst.readline()
    outs = cStringIO.StringIO()
    mtr = Mutator()
    mtr.apply(stra, inst, outs)
    print ">" + stra.id
    seq = outs.getvalue()
    for i in range(0, len(seq), 80):
        print seq[i:i+80]

def stats():
    '''
    View stats on a set of mutation strains
        yarm stats <strains>
    '''
    strains = pickle.load(open(sys.argv[2],'r'))
    for stra in strains:
        sc = len([x for x in stra.mutations if x.typ == 's'])
        ic = len([x for x in stra.mutations if x.typ == 'i'])
        dc = len([x for x in stra.mutations if x.typ == 'd'])
        tc = sc + ic + dc
        print "%s %i mut: %i, s: %i, i: %i, d: %i" \
            % (stra.id, stra.length, tc, sc, ic, dc)

def dump():
    '''
    Dump all information for strain(s)
        yarm dump <strains>
    '''
    print pickle.load(open(sys.argv[2],'r'))

def invert():
    '''
    Same as `dump`, but invert the coordinates to match reads mapped to strain
        [map=0] yarm invert <strains> [index]
    '''
    strains = pickle.load(open(sys.argv[2],'r'))
    n = int(defarg(3,0))
    asmap = bool(getenv('map'))
    if asmap:
        for i in zip(strains[n].mutations, _invert(strains[n]).mutations):
            print i
    else:
        print _invert(strains[n])

def compare():
    '''
    Compares a ground-truth strain to a variations file
        [typ=[s|i|d]] [id=] [verbose=0] yarm compare <strains> <variants>
    '''
    strains = pickle.load(open(sys.argv[2],'r'))
    ms = MutatedStrain.fromGff(open(sys.argv[3],'r'))

    sid = defenv('id', ms.id)
    typ = getenv('typ')
    groundTruth = [stra for stra in strains if stra.id == sid]
    try:
        [groundTruth] = groundTruth
        seqlen = groundTruth.length
        groundTruth = _invert(groundTruth)
        groundTruth.mutations.sort()
    except ValueError:
        raise LookupError('Unable to locate %s' % ms.id)

    if typ:
        act = set([x for x in ms.mutations if x.typ==typ])
        exp = set([x for x in groundTruth.mutations if x.typ==typ])
    else:
        act = set(ms.mutations)
        exp = set(groundTruth.mutations)

    # Locate putative false positive indel calls, treat them as true positives
    # if they're located close to an expected call of the same type.
    putativeFalseIndels = [x for x in act - exp if x.typ == 'd' or x.typ == 'i']
    trueIndels = groundTruth.mutations

    # Locate closest true indels. If less than 5 bases away and same call, set
    # to true indel call. NOTE: this does not adjust multiple close calls, only
    # the closest.
    for indel in putativeFalseIndels:
        idx = bisect(trueIndels, indel)
        # Watch for boundaries
        if idx == 0:
            idx = 1
        if idx == len(trueIndels):
            idx -= 1

        left = trueIndels[idx-1]
        right = trueIndels[idx]
        nearest = min([left,right], key=lambda m: abs(m.offs-indel.offs))
        if abs(nearest.offs - indel.offs) < 5 and nearest.mut == indel.mut and nearest.typ == indel.typ:
            act.discard(indel)
            act.add(nearest)
            # ammend the gff map to reflect
            ms.gffinf[nearest] = ms.gffinf[indel]

    # confusion matrix
    tp = act & exp
    tpl = len(tp)
    fp = act - exp
    fpl = len(fp)
    fn = exp - act
    tn = seqlen - len(act ^ exp) - tpl
    tpr = tpl/float(tpl+len(fn))
    fpr = fpl/float(fpl+tn)
    fdr = fpl/float(fpl+tpl)

    #print "P: %i N: %i " % (len(exp), seqlen-len(exp)),
    print "TP: %i FP: %i FN: %i TN: %i " % (tpl, fpl, len(fn), tn),
    print "TPR: %f FPR: %f FDR: %f" % (tpr, fpr, fdr)
    if bool(getenv('verbose')):
        print "TP\n".join(["%s\t%s\t"%(str(x),ms.gffinf[x].conf) for x in tp])+'TP'
        print "FP\n".join(["%s\t%s\t"%(str(x),ms.gffinf[x].conf) for x in fp])+'FP'
        print "FN\n".join([str(x) for x in fn])+'FN'

def context():
    '''
    Print out some sequence context around a given mutation offset
        [size=10] yarm context <fasta> <1-based offset>
    '''
    fasta = open(sys.argv[2], 'r')
    # convert offset to 0-based
    offs = int(sys.argv[3]) - 1
    size = int(defenv('size',10))
    fasta.readline()
    sstart = fasta.tell()
    linew = len(fasta.readline()) - 1
    fasta.seek(sstart,0)
    start = offs - size
    fasta.seek(start+(start/linew),1)
    nseq = fasta.read(size*2+1)
    seq = nseq.replace('\n','')
    seq = seq + fasta.read(len(nseq)-len(seq))
    print " ".join([seq[:size], seq[size], seq[size+1:]])

def _invert(stra):
    '''
    Takes a strain and `invert`s the offsets from wild to mutated
    coordinates in 1-based coordinate space.
    '''
    inv = MutatedStrain(stra.id)
    inv.length = stra.length
    ic = 0
    for orig in stra.mutations:
        # flip insertions/deletions, leave substitutions alone
        typ = 'i' if orig.typ == 'd' else 'd' if orig.typ == 'i' else orig.typ
        wild = orig.wild

        # offsets displayed in 1-based coordinates
        offs = orig.offs + ic + 1

        # in this context, deletions are always on the wild strand
        if typ == 'd': wild = '-'

        # offsets for insertions are assigned to the previous base
        if typ == 'i': offs = offs - 1

        inv.mutations.append(Mutation(offs, typ, wild, orig.mut))

        # (de|in)crement as necessary.
        if orig.typ == 'i':
            ic = ic + 1
        elif orig.typ == 'd':
            ic = ic - 1

    return inv

def defarg(i, v):
    return sys.argv[i] if len(sys.argv) == i + 1 else v

def defenv(k, v):
    return getenv(k) if getenv(k) else v

def usage():
    print 'Yet Another Reference Mutator (yarm)'
    print '%s [%s]' % (sys.argv[0], "|".join([k for k in cmds.keys() if k]))

cmds = {None: usage, 'create': create, 'apply': apply,
        'stats': stats, 'dump': dump, 'invert': invert,
        'compare': compare, 'context': context}
hip = ("Awww snap", "Check the 411", "For shizza", "You trippin")

def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else None
    try:
        cmds[cmd]()
    except KeyError:
        usage()
    except IndexError, e:
        print cmds[cmd].__doc__
        if bool(getenv('stack')): raise
    except Exception, e:
        sys.stderr.write("%s! %s\n" % (random.choice(hip), e))
        if bool(getenv('stack')): raise

if __name__ == '__main__':
    main()
