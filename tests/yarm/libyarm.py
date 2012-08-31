# Support library
# Yet Another Reference Mutator (yarm), the purpose of which is to help 
# validate existing/new variant algorithms.

from os import getenv
from collections import namedtuple

Mutation = namedtuple('Mutation',['offs','typ','wild','mut'])
GffAnnot = namedtuple('GffAnnot',['conf','cove'])

class MutatedStrain(object):
    '''
    Information about a mutated strain, which can be applied to a reference.
    '''
    def __init__(self, id):
        self._id = id
        self._mutations = []

    @staticmethod
    def fromGff(ins):
        '''
        Parse a GFF formatted input stream into a MutatedStrain. 
        '''
        import re

        ms = MutatedStrain("gff")
        ms.gffinf = dict()
         
        mc = 0 
        for l in ins:
            if l[0] == '#': continue
            f = l.split()
            offs = int(f[3])
            typ = f[2][0].lower()
            try:
                wild = re.search('variantSeq=(\w*)', f[8]).group(1)
            except:
                wild = '-'

            try:
                mut = re.search('reference=(\w*)', f[8]).group(1)
                if typ == 'i': mut = '-'
            except:
                mut = '-'

            conf = re.search('confidence=(\d+)', f[8]).group(1)
            cove = re.search('coverage=(\d+)', f[8]).group(1)

            # Special cases: collapse indel calls into substitutions 
            # NOTE: assuming we're sorted here
            # case 1 (d->i, same index)
            # deletion  10201   reference=G
            # insertion 10201   variantSeq=T
            # -- becomes --
            # Mutation(offs=10201, typ='s', wild='T', mut='G')
            # case 2 (i->d, close index)
            # insertion 5875    variantSeq=C
            # deletion  5877    reference=A
            # -- becomes --
            # Mutation(offs=5877, typ='s', wild='C', mut='A')
            try:
                prev = ms.mutations[mc-1]
                if offs - prev.offs < 3:
                    if prev.typ == 'd' and typ == 'i':
                        mod = prev._replace(typ='s', wild=wild)
                        del ms.gffinf[prev]
                        ms.gffinf[mod] = GffAnnot(conf,cove)
                        ms.mutations[mc-1] = mod
                        continue

                    if prev.typ == 'i' and typ == 'd':
                        mod = prev._replace(offs=offs, typ='s', mut=mut)
                        del ms.gffinf[prev]
                        ms.gffinf[mod] = GffAnnot(conf,cove)
                        ms.mutations[mc-1] = mod
                        continue
            except Exception, e:
                pass
            
            mutobj = Mutation(offs, typ, wild, mut)
            ms.mutations.append(mutobj)
            ms.gffinf[mutobj] = GffAnnot(conf,cove)
            mc = mc + 1
    
        return ms

    @property
    def id(self):
        '''ID of this mutated strain'''
        return self._id
    
    @property
    def length(self):
        '''The length, or total # of bases'''
        return self._length

    @length.setter
    def length(self, value):
        self._length = value

    @property
    def mutations(self):
        '''List of mutations'''
        return self._mutations

    def __repr__(self):
        return "\n%s %i\n%s\n" % (self.id, self.length, 
               "\n".join([self._repstr(x) for x in self.mutations]))

    def _repstr(self, x):
        '''We display offsets in 1-based coordinate space'''
        x._replace(offs=x.offs+1)
        return str(x)


class Mutator(object):
    '''
    Responsible for generating a `mutation strain` based on a `wild-type` 
    reference. The strain can later be applied to a reference.
    '''
    def __init__(self):
        self.bases = 'ATGC'
        # (s)ubstituion, (i)nsertion, (d)eletion
        self.types = getenv('typ') if getenv('typ') else 'sid'
        # in bases
        self.window = int(getenv('window')) if getenv('window') else 100

    def mutate(self, id, instream):
        '''
        Given a stream of bases, build a mutated strain.
        '''
        import random
        strain = MutatedStrain(id)
        seqOffs = 0
        seq = instream.read(self.window)
        # set true if desire mutations at a fixed interval
        fixed = bool(getenv('fixed'))
        mutLen = 0
        while seq:
            seq = seq.replace('\n','')
            seqlen = len(seq)
            if seqlen < 30: break
            # pick a location that avoids the edges of seq
            loc = seqlen/2 if fixed else random.randint(5, seqlen-5)
            wild = seq[loc]
            if wild == 'N': continue
            typ = random.choice(self.types)
            if    typ == 'd':
                mut = '-'
                mutLen = mutLen - 1
            else:
                if typ == 'i': 
                    wild = '-'
                    mutLen = mutLen + 1
                mut = random.choice(self.bases.replace(wild,''))

            strain.mutations.append(Mutation(seqOffs+loc, typ, wild, mut))

            seqOffs = seqOffs + seqlen
            seq = instream.read(self.window)

        strain.length = seqOffs + seqlen + mutLen

        return strain
    
    def apply(self, strain, instream, outstream):
        '''
        Given a stream of bases, apply the mutations in the given strain while
        writing it to the outputstream.
        ''' 
        jc = ' ' if bool(getenv('debug')) else ''
        seqOffs = 0
        smuts = strain.mutations
        mut = smuts.pop(0)
        seq = instream.read(self.window)
        while seq:
            seq = seq.replace('\n','')
            seqlen = len(seq)
            mseq = seq 
            if mut and mut.offs < seqOffs + seqlen:
                # apply a mutation 
                idx = mut.offs - seqOffs
                if   mut.typ == 's':
                    mseq = jc.join([mseq[:idx], mut.mut, mseq[idx+1:]])
                elif mut.typ == 'i':
                    mseq = jc.join([mseq[:idx], mut.mut, mseq[idx:]])
                else:
                    mseq = jc.join([mseq[:idx], mseq[idx+1:]])

                mut = smuts.pop(0) if len(smuts) > 0 else None

            outstream.write(mseq)
            seqOffs = seqOffs + seqlen
            seq = instream.read(self.window)
