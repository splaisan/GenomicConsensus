Get the arrow evidence dump and make sure it can be grokked.

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/all4mer/out.aligned_subreads.bam
  $ export REFERENCE=$DATA/all4mer/All4mer.V2.01_Insert.fa

Run arrow w/ evidence dump

  $ arrow --dumpEvidence=all $INPUT -r $REFERENCE -o v.gff -o css.fa -o css.fq

Inspect the output...

  $ find evidence_dump | sort
  evidence_dump
  evidence_dump/All4mer.V2.01_Insert
  evidence_dump/All4mer.V2.01_Insert/0-260
  evidence_dump/All4mer.V2.01_Insert/0-260/arrow-scores.h5
  evidence_dump/All4mer.V2.01_Insert/0-260/consensus.fa


Try to load it up using the API...

  $ python << EOF
  > from GenomicConsensus.arrow.evidence import ArrowEvidence
  > ev = ArrowEvidence.load("evidence_dump/All4mer.V2.01_Insert/0-260")
  > assert 8*len(ev.consensus)==len(ev.colNames)==2080
  > assert ev.delta.shape == (len(ev.rowNames), len(ev.colNames))==(95,2080)
  > EOF
