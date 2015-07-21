
Test Quiver invocation using resolved tool contract (basically identical to quiver-hcv.t with different input method)

  $ export DATA=$TESTDIR/../data
  $ export TC=$DATA/resolved_tool_contract.json
  $ cp -r $DATA/hcv .

  $ python -m GenomicConsensus.main --emit-tool-contract | grep -c GFF
  2

  $ python -m GenomicConsensus.main --resolved-tool-contract $TC
  Attempting to Load resolved tool contract from ['--resolved-tool-contract', *] (glob)

  $ cat v.gff
  ##gff-version 3
  ##pacbio-variant-version 2.1
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##source-alignment-file * (glob)
  ##source-reference-file * (glob)
  ##sequence-region 5primeEnd 1 156
  ##sequence-region 3primeEnd 1 386
  3primeEnd\t.\tdeletion\t296\t296\t.\t.\t.\treference=G;variantSeq=.;coverage=92;confidence=4 (esc)

  $ cat css.fasta
  >5primeEnd|quiver
  GGAACCGGTGAGTACACCGGAATTGCCAGGACGACCGGGTCCTTTCGTGGATAAACCCGC
  TCAATGCCTGGAGATTTGGGCGTGCCCCCGCAAGACTGCTAGCCGAGTAGTGTTGGGTCG
  CGAAAGGCCTTGTGGTACTGCCTGATAGGGTGCTTG
  >3primeEnd|quiver
  TACCTGGTCATAGCCTCCGTGAAGGCTCTCAGGCTCGCTGCATCCTCCGGGACTCCCTGA
  CTTTCACAGATAACGACTAAGTCGTCGCCACACACGAGCATGGTGCAGTCCTGGAGCCCA
  GCGGCTCGACAGGCTGCTTTGGCCTTGATGTAGCAGGTGAGGGTGTTACCACAGCTGGTC
  GTCAGTACGCCGCTCGCGCGGCACCTGCGATAGCCGCAGTTTTCCCCCCTTGAATTAGTA
  AGAGGGCCCCCGACATAGAGCCTCTCGGTGAGGGACTTGATGGCCACGCGGGCTTGGGGT
  CCAGGTCACAACATTGGTAAATTGCCTCCTCTGTACGGATATCGCTCTCAGTGACTGTGG
  AGTCAAAGCAGCGGGTATCATACGA
