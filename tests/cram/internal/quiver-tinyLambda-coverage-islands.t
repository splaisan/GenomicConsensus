Here's a test of the new amplicons support in Quiver.  Previously
Quiver would not see coverage spanning its fixed 500bp windows, and
would no-call the entire genome.  Now the windows extents are
determined adaptively based on where the coverage actually is.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/tinyLambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/tinyLambda/lambdaNEB.fa

  $ quiver -p unknown --quiet -j${JOBS-8} --noEvidenceConsensusCall=nocall \
  >        $INPUT -r $REFERENCE                                 \
  >        -o variants.gff -o css.fa -o css.fq


These variant calls actually look reasonable given the reads, but the
confidences are too high.  Fix this.

  $ grep -v '#' variants.gff
  lambda_NEB3011\t.\tinsertion\t24781\t24781\t.\t.\t.\treference=.;variantSeq=T;coverage=6;confidence=57 (esc)
  lambda_NEB3011\t.\tdeletion\t24878\t24878\t.\t.\t.\treference=A;variantSeq=.;coverage=16;confidence=43 (esc)
  lambda_NEB3011\t.\tinsertion\t30882\t30882\t.\t.\t.\treference=.;variantSeq=C;coverage=5;confidence=44 (esc)

  $ fastacomposition css.fa
  css.fa A 282 C 266 G 305 N 47361 T 281

  $ nucmer -mum $REFERENCE css.fa 2>/dev/null

  $ show-aligns out.delta lambda_NEB3011 'lambda_NEB3011|quiver'
  * (glob)
  
  ============================================================
  -- Alignments between lambda_NEB3011 and lambda_NEB3011|quiver
  
  -- BEGIN alignment [ +1 6531 - 6718 | +1 6531 - 6718 ]
  
  
  6531       ctgccgtgcttaagggcaaatacaccatgaccggtgaagccttcgatcc
  6531       ctgccgtgcttaagggcaaatacaccatgaccggtgaagccttcgatcc
                                                              
  
  6580       ggttgaggtggatatgggccgcagtgaggagaataacatcacgcagtcc
  6580       ggttgaggtggatatgggccgcagtgaggagaataacatcacgcagtcc
                                                              
  
  6629       ggcggcacggagtggagcaagcgtgacaagtccacgtatgacccgaccg
  6629       ggcggcacggagtggagcaagcgtgacaagtccacgtatgacccgaccg
                                                              
  
  6678       acgatatcgaagcctacgcgctgaacgccagcggtgtggtg
  6678       acgatatcgaagcctacgcgctgaacgccagcggtgtggtg
                                                      
  
  
  --   END alignment [ +1 6531 - 6718 | +1 6531 - 6718 ]
  -- BEGIN alignment [ +1 7266 - 7562 | +1 7266 - 7561 ]
  
  
  7266       cctgacggggacgaaagaagaactggcgctccgtgtggcagagctgaaa
  7266       cctgacggggacgaaagaagaactggcgctccgtgtggcagagctgaaa
                                                              
  
  7315       gaggagcttgatgacacggatgaaactgccggtcaggacacccctctca
  7315       gaggagcttgatgacacggatgaaactgccggtcaggacacccctctca
                                                              
  
  7364       gccgggaaaatgtgctgaccggacatgaaaatgaggtgggatcagcgca
  7364       gccgggaaaatgtgctgaccggacatgaaaatga.gtgggatcagcgca
                                               ^              
  
  7413       gccggataccgtgattctggatacgtctgaactggtcacggtcgtggca
  7412       gccggataccgtgattctggatacgtctgaactggtcacggtcgtggca
                                                              
  
  7462       ctggtgaagctgcatactgatgcacttcacgccacgcgggatgaacctg
  7461       ctggtgaagctgcatactgatgcacttcacgccacgcgggatgaacctg
                                                              
  
  7511       tggcatttgtgctgccgggaacggcgtttcgtgtctctgccggtgtggc
  7510       tggcatttgtgctgccgggaacggcgtttcgtgtctctgccggtgtggc
                                                              
  
  7560       agc
  7559       agc
                
  
  
  --   END alignment [ +1 7266 - 7562 | +1 7266 - 7561 ]
  -- BEGIN alignment [ +1 24760 - 25167 | +1 24759 - 25166 ]
  
  
  24760      tgaaatgatgaagagctctgtgtt.tgtcttcctgcctccagttcgccg
  24759      tgaaatgatgaagagctctgtgttttgtcttcctgcctccagttcgccg
                                     ^                        
  
  24808      ggcattcaacataaaaactgatagcacccggagttccggaaacgaaatt
  24808      ggcattcaacataaaaactgatagcacccggagttccggaaacgaaatt
                                                              
  
  24857      tgcatatacccattgctcacgaaaaaaaatgtccttgtcgatataggga
  24857      tgcatatacccattgctcacgaaaaaa.atgtccttgtcgatataggga
                                        ^                     
  
  24906      tgaatcgcttggtgtacctcatctactgcgaaaacttgacctttctctc
  24905      tgaatcgcttggtgtacctcatctactgcgaaaacttgacctttctctc
                                                              
  
  24955      ccatattgcagtcgcggcacgatggaactaaattaataggcatcaccga
  24954      ccatattgcagtcgcggcacgatggaactaaattaataggcatcaccga
                                                              
  
  25004      aaattcaggataatgtgcaataggaagaaaatgatctatattttttgtc
  25003      aaattcaggataatgtgcaataggaagaaaatgatctatattttttgtc
                                                              
  
  25053      tgtcctatatcaccacaaaatggacatttttcacctgatgaaacaagca
  25052      tgtcctatatcaccacaaaatggacatttttcacctgatgaaacaagca
                                                              
  
  25102      tgtcatcgtaatatgttctagcgggtttgtttttatctcggagattatt
  25101      tgtcatcgtaatatgttctagcgggtttgtttttatctcggagattatt
                                                              
  
  25151      ttcataaagcttttcta
  25150      ttcataaagcttttcta
                              
  
  
  --   END alignment [ +1 24760 - 25167 | +1 24759 - 25166 ]
  -- BEGIN alignment [ +1 30837 - 30950 | +1 30835 - 30945 ]
  
  
  30837      ttttatccggaaactgctgtctggctttttttgatttcagaattag.cc
  30835      ttttatccggaaactgctgtctggcttttt.tgatttcagaa.tagccc
                                           ^           ^   ^  
  
  30885      tgacgggcaatgctgcgaagggcgttttcctgctgaggtgtcattgaac
  30882      tgacgcg.gatgctgcgaagggcgttttcctgctgagg.gtcattgaac
                  ^ ^^                             ^          
  
  30934      aagtcccatgtcggcaa
  30929      aagtcccatgtcggcaa
                              
  
  
  --   END alignment [ +1 30837 - 30950 | +1 30835 - 30945 ]
  -- BEGIN alignment [ +1 43908 - 44037 | +1 43902 - 44030 ]
  
  
  43908      aatttcattcgccaaaaagcccgatgatgagcgactcaccacgggccac
  43902      aatttcattcgccaaaaagc.cgatgatgagcgactcaccacgggccac
                                 ^                            
  
  43957      ggcttctgactctctttccggtactgatgtgatggctgctatggggatg
  43950      ggcttctgactctctttccggtactgatgtgatggctgctatggggatg
                                                              
  
  44006      gcgcaatcacaagccggattcggtatggctgc
  43999      gcgcaatcacaagccggattcggtatggctgc
                                             
  
  
  --   END alignment [ +1 43908 - 44037 | +1 43902 - 44030 ]
  
  ============================================================




























































