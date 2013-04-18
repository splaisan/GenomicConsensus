Here's a test of the new amplicons support in Quiver.  Previously
Quiver would not see coverage spanning its fixed 500bp windows, and
would no-call the entire genome.  Now the windows extents are
determined adaptively based on where the coverage actually is.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/tinyLambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/tinyLambda/lambdaNEB.fa

  $ quiver --quiet -j${JOBS-8} --noEvidenceConsensusCall=nocall \
  >        $INPUT -r $REFERENCE                                 \
  >        -o variants.gff -o css.fa -o css.fq

  $ grep -v '#' variants.gff
  [1]

  $ fastacomposition css.fa
  css.fa A 287 C 267 G 309 N 47361 T 283

  $ nucmer -mum $REFERENCE css.fa 2>/dev/null

  $ show-aligns out.delta lambda_NEB3011 'lambda_NEB3011|quiver'
  * (glob)
  
  ============================================================
  -- Alignments between lambda_NEB3011 and lambda_NEB3011|quiver
  
  -- BEGIN alignment [ +1 6531 - 6718 | +1 6531 - 6721 ]
  
  
  6531       ctgccgtgcttaagggcaaatacaccatgaccggtgaagccttcgatcc
  6531       ctgccgtgcttaagggcaaatacaccatgaccggtgaagccttcgatcc
                                                              
  
  6580       ggttgaggtg.gatatgggccgcagtgaggagaataacatcacgcagtc
  6580       ggttgaggtgggatatgggccgcagtgaggagaataacatcacgcagtc
                       ^                                      
  
  6628       cggcggcacggagtggagcaagcgtgacaagtccacgtatgaccc.gac
  6629       cggcggcacggagtggagcaagcgtgacaagtccacgtatgacccggac
                                                          ^   
  
  6676       cgacgatatcgaagcctacgcgctga.acgccagcggtgtggtg
  6678       cgacgatatcgaagcctacgcgctgaaacgccagcggtgtggtg
                                       ^                 
  
  
  --   END alignment [ +1 6531 - 6718 | +1 6531 - 6721 ]
  -- BEGIN alignment [ +1 7266 - 7562 | +1 7269 - 7566 ]
  
  
  7266       cctgacggggacgaaagaagaactggcgctccgtgtggcagagctgaaa
  7269       cctgacggggacgaaagaagaactggcgctccgtgtggcagagctgaaa
                                                              
  
  7315       gaggagcttgatgacacggatgaaactgccggtcaggacacccctctca
  7318       gaggagcttgatgacacggatgaaactgccggtcaggacacccctctca
                                                              
  
  7364       gccgggaaaatgtgctgaccggacatgaaa.atgaggtgggatcagcgc
  7367       gccgggaaaatgtgctgaccggacatgaaaaatgaggtgggatcagcgc
                                           ^                  
  
  7412       agccggataccgtgattctggatacgtctgaactggtcacggtcgtggc
  7416       agccggataccgtgattctggatacgtctgaactggtcacggtcgtggc
                                                              
  
  7461       actggtgaagctgcatactgatgcacttcacgccacgcgggatgaacct
  7465       actggtgaagctgcatactgatgcacttcacgccacgcgggatgaacct
                                                              
  
  7510       gtggcatttgtgctgccgggaacggcgtttcgtgtctctgccggtgtgg
  7514       gtggcatttgtgctgccgggaacggcgtttcgtgtctctgccggtgtgg
                                                              
  
  7559       cagc
  7563       cagc
                 
  
  
  --   END alignment [ +1 7266 - 7562 | +1 7269 - 7566 ]
  -- BEGIN alignment [ +1 24760 - 25167 | +1 24764 - 25173 ]
  
  
  24760      tgaaatgatgaagagctctgtgtt.tgtcttcctgcctccagttcgccg
  24764      tgaaatgatgaagagctctgtgttttgtcttcctgcctccagttcgccg
                                     ^                        
  
  24808      ggcattcaacataaaaactgatagcacccggagttccggaaacg..aaa
  24813      ggcattcaacataaaaactgatagcacccggagttccggaaacggaaaa
                                                         ^^   
  
  24855      tttgcatatacccattgctcacgaaaaaaaatgtccttgtcgatatagg
  24862      tttgcatatacccattgctcacgaaaaaa.atgtccttgtcgatatagg
                                          ^                   
  
  24904      gatgaatcgcttggtgtacctcatctactgcgaaaacttgacctttctc
  24910      gatgaatcgcttggtgtacctcatctactgcgaaaacttgacctttctc
                                                              
  
  24953      tcccatattgcagtcgcggcacgatggaactaaattaataggcatcacc
  24959      tcccatattgcagtcgcggcacgatggaactaaattaataggcatcacc
                                                              
  
  25002      gaaaattcaggataatgtgcaataggaagaaaatgatctatattttttg
  25008      gaaaattcaggataatgtgcaataggaagaaaatgatctatattttttg
                                                              
  
  25051      tctgtcctatatcaccacaaaatggacatttttcacctgatgaaacaag
  25057      tctgtcctatatcaccacaaaatggacatttttcacctgatgaaacaag
                                                              
  
  25100      catgtcatcgtaatatgttctagcgggtttgtttttatctcggagatta
  25106      catgtcatcgtaatatgttctagcgggtttgtttttatctcggagatta
                                                              
  
  25149      ttttcataaagcttttcta
  25155      ttttcataaagcttttcta
                                
  
  
  --   END alignment [ +1 24760 - 25167 | +1 24764 - 25173 ]
  -- BEGIN alignment [ +1 30837 - 30950 | +1 30842 - 30954 ]
  
  
  30837      ttttatccggaaactgctgtctggctttttttgatttcagaattagc.c
  30842      ttttatccggaaactgctgtctggcttttt.tgatttcagaattagccc
                                           ^                ^ 
  
  30885      tgacgggcaatgctgcgaagggcgttttcctgctgaggtgtcattgaac
  30890      tgacgcg.gatgctggcaagggcgttttcctgctgaggtgtcattgaac
                  ^ ^^      ^^                                
  
  30934      aagtcccatgtcggcaa
  30938      aagtcccatgtcggcaa
                              
  
  
  --   END alignment [ +1 30837 - 30950 | +1 30842 - 30954 ]
  -- BEGIN alignment [ +1 43907 - 44037 | +1 43911 - 44042 ]
  
  
  43907      aaatttcattcgccaaaaagcccgatgatgagcgactcaccacgggcca
  43911      aaatttcattcgccaaaaagcccgatgatgagcgactcaccacgggcca
                                                              
  
  43956      cggcttctgactctctttccggtactgatgtgatggctgctatggggat
  43960      cggcttctgactctctttccggtactgatgtgatggctgctatggggat
                                                              
  
  44005      ggcgca.atcacaagccggattcggtatggctgc
  44009      ggcgcaaatcacaagccggattcggtatggctgc
                   ^                           
  
  
  --   END alignment [ +1 43907 - 44037 | +1 43911 - 44042 ]
  
  ============================================================




























































