
Test conversion GFF -> VCF

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/converters/variants.gff.gz

  $ gffToVcf --globalReference=Staphylococcus_aureus_USA300_TCH1516 $INPUT
  ##fileformat=VCFv3.3
  ##fileDate=* (glob)
  ##source=* (glob)
  ##reference=Staphylococcus_aureus_USA300_TCH1516
  ##INFO=NS,1,Integer,"Number of Samples with Data"
  ##INFO=DP,1,Integer,"Total Depth of Coverage"
  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t701415\t.\tG\tD1\t46.00\t0\tNS=1;DP=97 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t970970\t.\tG\tD1\t48.00\t0\tNS=1;DP=97 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1065968\t.\tG\tD1\t49.00\t0\tNS=1;DP=87 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1081288\t.\tG\tD1\t40.00\t0\tNS=1;DP=81 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1315975\t.\tC\tD1\t41.00\t0\tNS=1;DP=100 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1342770\t.\t.\tIA\t49.00\t0\tNS=1;DP=27 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1439019\t.\tC\tD1\t49.00\t0\tNS=1;DP=88 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1456850\t.\tC\tD1\t48.00\t0\tNS=1;DP=84 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1623535\t.\tC\tD1\t47.00\t0\tNS=1;DP=99 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1998595\t.\tG\tD1\t47.00\t0\tNS=1;DP=75 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t2002376\t.\tC\tD1\t48.00\t0\tNS=1;DP=73 (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t2179435\t.\tC\tD1\t48.00\t0\tNS=1;DP=52 (esc)
