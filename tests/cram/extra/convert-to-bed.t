
Test conversion variants GFF -> BED.

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/converters/variants.gff.gz

  $ gffToBed --name=variants \
  >          --description="PacBio variant calls" \
  > variants $INPUT
  track name=variants description="PacBio variant calls" useScore=0
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t701414\t701415\t701415del\t46.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t970969\t970970\t970970del\t48.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1065967\t1065968\t1065968del\t49.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1081287\t1081288\t1081288del\t40.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1315974\t1315975\t1315975del\t41.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1342769\t1342770\t1342770_1342771insA\t49.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1439018\t1439019\t1439019del\t49.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1456849\t1456850\t1456850del\t48.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1623534\t1623535\t1623535del\t47.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t1998594\t1998595\t1998595del\t47.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t2002375\t2002376\t2002376del\t48.000\t. (esc)
  gi|160367075|gb|CP000730.1| Staphylococcus aureus subsp. aureus USA300_TCH1516, complete genome\t2179434\t2179435\t2179435del\t48.000\t. (esc)
