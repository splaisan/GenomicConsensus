  $ export PATH=$TESTDIR/../yarm:$PATH
  $ export DATA=$TESTDIR/../data

Substitutions

  $ for i in $(seq 2 2 10);do typ=s id=9_ecoli_mutated yarm.py compare $DATA/yarm/strain.pkl $DATA/yarm/9_ecoli_mutated_${i}.gff;done
  TP: 47 FP: 2 FN: 27 TN: 4639646  TPR: 0.635135 FPR: 0.000000 FDR: 0.040816
  TP: 64 FP: 0 FN: 10 TN: 4639648  TPR: 0.864865 FPR: 0.000000 FDR: 0.000000
  TP: 68 FP: 0 FN: 6 TN: 4639648  TPR: 0.918919 FPR: 0.000000 FDR: 0.000000
  TP: 68 FP: 1 FN: 6 TN: 4639647  TPR: 0.918919 FPR: 0.000000 FDR: 0.014493
  TP: 67 FP: 0 FN: 7 TN: 4639648  TPR: 0.905405 FPR: 0.000000 FDR: 0.000000


Insertions

  $ for i in $(seq 2 2 10);do typ=i id=9_ecoli_mutated yarm.py compare $DATA/yarm/strain.pkl $DATA/yarm/9_ecoli_mutated_${i}.gff;done
  TP: 56 FP: 19 FN: 30 TN: 4639617  TPR: 0.651163 FPR: 0.000004 FDR: 0.253333
  TP: 65 FP: 1 FN: 21 TN: 4639635  TPR: 0.755814 FPR: 0.000000 FDR: 0.015152
  TP: 64 FP: 1 FN: 22 TN: 4639635  TPR: 0.744186 FPR: 0.000000 FDR: 0.015385
  
      Compares a ground-truth strain to a variations file
          [typ=[s|i|d]] [id=] [verbose=0] yarm compare <strains> <variants>
      
  TP: 68 FP: 1 FN: 18 TN: 4639635  TPR: 0.790698 FPR: 0.000000 FDR: 0.014493


Deletions

  $ for i in $(seq 2 2 10);do typ=d id=9_ecoli_mutated yarm.py compare $DATA/yarm/strain.pkl $DATA/yarm/9_ecoli_mutated_${i}.gff;done
  TP: 46 FP: 17 FN: 24 TN: 4639635  TPR: 0.657143 FPR: 0.000004 FDR: 0.269841
  TP: 55 FP: 6 FN: 15 TN: 4639646  TPR: 0.785714 FPR: 0.000001 FDR: 0.098361
  TP: 54 FP: 4 FN: 16 TN: 4639648  TPR: 0.771429 FPR: 0.000001 FDR: 0.068966
  TP: 53 FP: 5 FN: 17 TN: 4639647  TPR: 0.757143 FPR: 0.000001 FDR: 0.086207
  TP: 51 FP: 4 FN: 19 TN: 4639648  TPR: 0.728571 FPR: 0.000001 FDR: 0.072727

