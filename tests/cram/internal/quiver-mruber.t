
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/mruber/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/mruber/Mruber_DSM_1279.fasta
  $ quiver -p C2 -j${JOBS-8} $INPUT -r $REFERENCE -o variants.gff -o css.fasta

Inspect the variant calls.

  $ grep -v "#" variants.gff | sed 's/\t/ /g'
  M.ruber . substitution 357364 357364 . . . reference=C;variantSeq=T;coverage=100;confidence=47
  M.ruber . insertion 640716 640716 . . . reference=.;variantSeq=C;coverage=100;confidence=49
  M.ruber . insertion 1320669 1320669 . . . reference=.;variantSeq=C;coverage=100;confidence=49
  M.ruber . deletion 1878514 1878953 . . . reference=AGGGCGTACTTCTTTTCGGGTGCAGATGCGTAGGCATCGTAGTTGAACAGGGTTTTGACCGCCATTGAGCACTCCTTTTACGGTTCCACAATGAGTTTGCTGATCATGTTGGCGTGGCCGATGCCGCAGTATTCGTTGCAGATGATGGGATACTCACCGGGTTTGCTGAAGGTGTAGCTGACCTTGGCAATTTCCCCCGGTATCACCTGTACGTTGATGTTGGTGTTGTGTACGTGGAAGCTGTGCTGCACATCGGGTGAGGTGATATAGAAGGTTACCTTCCTGCCCACCTTGAACCGCATCTCCGCTGGCAGGTAGCCAAAGGCAAAGGCCTGCACATAGGCCACGTACTCGTTGCCGACCTGCTCAACCCGTGGGTTGGCAAAGTCTCCCTCGGTGCGCACCTTGGTGGCGTCGATGCGGCCTGCCCCCACCGGGTT;variantSeq=.;coverage=100;confidence=50
  M.ruber . insertion 1987969 1987969 . . . reference=.;variantSeq=G;coverage=100;confidence=49
  M.ruber . insertion 2010700 2010700 . . . reference=.;variantSeq=T;coverage=100;confidence=48
  M.ruber . insertion 2070035 2070035 . . . reference=.;variantSeq=A;coverage=100;confidence=49
  M.ruber . insertion 2827713 2827713 . . . reference=.;variantSeq=T;coverage=100;confidence=49
  M.ruber . deletion 2841287 2841301 . . . reference=AAGCACGCCGAGGGA;variantSeq=.;coverage=100;confidence=49


The variant calls have all been Sanger validated!

	|                  | Confirmed | Confirmed  |
	| Variant call     | by eye?   | by Sanger? |
	|------------------+-----------+------------|
	| 357364 C>T       | YES       | YES        |
	| 640716 InsC      | YES       | YES        |
	| 1320669 InsC     | YES       | YES        |
	| 1878514 Del440bp | YES       | YES        |
	| 1987969 InsG     | YES       | YES        |
	| 2010700 InsT     | YES       | YES        |
	| 2070035 InsA     | YES       | YES        |
	| 2827713 InsT     | YES       | YES        |
	| 2841287 Del15bp  | YES       | YES        |


Look at the consensus output.

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null

There are two gaps corresponding to the structural deletions:

  $ show-diff -H out.delta | sed 's/\t/ /g'
  M.ruber GAP 1878514 1878953 440 0 440
  M.ruber GAP 2851299 2831302 -19996 -19981 -15

... and there are some SNPS.  Five of them are at the coverage desert
before the large deletion, seven are accounted for in the
variants.gff, and the remaining three are low-confidence miscalls.

  $ show-snps -H -C out.delta
    233298   C .   233297    |   124066   233297  |  1  1  M.ruber\tM.ruber|quiver (esc)
    357364   C T   357363    |   124066   357363  |  1  1  M.ruber\tM.ruber|quiver (esc)
    640719   . C   640719    |   283355   640719  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1320671   . C   1320672   |   299698  1320671  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1620369   C .   1620369   |   258131  1476643  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1878500   . C   1878501   |        0  1218511  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1878500   . C   1878502   |        0  1218510  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1878500   . C   1878503   |        0  1218509  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1878500   . A   1878504   |        0  1218508  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1878500   . G   1878505   |        0  1218507  |  1  1  M.ruber\tM.ruber|quiver (esc)
   1987973   . G   1987539   |    22731  1109473  |  1  1  M.ruber\tM.ruber|quiver (esc)
   2010704   . T   2010271   |    22731  1086741  |  1  1  M.ruber\tM.ruber|quiver (esc)
   2070035   . A   2069603   |    59331  1027409  |  1  1  M.ruber\tM.ruber|quiver (esc)
   2827716   . T   2827285   |    13583   269727  |  1  1  M.ruber\tM.ruber|quiver (esc)
