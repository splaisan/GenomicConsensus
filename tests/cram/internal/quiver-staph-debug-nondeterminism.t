
  $ export BAM=/mnt/secondary/Share/Quiver/TestData/staph/aligned_reads.bam
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/staph/S_aureus_USA300_TCH1516.fasta

  $ quiver -w0:2140000-2150000 -vv -j1  $BAM -r $REFERENCE -o variants.gff  2> /home/UNIXHOME/dalexander/QuiverLogs/$(date +"%m_%d_%Y__%H:%M:%S").log


  $ grep -v "#" variants.gff | sed 's/\t/ /g'
