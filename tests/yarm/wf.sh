#!/bin/bash

# This is a simple workflow script that executes some stuff related to variant
# calling.  The output from this can be evaluated using the yarm tool.

# Essentially a decorator function for 'logging' and timing execution of each
# availible command in the workflow
run () {
    desc=$1
    shift
    cmd="$@"
    echo $cmd
    /usr/bin/time -f"${desc} time: %E real, %U user, %S sys" $cmd || exit
}

# use the last successful mainline build
SEYMOUR_HOME=/mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_LastSuccessfulBuild
export SEYMOUR_HOME

# Prepare to launch
. $SEYMOUR_HOME/etc/setup.sh

# In case some, other, poor soul ends up running it.
if [ $# -lt 2 ]
then
    echo "USAGE: $0 <prefix> [range]"
    exit 1
fi

prefix=$1
shift

refpath=$prefix

# works for ../blah/blee, ./blah/blee, ~/blah/blee, /blah/blee
# returns blee
refname=${refpath##[/~.]*/}

# Set of available commands.  They're basically defined in the order they should
# be run, though one can execute any arbitrary range via the [range] command line
# option using standard bash expansion e.g:
# {3..5}: executes commands 3,4,5
#      6: just execute commnd 6
# <none>: executes all commands
# The parameters settings are taken from the smrtpipe modules and should probably
# be updated from time to time.
dsc[0]="blasr"
cmd[0]="compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=3 --noXML --h5mode=w --h5fn=${prefix}.cmp.h5 --minAccuracy=0.75 --minLength=50  -x -minMatch 12 -x -bestn 1 -x -minPctIdentity 70.0 -x -sa ${refpath}/sequence/${refname}.fasta.sa --regionTable=${prefix}filtered.fofn ${prefix}input.fofn ${refpath}"

dsc[1]="loadpls"
cmd[1]="loadPulses ${prefix}input.fofn ${prefix}.cmp.h5 -metrics InsertionQV,IPD,DeletionQV,PulseWidth,QualityValue"

dsc[2]="cmpsort"
cmd[2]="cmpH5Sort.py ${prefix}.cmp.h5"

dsc[3]="gencons"
cmd[3]="variantCaller.py -j4 --algorithm=plurality -o ${prefix}_gencons.gff -r ${refpath}/sequence/${refname}.fasta ${prefix}.cmp.h5"

# NOTE: Evicons is deprecated now. It doesn't exist in recent builds (at this location)
#dsc[4]="evicons"
#cmd[4]="java -javaagent:${SEYMOUR_HOME}/common/lib/java/jyield-0.0.6.jar -jar ${SEYMOUR_HOME}/analysis/lib/java/secondary-analysis-evicons.jar -baseMap 3,4,2,1 -postMin 0.0 -nProc 4 -hdf5Reference ref000001 -variantsFile ${prefix}_evicons.gff -consensusFile out/${prefix}consensus -confidenceFile out/${prefix}confidence -subAlignment -refStart 0 -refEnd 4639651 -hdf5Output out/${prefix}garbage.h5 ${prefix}.cmp.h5"

dsc[4]="bamify"
cmd[4]="SAMIO.py -i -o ${prefix}.sam --bam ${prefix}.cmp.h5" 

dsc[5]="gatk-countcov"
cmd[5]="java -Xmx4g -Djava.io.tmpdir=/scratch -jar $SEYMOUR_HOME/analysis/lib/java/GenomeAnalysisTK.jar -R ${refpath}/sequence/${refname}.fasta -I ${prefix}.bam -knownSites empty.vcf -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -recalFile ${prefix}.recal -nt 6"

dsc[6]="gatk-recal"
cmd[6]="java -Xmx4g -Djava.io.tmpdir=/scratch -jar $SEYMOUR_HOME/analysis/lib/java/GenomeAnalysisTK.jar -R ${refpath}/sequence/${refname}.fasta -I ${prefix}.bam -T TableRecalibration -o ${prefix}_recal.bam -recalFile ${prefix}.recal"

dsc[7]="gatk-unigen"
cmd[7]="java -Xmx4g -Djava.io.tmpdir=/scratch -jar $SEYMOUR_HOME/analysis/lib/java/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${refpath}/sequence/${refname}.fasta -I ${prefix}_recal.bam -o ${prefix}.vcf --output_mode EMIT_VARIANTS_ONLY --max_deletion_fraction 0.55 -nt 6"

ran_gatk=false
for i in ${@-$(seq 0 7)}
do
    if [ "$i" -gt 7 ]
    then
        break
    fi

    run "${dsc[i]}" "${cmd[i]}"

    if [ "$i" -eq 7 ]
    then
        ran_gatk=true
    fi
done

# Extra post-processing only run this if we ran gatk
if $ran_gatk
then
    echo "vcfToGff.py ${prefix}.vcf > ${prefix}_gatk.gff"
    /usr/bin/time -f"vcf2gff time: %E real, %U user, %S sys" vcfToGff.py ${prefix}.vcf > ${prefix}_gatk.gff
fi
