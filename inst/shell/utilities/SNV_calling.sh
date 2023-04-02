#!/bin/bash
while [[ $# -gt 0 ]]
do
    key=$1
    case $key in
        --gtf)
            shift 
            gtf=$1
            shift 
            ;;
        --ref)
            shift 
            ref=$1
            shift 
            ;;
        --input-bam-basename) # including path for bam files, without .bam suffix
            shift 
            bamfilename=$1
            filebasename=`basename $bamfilename`
            path_to_bam=`dirname $1` # this should be basefolder/bam
            basefolder=${path_to_bam::-4}
            shift
            ;;
        --picard) 
            shift
            picard=$1
            shift
            ;;
        --refsnp) 
            shift
            refsnp=$1
            shift
            ;;
        --refindel)
            shift
            refindel=$1
            shift
            ;;
        *)
            echo "Unvalid input for command snv_calling!" >&2
            exit 1
            ;;
    esac
done

# folder path
mkdir -p $path_to_bam
path_to_snv=${basefolder}/snv/
mkdir -p $path_to_snv
path_to_cell_level_snv="$path_to_snv"'cell_level_snv/'
mkdir -p $path_to_cell_level_snv
path_to_time_stats=${basefolder}/time_stats/
mkdir -p $path_to_time_stats
gatk_time_stats="$path_to_time_stats"'time_GATK.csv'
if [ ! -f gatk_time_stats ]; then
    echo "Filename, Time" >> $gatk_time_stats
fi

# check if refsnp and refindel have index
refsnp_idx=${refsnp}.idx
refindel_idx=${refindel}.idx
if [ ! -f "$refsnp_idx" ]; then
    echo "You need to supply index file for $refsnp in same folder!" >&2
    exit 1
fi

if [ ! -f "$refindel_idx" ]; then
    echo "You need to supply index file for $refindel in same folder!" >&2
    exit 1
fi

# For GATK, it requires that reference folder contains .fa .fai .dict
# check if .fa, .fai, .dict are in reference folder
ref_fai="$ref".fai
ref_dict=`echo $ref | sed -e 's/.fa/.dict/g'`
if [ ! -f $ref_fai ]; then
    echo "genome.fa.fai doesn't exist, this code just build a new one."
    samtools faidx "$ref"
fi
if [ ! -f $ref_dict ]; then
    echo "genome.dict doesn't exist, this code just build a new one."
    java -jar $picard CreateSequenceDictionary R=$ref
fi
# Start real work
start_GATK=`date +%s`
# add read group
in_bam="$path_to_bam"/"$filebasename"Aligned.sortedByCoord.out.bam
addrg_bam="$path_to_cell_level_snv""$filebasename"'_sort_rg.bam'
myRGID="$filebasename"'_RGID'
myRGLB=$filebasename
myRGPU=$filebasename
myRGSM=$filebasename
java -jar $picard AddOrReplaceReadGroups \
    I=$in_bam \
    O=$addrg_bam \
    RGID=$myRGID \
    RGLB=$myRGLB \
    RGPL=illumina \
    RGPU=$myRGPU \
    RGSM=$myRGSM
# markduplicates
dedup_bam="$path_to_cell_level_snv""$filebasename"'_sort_rg_dedup.bam'
metrics="$path_to_cell_level_snv""$filebasename"'_metrics.txt'
java -jar $picard MarkDuplicates \
    INPUT=$addrg_bam \
    OUTPUT=$dedup_bam \
    METRICS_FILE=$metrics
# split N CIGAR for RNA
split_bam="$path_to_cell_level_snv""$filebasename"'_sort_rg_dedup_split.bam'
gatk SplitNCigarReads -R $ref -I $dedup_bam -O $split_bam 
# Base Quality Recalibration
recal_table="$path_to_cell_level_snv""$filebasename"'_recal_table.txt'
gatk BaseRecalibrator -R $ref -I $split_bam -O $recal_table \
    --known-sites $refsnp --known-sites $refindel
# PrintReads
recal_reads_bam="$path_to_cell_level_snv""$filebasename"'_sort_rg_dedup_recal.bam'
gatk PrintReads -R $ref -I $split_bam -O $recal_reads_bam 
# ApplyBQSR
gatk_bam="$path_to_cell_level_snv""$filebasename"'_gatk.bam'
gatk ApplyBQSR -I $recal_reads_bam -O $gatk_bam --bqsr-recal-file $recal_table
# HaplotypeCaller
raw_variants="$path_to_cell_level_snv""$filebasename"'_raw_variants.vcf'
gatk HaplotypeCaller -R $ref -I $gatk_bam  --dont-use-soft-clipped-bases \
    -stand-call-conf 20 --native-pair-hmm-threads 8 -O $raw_variants
# VariantFiltration
filtered_variants="$path_to_cell_level_snv""$filebasename"'_filtered.vcf'
gatk VariantFiltration -R $ref -V $raw_variants --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O $filtered_variants
# Filter out PASS item
filtered_pass_variants="$path_to_cell_level_snv""$filebasename"'_filtered_pass.vcf'
cat $filtered_variants | grep -e "#\|PASS"  > $filtered_pass_variants
stop_GATK=`date +%s`
echo $filebasename","$((stop_GATK-start_GATK)) >> $gatk_time_stats
