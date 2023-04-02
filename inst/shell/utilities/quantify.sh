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
        --input-bam-basename) # path + basename
            shift
            bamfilename=$1
            path_to_bam=`dirname $1` # path to bam
            basefolder=${path_to_bam::-4}
            filebasename=`basename $1`
            shift
            ;;
        *)
            echo "Unvalid input for quantify!" >&2
            exit 1
            ;;
    esac
done

path_to_time_stats=${basefolder}/time_stats/
mkdir -p $path_to_time_stats
subfix='Aligned.sortedByCoord.out.bam'
bamfile=${bamfilename}${subfix}
threads=20
path_to_quantification_stats=${basefolder}/quantification_stats/cell_level_files/
mkdir -p $path_to_quantification_stats
output=${path_to_quantification_stats}${filebasename}.featurecounts
start_featurecounts=`date +%s`
featureCounts -T $threads -p -a $gtf -o $output $bamfile 
stop_featurecounts=`date +%s`
time_file="$path_to_time_stats"time_featurecounts.csv
# check if time file of feature counts exist, if not, add header
if [ ! -f "$time_file" ]; then
    echo 'Filename,Time' >> $time_file
fi
echo ${filebasename}","$((stop_featurecounts-start_featurecounts)) >> $time_file
