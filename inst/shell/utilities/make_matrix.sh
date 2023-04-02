#!/bin/bash

CURRENTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
snv_matrix_py=${CURRENTPATH}/snv_matrix.py

star_statistic(){
    while [[ $# -gt 0 ]]
    do
        key=$1
        case $key in 
            --input-dir) # basefolder
                shift
                basefolder=$1
                basefolder=`realpath $basefolder`
                shift
                ;;
            *)
                echo "Unvalid input for command star_statistic!" >&2
                exit 1
                ;;
        esac
    done

    path_to_output=${basefolder}/mapping_stat/
    mkdir -p  $path_to_output
    outputfile="$path_to_output"'mapping_stat.csv'
    if [ -f $outputfile ]; then
        echo "rewrite a new CSV file!"
        rm $outputfile
    fi
    echo "Filename,StartTime,FinishTime,ReadLength,NumReads,Unique,Unique_Percent,Multi,Multi_Percent" >> $outputfile

    path_to_bam=${basefolder}/bam/
    for i in "$path_to_bam"*Log.final.out;
    do
        Filename=`echo $i | awk -F/ '{print $NF}'|awk '{gsub("Log.final.out", "");print}'`
        StartTime=`cat $i | grep "Started mapping on" | awk -F"|" '{print $NF}'`
        FinishTime=`cat $i | grep "Finished on" | awk -F"|" '{print $NF}'`
        ReadLength=`cat $i | grep "Number of input reads" | awk -F"|" '{print $NF}'`
        NumReads=`cat $i | grep "Average input read length" | awk -F"|" '{print $NF}'` 
        Unique=`cat $i | grep "Uniquely mapped reads number" | awk -F"|" '{print $NF}'` 
        Unique_Percent=`cat $i | grep "Uniquely mapped reads %" | awk -F"|" '{print $NF}'` 
        Multi=`cat $i | grep "Number of reads mapped to multiple loci" | awk -F"|" '{print $NF}'` 
        Multi_Percent=`cat $i | grep "% of reads mapped to multiple loci" | awk -F"|" '{print $NF}'` 
        echo $Filename","$StartTime","$FinishTime","$ReadLength","$NumReads","$Unique","$Unique_Percent","$Multi","\
            $Multi_Percent>> $outputfile
    done
}

featurecounts_statistic(){
    while [[ $# -gt 0 ]]
    do
        key=$1
        case $key in 
            --input-dir) # basefolder
                shift
                basefolder=$1
                basefolder=`realpath $basefolder`
                shift
                ;;
            *)
                echo "Unvalid input for command star_statistic!" >&2
                exit 1
                ;;
        esac
    done

    path_to_input=${basefolder}/quantification_stats/cell_level_files/
    path_to_output=${basefolder}/quantification_stats/gene_expression_matrix/
    outputfile_gene_expression="$path_to_output"gene_expression_matrix.csv
    outputfile_featurecounts_summary="$path_to_output"featurecounts_summary_matrix.csv
    # to check if input folder exist
    if [ ! -d "$path_to_input" ]; then
        echo 'Folder '"$path_to_output"' does not exist! Can not make gene expression profile.'
        exit 1
    fi
    # make output folder
    mkdir -p $path_to_output

    internalfile=${path_to_output}expressionMatrix.tmp
    indexfile=${path_to_output}wholeindex.tmp
    cnt=0
    suffix=".featurecounts"
    for i in "$path_to_input"/*"$suffix" ;
    do
        filename=`basename $i`
        filebasename=${filename%"$suffix"}
        ((cnt+=1))
        if [[ "$cnt" == 1 ]];then
            cat $i | awk -v OFS="," 'FNR > 2 {print $1}' > $indexfile
            cat $i | awk -v OFS="," -v Fname=${filebasename} 'BEGIN {print "Geneid," Fname} FNR > 2 {print $1,$NF}' > $outputfile_gene_expression
        else
            currentfilename=`basename $i`
            tmpindexfile=${path_to_output}${currentfilename}.index.tmp
            cat $i | awk -v OFS="," 'FNR > 2 {print $1}' > $tmpindexfile
            if cmp -s "$indexfile" "$tmpindexfile"; then
                tmpcountfile=${path_to_output}${currentfilename}.count.tmp
                cat $i | awk -v OFS="," -v Fname=${filebasename} 'BEGIN {print Fname} FNR > 2 {print $NF}' > $tmpcountfile
                paste -d, ${outputfile_gene_expression} ${tmpcountfile} > $internalfile
                mv $internalfile $outputfile_gene_expression
                rm $tmpindexfile
                rm $tmpcountfile
            else
                echo "ERROR: It seems that you use different annotation files for different cells. Please uniform it or make gene expression matrix manually!" >&2
                exit 1
            fi
        fi
    done

    rm $indexfile
    
    echo "Gene expression matrix has been built!"

    # feature summary matrix
    if [ ! -f "$outputfile_featurecounts_summary" ]; then
        echo 'Filename,Assigned,Unassigned_Unmapped,Unassigned_Read_Type,Unassigned_Singleton,Unassigned_MappingQuality,Unassigned_Chimera,Unassigned_FragmentLength,Unassigned_Duplicate,Unassigned_MultiMapping,Unassigned_Secondary,Unassigned_NonSplit,Unassigned_NoFeatures,Unassigned_Overlapping_Length,Unassigned_Ambiguity' >> "$outputfile_featurecounts_summary"
    fi
    for i in "$path_to_input"*.featurecounts.summary ;
    do
        Filename=`echo $i | awk -F/ '{print $NF}'|awk '{gsub(".featurecounts.summary", "");print}'`
        Assigned=`cat $i | grep "Assigned" | awk '{print $NF}'`
        Unassigned_Unmapped=`cat $i | grep "Unassigned_Unmapped" | awk '{print $NF}'`
        Unassigned_Read_Type=`cat $i | grep "Unassigned_Read_Type" | awk '{print $NF}'`
        Unassigned_Singleton=`cat $i | grep "Unassigned_Singleton" | awk '{print $NF}'`
        Unassigned_MappingQuality=`cat $i | grep "Unassigned_MappingQuality" | awk '{print $NF}'`
        Unassigned_Chimera=`cat $i | grep "Unassigned_Chimera" | awk '{print $NF}'`
        Unassigned_FragmentLength=`cat $i | grep "Unassigned_FragmentLength" | awk '{print $NF}'`
        Unassigned_Duplicate=`cat $i | grep "Unassigned_Duplicate" | awk '{print $NF}'`
        Unassigned_MultiMapping=`cat $i | grep "Unassigned_MultiMapping" | awk '{print $NF}'`
        Unassigned_Secondary=`cat $i | grep "Unassigned_Secondary" | awk '{print $NF}'`
        Unassigned_NonSplit=`cat $i | grep "Unassigned_NonSplit" | awk '{print $NF}'`
        Unassigned_NoFeatures=`cat $i | grep "Unassigned_NoFeatures" | awk '{print $NF}'`
        Unassigned_Overlapping_Length=`cat $i | grep "Unassigned_Overlapping_Length" | awk '{print $NF}'`
        Unassigned_Ambiguity=`cat $i | grep "Unassigned_Ambiguity" | awk '{print $NF}'`
        echo $Filename","$Assigned","$Unassigned_Unmapped","$Unassigned_Read_Type","$Unassigned_Singleton","$Unassigned_MappingQuality","$Unassigned_Chimera","$Unassigned_FragmentLength","$Unassigned_Duplicate","$Unassigned_MultiMapping","$Unassigned_Secondary","$Unassigned_NonSplit","$Unassigned_NoFeatures","$Unassigned_Overlapping_Length","$Unassigned_Ambiguity >> "$outputfile_featurecounts_summary"
    done

    echo "Feature summary matrix has been built!"

}

snv_statistics(){
threadN=1 # default thread
while [[ $# -gt 0 ]]
do
    key=$1
    case $key in
        --input-dir) # basefolder/snv/cell_level_snv
            shift
            inputfolder=$1
            inputfolder=`realpath $inputfolder`
            shift
            ;;
        --output-dir)
            shift
            outputfolder=$1
            outputfolder=`realpath $outputfolder`
            shift
            ;;
        *)
            echo "Unvalid argument!" >&2
            exit 1
            ;;
    esac
done

# check argument
if [ -z ${inputfolder+x} ]; then
    echo "ERROR in building SNV matrix: You doesn't supply --input-dir for $0!" >&2
    exit 1
fi
if [ -z ${outputfolder+x} ]; then
    echo "ERROR in building SNV matrix: You doesn't supply --output-dir for $0!" >&2
    exit 1
fi

python $snv_matrix_py $inputfolder $outputfolder

}

subcommand=$1
case $subcommand in 
    star_statistic)
        shift
        star_statistic $@
        ;;
    featurecounts_statistic)
        shift
        featurecounts_statistic $@
        ;;
    snv_statistics)
        shift
        snv_statistics $@
        ;;
    *)
        echo "ERROR: bugs in $0! check it!" >&2
        exit 1
        ;;
esac
