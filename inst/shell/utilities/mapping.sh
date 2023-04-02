#!/bin/bash
star_index(){
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
            --output)
                shift 
                path_to_star_index=`realpath $1`
                path_to_star_index=${path_to_star_index}/star_index/
                shift
                ;;
            *)
                echo "Unvalid input for command star_index!" >&2
                exit 1
                ;;
        esac
    done    
    mkdir -p "$path_to_star_index"
    STAR --runThreadN 30  --runMode genomeGenerate \
        --genomeDir $path_to_star_index \
        --genomeFastaFiles $ref \
        --sjdbGTFfile $gtf \
        --sjdbOverhang 100  
    mv Log.out $path_to_star_index
}


star_mapping(){
    while [[ $# -gt 0 ]] 
    do
        key=$1
        case $key in
            --input-fastq-basename) # path + basename without suffix
                shift 
                fastqfilename=$1 
                filebasename=`basename $fastqfilename`
                shift 
                ;;
            --index-folder)
                shift 
                path_to_star_index=`realpath $1` # should be path/star_index
                basefolder=${path_to_star_index::-11}
                path_to_bam=${basefolder}/bam/
                mkdir -p $path_to_bam
                shift
                ;;
            *)
                echo "Unvalid input for command star_mapping!" >&2
                exit 1
                ;;
        esac
    done    
    path_to_time_stats=${basefolder}/time_stats/
    mkdir -p $path_to_time_stats
    fastq1="$fastqfilename"_1.fastq
    fastq2="$fastqfilename"_2.fastq
    output="$path_to_bam""$filebasename"
    start_STAR=`date +%s`
    STAR --runThreadN 10 --genomeDir "$path_to_star_index" \
        --readFilesIn $fastq1 $fastq2 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $output \
        --quantMode TranscriptomeSAM 
    samtools index -@ 10 ""$output"Aligned.sortedByCoord.out.bam"

    # this part is for quantmode, we need to sort and index it
    samtools sort -@ 10 ""$output"Aligned.toTranscriptome.out.bam" \
         -T ""$output"tmp" \
         -o ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
    samtools index -@ 10 ""$output"Aligned.toTranscriptome.sortedByCoord.out.bam"
    stop_STAR=`date +%s`
    time_file="$path_to_time_stats"time_STAR.csv
    # check if time file exist. if not, add header
    if [ ! -f $time_file ]; then
        echo 'Filename,Time' >> $time_file
    fi
    echo $filebasename","$((stop_STAR-start_STAR)) >> $time_file 
}


subcommand=$1
case $subcommand in 
    star_index)
        shift 
        star_index $@
        ;;
    star_mapping)
        shift 
        star_mapping $@
        ;;
    *)
        echo "ERROR: bugs in $0! check it!" >&2
        exit 1
        ;;
esac
