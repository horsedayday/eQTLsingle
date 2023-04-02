#!/bin/bash

# define subcommand
CURRENTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MAPPINGFUNCTION=${CURRENTPATH}/utilities/mapping.sh
QUANTIFYFUNCTION=${CURRENTPATH}/utilities/quantify.sh
SNVCALLINGFUNCTION=${CURRENTPATH}/utilities/SNV_calling.sh
MAKEMATRIXFUNCTION=${CURRENTPATH}/utilities/make_matrix.sh
CHECKCOMMANDFUNCTION=${CURRENTPATH}/utilities/check_command.sh

# default variable setting
threadN=1 # default thread
IF_MAPPING=true
IF_CALLING=true

# help function
myhelp(){
    echo "Default operation are:"
    echo "  (1) Build index for mapping software STAR"
    echo "  (2) Mapping from Fastq files to Bam files; Building gene expression matrix"
    echo "  (3) SNV calling from Bam files"
    echo 
    echo "If you want to use the full function, i.e., (1)(2)(3), you need to set the following arguments:"
    printf "    %-30s  %-35s\n" "--ref <file>" "Reference genome" 
    printf "    %-30s  %-35s\n" "--gtf <file>" "Genome annotation" 
    printf "    %-30s  %-35s\n" "--input-fastq-dir <folder>" "Folder of Fastq files" 
    printf "    %-30s  %-35s\n" "--output-dir <folder>" "Output folder to store index of STAR, Bam files, gene expression matrix, SNV calling results, etc." 
    printf "    %-30s  %-35s\n" "--picard <file>" "Path to picard.jar file" 
    printf "    %-30s  %-35s\n" "--refsnp <file>" "Path to reference on SNP, e.g., dbsnp_138.hg19.vcf" 
    printf "    %-30s  %-35s\n" "--refindel <file>" "Path to reference on indel, e.g., 1000G_phase1.indels.hg19.sites.vcf" 
    printf "    %-30s  %-35s\n" "[--parallel <num>]" "Number of files processed in parallel, optional argument, default value is 1" 
    echo 
    echo "If you only want to use part of utilities, you can use ignore argument to ignore specific operation:"
    printf "    %-30s  %-35s\n" "--ignore-build-index" "To ignore the index building step" 
    printf "    %-30s  %-35s\n" "--ignore-mapping" "To ignore the mapping and quantification step" 
    printf "    %-30s  %-35s\n" "--ignore-SNV-calling" "To ignore the SNV calling step" 
    echo 
    echo "For example,"
    echo "  * If you get full results including gene expression matrix and snv matrix from raw reads (fastq files), your command is like:"
    echo "        $0 --ref ref_file --gtf gtf_file --input-fastq-dir fastq_folder --outut-dir output_folder --picard picard_file --refsnp ref_snp_file --refindel ref_indel_file --parallel 2"
    echo "  * If you just want to mapping (including building index), your command is like:"
    echo "        $0 --ref ref_file --gtf gtf_file --input-fastq-dir fastq_folder --outut-dir output_folder --ignore-SNV-calling"
    echo "  * If you just want to do SNV calling, your command is like the command below. Note that, you need to put bam files in output_folder/bam"
    echo "        $0 --ref ref_file --gtf gtf_file --outut-dir output_folder --picard picard_file --refsnp ref_snp_file --refindel ref_indel_file --ignore-build-index --ignore-mapping"
    exit 1
}

if [[ $# == 0 ]]; then
    myhelp
    exit 1
fi

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
        --input-fastq-dir) 
            shift
            path_to_fastq=`realpath $1`
            path_to_fastq=${path_to_fastq}/
            shift
            ;;
        --output-dir) 
            shift
            basefolder=$1
            basefolder=`realpath $basefolder`
            path_to_bam=${basefolder}/bam/
            path_to_index=${basefolder}/star_index/
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
        --parallel)
            shift
            threadN=$1
            shift
            ;;
        --ignore-mapping)
            shift
            IF_MAPPING=false
            ;;
        --ignore-SNV-calling)
            shift
            IF_CALLING=false
            ;;
        --help | -h)
            shift
            myhelp
            exit 1
            ;;
        *)
            echo "Unvalid argument!" >&2
            exit 1
            ;;
    esac
done

# parallel processing
tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -rf $tempfifo

for ((i=1; i<=$threadN; i++))
do
        echo >&1000
done

# check arguments for mapping and building index
if $IF_MAPPING; then
    if [[ ! -e "$ref" ]]; then
        echo "ERROR in building index for mapping: The genome reference you give: $ref doesn't exist." >&2
        exit 1
    fi
    if [ -z ${path_to_fastq+x} ]; then
        echo "ERROR in mapping: You doesn't supply --input-fastq-dir!" >&2
        exit 1
    fi
    if [[ ! -e "$gtf" ]]; then
        echo "ERROR in mapping: The genome annotation you give: $gtf doesn't exist. Please give a valid file with --gtf" >&2
        exit 1
    fi
    if [ -z ${basefolder+x} ]; then 
        echo "ERROR in building index for mapping: You doesn't supply --output-dir!" >&2
        exit 1
    fi
    # check software
    if ! $CHECKCOMMANDFUNCTION STAR; then
        exit 1
    fi
    if ! $CHECKCOMMANDFUNCTION samtools; then
        exit 1
    fi
    if ! $CHECKCOMMANDFUNCTION featureCounts; then
        exit 1
    fi
    echo
    echo "*****  Building index for mapping!  *****"
    echo
    $MAPPINGFUNCTION star_index --gtf $gtf --ref $ref --output $basefolder
fi

if $IF_CALLING; then
    if [ -z ${basefolder+x} ];then
        echo "ERROR in SNV calling: You doesn't supply --output-dir! The code need to search bam files in output-dir/bam/ " >&2
        exit 1
    fi
    if [[ ! -e "$gtf" ]]; then
        echo "ERROR in SNV calling: The genome annotation you give: $gtf doesn't exist. Please give a valid file with --gtf" >&2
        exit 1
    fi
    if [[ ! -e "$ref" ]]; then
        echo "ERROR in SNV calling: The genome reference you give: $ref doesn't exist. Please give a valid file with --ref" >&2
        exit 1
    fi
    if [[ ! -e "$picard" ]]; then
        echo "ERROR in SNV calling: The picard file you give: $picard doesn't exist. Please give a valid file with --picard" >&2
        exit 1
    fi
    if [[ ! -e "$refsnp" ]]; then
        echo "ERROR in SNV calling: The snp reference you give: $refsnp doesn't exist. Please give a valid file with --refsnp" >&2
        exit 1
    fi
    if [[ ! -e "$refindel" ]]; then
        echo "ERROR in SNV calling: The indel reference you give: $refindel doesn't exist. Please give a valid file with --refindel" >&2
        exit 1
    fi
    # check software
    if ! $CHECKCOMMANDFUNCTION java; then
        exit 1
    fi
    if ! $CHECKCOMMANDFUNCTION gatk; then
        exit 1
    fi
    if ! $CHECKCOMMANDFUNCTION picard $picard; then
        exit 1
    fi
    if ! $CHECKCOMMANDFUNCTION python; then
        exit 1
    fi
fi

# start mapping and snv calling
if $IF_MAPPING; then 
    for i in "$path_to_fastq"*_1.fastq;
    do
        read -u1000
        {
            filebasename=`echo $i |awk -F/ '{print $NF}' |  awk 'gsub("_1.fastq","")'`
            fastqfilename=${path_to_fastq}/${filebasename}
            echo
            echo "*****  Mapping for ${filebasename} ...  *****"
            echo
            $MAPPINGFUNCTION star_mapping --input-fastq-basename $fastqfilename --index-folder $path_to_index
            bamfilename=${path_to_bam}${filebasename}
            echo
            echo "*****  Quantifing for ${filebasename} ...  *****"
            echo
            $QUANTIFYFUNCTION --gtf $gtf --input-bam-basename $bamfilename 
            if $IF_CALLING; then
                echo
                echo "*****  SNV Calling for ${filebasename} ...  *****"
                echo
                $SNVCALLINGFUNCTION --gtf $gtf --ref $ref --picard $picard --refsnp $refsnp --refindel $refindel --input-bam-basename $bamfilename
            fi
            echo >&1000
        } &
    done
    wait
    ## generate statistic for mapping and quantification
    echo
    echo "*****  Making gene expression matrix!  *****"
    echo
    $MAKEMATRIXFUNCTION star_statistic --input-dir ${basefolder}
    $MAKEMATRIXFUNCTION featurecounts_statistic --input-dir ${basefolder}


elif $IF_CALLING; then # only snv calling
    BAMSUFFIX="Aligned.sortedByCoord.out.bam"
    for i in "$path_to_bam"*"$BAMSUFFIX";
    do
        read -u1000
        {
            filebasename=`echo $i |awk -F/ '{print $NF}' |  awk -v bamsuffix="$BAMSUFFIX" 'gsub(bamsuffix,"")'`
            bamfilename=${path_to_bam}${filebasename}
            echo
            echo "*****  SNV Calling for ${filebasename} ...  *****"
            echo
            $SNVCALLINGFUNCTION --gtf $gtf --ref $ref --picard $picard --refsnp $refsnp --refindel $refindel --input-bam-basename $bamfilename
            echo >&1000
        } &
    done
    wait
fi

## generate statistic for SNV matrix
if $IF_CALLING; then
    echo
    echo "*****  Making SNV matrix!  *****"
    echo
    path_to_snv_cell_level=${basefolder}/snv/cell_level_snv/
    path_to_snv_matrix=${basefolder}/snv/snv_matrix/
    mkdir -p $path_to_snv_matrix
    $MAKEMATRIXFUNCTION snv_statistics --input-dir $path_to_snv_cell_level --output-dir $path_to_snv_matrix 
fi
