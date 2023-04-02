#!/bin/bash
CURRENTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CHECKCOMMANDFUNCTION=${CURRENTPATH}/utilities/check_command.sh

# default variable setting
threadN=1 # default thread
IF_CELL_LIST=false

# help function
myhelp(){
    echo " Arguments:"
    printf "    %-30s  %-35s\n" "--input-dir <folder>" "Folder of Bam files (processed by GATK)"
    printf "    %-30s  %-35s\n" "--output-dir <folder>" "Output folder for read counts result on given locations"
    printf "    %-30s  %-35s\n" "--SNV-list <file>" "SNV locations needed to analysis. Format: chromosome locus locus."
    printf "    %-30s  %-35s\n" "[--valid-cell-list <file>]" "Cells we need to analyze. Each line represents the name for a file. If you don't set it, the code will do analysis on all files in the given input folder"
    printf "    %-30s  %-35s\n" "[--parallel <num>]" "Number of files processed in parallel, optional argument, default value is 1"
    echo
    echo " Your command is like:"
    echo "        $0 --input-dir input_dir --output-dir output_dir --SNV-list SNV_list_file --valid-cell-list valid_cell_list_file --parallel 2"
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
        --ref)
            shift
            ref=$1
            shift
            ;;
        --input-dir) # should be basefolder/snv/cell_level_snv/
            shift
            path_to_gatk_bam=`realpath $1`
            shift
            ;;
        --valid-cell-list)
            shift
            valid_cell_list=$1
            IF_CELL_LIST=true
            shift
            ;;
        --SNV-list)
            shift
            snv_list=$1
            shift
            ;;
        --output-dir) 
            shift
            outputfolder=$1
            outputfolder=`realpath $outputfolder`
            shift
            ;;
        --parallel)
            shift
            threadN=$1
            shift
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

# check argument
if [ -z ${outputfolder+x} ];then
    echo "ERROR in $0: You doesn't supply --output-dir!" >&2
    exit 1
fi
if [ -z ${path_to_gatk_bam+x} ];then
    echo "ERROR in $0: You doesn't supply --input-dir!" >&2
    exit 1
fi
if [ -z ${ref+x} ];then
    echo "ERROR in $0: You doesn't supply --ref!" >&2
    exit 1
fi
if [[ ! -e "$ref" ]]; then
    echo "ERROR in $0: $ref doesn't exist. Please give a valid file with --ref" >&2
    exit 1
fi
if [ -z ${snv_list+x} ];then
    echo "ERROR in $0: You doesn't supply --SNV-list!" >&2
    exit 1
fi
if [[ ! -e "$snv_list" ]]; then
    echo "ERROR in $0: $SNV_list doesn't exist. Please give a valid file with --SNV-list" >&2
    exit 1
fi
if [ -z ${valid_cell_list+x} ];then
    echo
    echo "You doesn't supply --valid-cell-list, code will analyze all files in --input-dir!"
    echo
fi

# check software
if ! $CHECKCOMMANDFUNCTION bamreadcount;then
    exit 1 # no software
fi

# build folder
outputfolder_cell_level=${outputfolder}/cell_level_selected_snv_read_counts/
mkdir -p $outputfolder_cell_level
outputfolder_matrix=${outputfolder}/selected_snv_read_counts_matrix/
mkdir -p $outputfolder_matrix

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

# cellwise-operation
for i in "$path_to_gatk_bam""/*_gatk.bam";
do
    read -u1000
    {
        filename=`echo $i |awk -F/ '{print $NF}' |  awk 'gsub("_gatk.bam","")'`
        outputfile=${outputfolder_cell_level}${filename}.selected_SNV.readcounts
        if $IF_CELL_LIST; then # check if these cells are valid
            if grep -q ${filename} ${valid_cell_list}; then
                bam-readcount -f $ref -l $snv_list $i > ${outputfile}
                sed -e "s/$/\t${filename}/" -i ${outputfile}
            fi
        else
            bam-readcount -f $ref -l $snv_list $i > ${outputfile}
            sed -e "s/$/\t${filename}/" -i ${outputfile}
        fi
        echo >&1000
    } &
done
wait

selected_snv_read_counts_matrix_file=${outputfolder_matrix}selected_SNV_readcounts_matrix.txt
cat ${outputfolder_cell_level}*selected_SNV.readcounts > $selected_snv_read_counts_matrix_file
