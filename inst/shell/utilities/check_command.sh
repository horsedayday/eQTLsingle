#!/bin/bash

## STAR
check_STAR(){
    if ! command -v STAR &> /dev/null
    then
        echo "COMMAND * STAR * could not be found, please configure it first!" >&2
        exit 1
    fi
}

## samtools
check_samtools(){
    if ! command -v samtools &> /dev/null
    then
        echo "COMMAND * samtools * could not be found, please configure it first!" >&2
        exit 1
    fi
}

## featurecounts
check_featureCounts(){
    if ! command -v featureCounts &> /dev/null
    then
        echo "COMMAND * featureCounts * could not be found, please configure it first!" >&2
        exit 1
    fi
}

## java
check_java(){
    if ! command -v java &> /dev/null
    then
        echo "COMMAND * java * could not be found, please configure it first!" >&2
        exit 1
    fi
}

## gatk
check_gatk(){
    if ! command -v gatk &> /dev/null
    then
        echo "COMMAND * gatk * could not be found, please configure it first!" >&2
        exit 1
    fi
}

## picard
check_picard(){
    if [[ ! -e $1 ]]; then 
        echo "Error: the picard file you give doesn't exist!" >&2
        exit 1
    fi
    if [[ "$1" != *"picard.jar"* ]]; then
        echo "Error: you need to give a valid picard.jar file!" >&2 
        exit 1
     fi 
    check_java
}

## bam-readcount
check_bamreadcount(){
    if ! command -v bam-readcount &> /dev/null
    then
        echo "COMMAND * bam-readcount * could not be found, please configure it first! Build Instructions can be found in https://github.com/genome/bam-readcount ; Meanwhile, you can also install it with conda, details can be found in https://anaconda.org/bioconda/bam-readcount " >&2
        exit 1
    fi
}

## python (pandas)
check_python(){
    if ! command -v python &> /dev/null
    then 
        echo "COMMAND * python * could not be found, please install python first!" >&2
        exit 1
    fi
    if ! python -c "import pandas" &>/dev/null
    then 
        echo "Can not find the Module 'pandas' for python. Please install pandas first!" >&2
        exit 1
    fi
}
subcommand=$1
case $subcommand in
    "") 
        echo "You need to specifiy which command you want to check!" >&2
        exit 1
        ;;
    *) 
        shift 
        check_${subcommand} $@
        if [ $? = 127 ]; then
            echo "Error: $subcommand is not a know subcommand." >&2
            exit 1
        fi
        ;;
esac
