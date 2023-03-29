#!/bin/sh

if [ "$1" == "-h" ]; then
    echo "Usage: Pull reads from individual cells, trim them, remove poly-A and assemble each cell individually Run.\nProvide the following ordered inputs: `basename $0` <a tab delimited table of cell barcodes with associated fastq file> \n Module requirements: gcc/9.2.0, spades/3.15.0, seqtk/1.2 TrimGalore/0.6.5 "
    exit 0
fi

BARCODES=$1 # List of cell barcodes. E.g. cells cells with high probability of infection

echo "Pulling sequence headers by cell barcode"
#FASTQB=$(basename  $FASTQ)

while IFS="" read -r barcode FASTQ || [ -n "$barcode" ]
    do
        FASTQB=$(basename  $FASTQ)

        BARCODE_NO_WHITESPACE="$(echo -e "${barcode}" | tr -d '[:space:]')"
        echo "$BARCODE_NO_WHITESPACE"
        mkdir -p $BARCODE_NO_WHITESPACE
        cd $BARCODE_NO_WHITESPACE
        for FASTQR1 in $FASTQ/*_R1_*.fastq.gz # Iterate fastq files from all four lanes
            do
            echo $FASTQR1
            zcat $FASTQR1 | grep -B1 $BARCODE_NO_WHITESPACE | grep @ | tr -d @ | awk '{print $1;}' >> $FASTQB.headers.temp # Greb sequence headers by cell barcode recorded in R1
            FASTQR2=${FASTQR1/R1/R2} # Switch to R2
            echo "Filtering fastq sequences"
            echo $FASTQR2
            f="$(basename $FASTQR2 .fastq.gz)"
            OUT=$f"_"$BARCODE_NO_WHITESPACE".filtered.fastq.gz"
            echo $OUT
            seqtk subseq $FASTQR2 $PWD/$FASTQB.headers.temp | gzip > $OUT # Get sequences stored in R4 by headers
            rm $FASTQB.headers.temp
            done

        for fastq in *"filtered"*".fastq.gz"
            do
            if [ -f "$fastq" ];then
                echo "$fastq"
                trim_galore --phred33 -j 8 --length 36 -q 5 --stringency 1 --gzip --fastqc -e 0.1 $fastq
                sleep 1
            fi
            done

        for fastq in *"trim"*".fq.gz"
            do
            if [ -f "$fastq" ];then
                echo "$fastq"
                trim_galore --polyA -j 1 --length 36 --gzip $fastq
                sleep 1
            fi
            done
        ## spades
        mkdir -p assembly
        for FQ in *L001*.filtered_trimmed_trimmed.fq.gz
            do
                NAME="$(basename $FQ .fq.gz)"
                fq1=$FQ
                fq2="$(echo ${FQ/L001/L002})" 
                fq3="$(echo ${FQ/L001/L003})"
                fq4="$(echo ${FQ/L001/L004})"
                
                ## run script
                /apps/RH7U2/gnu/SPAdes/3.15.0/bin/rnaspades.py --s 1 $fq1 --s 2 $fq2 --s 3 $fq3 --s 4 $fq4 -k 21 33 -t 8 -o assembly
                sleep 1
            done
        rm *.fq.gz
        rm *.filtered*
            
        cd ..
    done  < $BARCODES