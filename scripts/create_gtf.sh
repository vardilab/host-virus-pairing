#/bin/sh

FASTA=$1

module load bioawk

# 01. simplify header
awk '{ print $1 }' <  $FASTA > "01_"$FASTA

# 02. get gene sequence lengths
cat "01_"$FASTA | bioawk -c fastx '{print ">" $name ORS length($seq)}' | paste - - | tr -d ">" > "02_"$FASTA".length.txt"

# 03. create gene intervals file
echo -e "chrom\tstart\tend\tstrand\tgene_name" > "03_"$FASTA".gene_intervals.txt"; awk '{ split($0, a, " "); printf "%s\t%s\t%s\t%s\t%s\n", a[1], 1, a[2], 1, a[1] }' "02_"$FASTA".length.txt" >> "03_"$FASTA".gene_intervals.txt"

# 04. create a pseudo gtf file
sed '/chrom\tstart\tend/d;' "03_"$FASTA".gene_intervals.txt" | sed 's/lcl|//g; s/\([0-9]\+\t\)-1/\1-/; s/\([0-9]\+\t\)1\(\t[a-zA-Z]\+\)/\1+\2/' | awk '{print $1"\tprotein_coding\tgene\t"$2"\t"$3"\t.\t"$4"\t.\tgene_id \""$5"\";\n"$1"\tprotein_coding\texon\t"$2"\t"$3"\t.\t"$4"\t.\tgene_id \""$5"\"; transcript_id \""$5"\";"}' | sed '1s/^/#!genome-build EhVM1\n#!genome-version EhVM1\n/' > "04_"$FASTA".gene_intervals.gtf"

echo Done!
