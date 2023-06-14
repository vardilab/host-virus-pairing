  
# Zooming into the rare virosphere reveals the native host of giant viruses
Written by Amir Fromm  
amir.fromm@weizmann.ac.il  
Code repository and instruction for repeating the results of the manuscript "Zooming into the rare virosphere reveals the native host of giant viruses" (2023). 

### Dependencies:  

cutadapt/4.2  
cd-hit/4.6.6  
sortmerna/4.3.6  
bioawk  
bcl2fastq/2.20.0.422  
cellranger/5.0.0  
rsem/1.3.1  
perl  
bowtie2/2.3.3  
seqtk/1.2  
TrimGalore/0.6.5   
gcc/9.2.0  
spades/3.15.0  
BLAST+/2.11.0  
BBMap/38.90

### Scripts:

------------------
### A step-by step recreation of the results on the manuscript
#### Create a reference of viral marker genes  
input: NCLDV marker genes file (link)  
scripts:   
create_gtf.sh  
  
a. De-duplicate gene file  
```
cd-hit-est -c 0.9 -T 16 -i GVDB.markergenes.fna -o GVDB.markergenes.90.fna  
```
b. create a pseudo-GTF file  
```
create_gtf.sh GVDB.markergenes.try.90.fna  
```
c. create a reference file  
```
cellranger mkref --genome=$REFNAME --fasta=$FASTA --genes=$GTF --memgb=4  
```
#### Map raw fastq files to the reference  
```
cellranger count --id=$ID --transcriptome=$REF --fastqs=$dir --project=$ID --chemistry=SC3Pv3  
```
#### Generate raw UMI counts from 10X single-cell data in a pickle format 
```
script: 00.00.raw_UMI_counts_10X_data.py  
  
usage: 00.00.raw_UMI_counts_10X_data.py [-h] --data_dir DATA_DIR  
                                        [--file_type FILE_TYPE]  

optional arguments:  
  -h, --help            show this help message and exit  
  --data_dir DATA_DIR   Path to direcory with 10X outputs: matrix.mtx.gz,  
                        features.tsv.gz, barcodes.tsv.gz)  
  --file_type FILE_TYPE  
                        Note if the input count matrix is in format mtx, csv,  
                        or tsv (default: mtx)  
  
example: python 00.00.raw_UMI_counts_10X_data.py --data_dir ./  

```
output: data_raw.pickle.gz  
  
#### Combine pickle files  
combine multiple UMI tables of different samples in a pickle format. 
```
script: combine_pickle.py  
usage:   
    combine_pickle.py --data_dir DATA_DIR  
                     --list_folders LIST   
                     --output_folder OUT   
                     --raw  
  
arguments:                      
    --data_dir  DATA_DIR    Path to directory containing the folders of the 10X outputs  
    --list_folders  LIST    List of folders containing the 10X files (in pickle format), comma saparated; write the name of the parent folder, not the '/outs/filtered_feature_bc_matrix' path  
    --output_folder OUT     Name of output folder for combined file  
    [--raw|--not_raw]       Specify if raw data or pre-processed  
  
                            Note, at this stage there is no metadata file therefore only --raw will be accepted  
  
example: combine_pickle.py --data_dir ./ --list_folders "B7T16,B7T17,B7T18" --output_folder Combined_B7R --raw  
The name of each cell will be the barcode+the name of the folder containing the 10X outputs (for example, ATGCA.B7T17).  
```

#### A subset of cells with high expression of viral transcripts is selected.  
```
usage:  
    choose_cells.py --data_dir DATA_DIR  
                        --sample_table TABLE  
                        --output_folder OUT  
                        --sum [int]  
                        --count [int]   
                        --exp [int]  
arguments:  
    --data_dir DATA_DIR     Path to directory containing the combined raw data  
    --sample_table TABLE    A tab-delimited table of sample names and the complete path to its raw fastq file  
    --output_folder OUT     Name of folder for output file  
    --sum [int]             Threshold value: The sum of UMIs in a cell should be equal or higher than this number, default:1  
    --count [int]           Threshold value: The total number of genes expressed in a cell should be equal or higher than this number, default:1  
    --exp [int]             Threshold value: The most highly expressed gene should have a number of UMIs equal or higher than this number, default:1  
  
example: choose_cells.py --data_dir /. --sample_table sample_table_fastq.tsv --output_folder /. --sum 10 --count 2 --exp 2   

```
#### Single-cell reads are extracted from each selected cell, trimmed, and assembled  

```
usage:  
    assemble_cells.sh BARCODES  
  
arguments:  
BARCODES    a tab delimited file of cell barcodes and the raw fastq file of their sample (output of choose_cells.py)  
```  
output: multiple folders, one for each cell, containing single-cell reads, and a folder of the assembly. 

#### Transcripts and read files are named after the cell barcode  
This needs to be executed in the parent folder containing per-cell folders.  
```
for some_path in */  
do  
    CELL=$(basename $some_path)  
    for file in $some_path/*_L001_R2_001_*.filtered.fastq.gz  
    do  
        SUBSTRING="${file##*/}"  
        ID="$(awk -F '_L' '{print $1}' <<< "$SUBSTRING" )"  
    done  
    cat $some_path/*.filtered_trimmed_trimmed.fq.gz > $some_path/$CELL.$ID.combined.fq.gz  
    sed "s/>.*/&_$CELL.$ID/" $some_path/assembly/transcripts.fasta > $some_path/assembly/transcripts.edit.fasta  
done  
```  
output: Trimmed single cell reads from each cell, asssembled contigs from each cell, named after the cell barcode.  
  
#### Transcript files from each cell are concatanated   
```  
cat parent_folder/*/assembly/transcripts.edit.fasta > all_cells.transcripts.edit.fasta  
```
#### Transcripts are filtered using sortmerna against the pr2 database  
NOTE: We de-duplicated the pr2 database using cd-hit at 99% identity.   
It is not obligatory but it accelerates the process.  
```    
sortmerna --ref  PR2.99.fasta \  
--reads all_cells.transcripts.edit.fasta \  
--aligned sortmerna/all_cells.transcripts.edit.filtered \  
-fastx sortmerna/all_cells.transcripts.edit.filtered  
```  
output: sortmerna/all_cells.transcripts.edit.filtered.fa  
  
#### Reads are aligned to PR2 and metaPR2 database using blast  
```    
blastn -outfmt 6 -evalue 1e-10 -perc_identity 99 -query all_cells.transcripts.edit.filtered.fa -subject 099.pr2.fasta -out transcripts.PR2.tsv  
blastn -outfmt 6 -evalue 1e-10 -perc_identity 99 -query all_cells.transcripts.edit.filtered.fa -subject 099.metapr2.fasta -out transcripts.metaPR2.tsv  
``` 
#### Blast results are summarized  
Results of 18s rRNA homology against the PR2 and metaPR2 databases are summarized.  
Pick only the best hit with identity of >= 99 percent and an alignment length of >100 bp  
```  
usage:  
    01.summarize_blast.py --data_dir DATA_DIR   
                        --database DB  
  
argument:  
    --data_dir DATA_DIR Path to directory containing the blast results file  
    --database DB       Name of database used (either metaPR2 or PR2) name has to be: all_cells.transcripts.edit. + DB + .tsv  
```  
output: a summary file for each result file selected.  
  
#### Find homology to virus by mapping raw reads to virus marker genes  
(frank)  
  
#### Link results of 18S rRNA search and virus homology search to pair hosts and viruses  
The plot is produced in a jupyter notebook. explanations and input requirements are in the notebook:  
sankey_wrapper.ipynb  
  
#### Creating a reference database from identified single-cells  
  
Assembled transcripts de-duplicated using cd-hit  
```    
cd-hit-est -c 0.95 -T 32 -i all_cells.transcripts.edit.fasta -o all_cells.transcripts.95.fasta  
```  
#### Low-complexity transcripts are removed from the transcript file  

input: deduplicated transcripts from all cells (cells.transcripts.95.fasta)  
output: deduplicated transcripts, with low-complexity transcripts removed.  
```  
bbduk.sh in=cells.transcripts.95.fasta out=cells.transcripts.95.bbduk.fasta outm=low_complexity.fq entropy=0.7  
```
NOTE: We also manually removed a long repetitive sequence from the output file.  

#### Transcripts are selected from the 72 linked cells in order to create a reference database   
  
inputs: deduplicated transcripts, with low-complexity transcripts removed from all cells (cells.transcripts.95.bbduk.fasta); text file with barcodes from the 72 cells (filtered_cells_list.txt)  
output: deduplicated transcripts, with low-complexity transcripts removed from the 72 identified cells  
```
filterbyname.sh in=cells.transcripts.95.bbduk.fasta out=cells.filtered.fasta include=t names=filtered_cells_list.txt substring  
```
#### A new reference is curated from selected cells, E. huxleyi and E .huxleyi virus  
cat transcriptome_ehuxleyi.EhVM1.fasta cells.filtered.fasta > combined_cells_ehux.fasta  
  
#### The new transcriptome is de-duplicated using cd-hit  
```
cd-hit-est -c 0.95 -T 16 -i combined_cells_ehux.fasta -o combined_cells_ehux.95.fasta  
```
#### A reference database is built from from the combined transcriptome  
  
input: combined database file  
scripts:   
create_gtf.sh  
  
a. create a pseudo-GTF file  
```  
create_gtf.sh combined_cells_ehux.95.fasta  
```  
b. create a reference file  
```
cellranger mkref --genome=$REFNAME --fasta=$FASTA --genes=$GTF --memgb=4  
```  
#### Map raw fastq files to the reference  
```
cellranger count --id=$ID --transcriptome=$REF --fastqs=$dir --project=$ID --chemistry=SC3Pv3  
```  
#### Generate raw UMI counts from 10X single-cell data into a pickle format  
See instructions in  section #6  
  
#### Preprocessing UMI counts: filter, normalize, scale  
```  
usage: 00.01.filter_normalize_scale_single_cell_data.py [-h] --data_dir DATA_DIR  
                                                        [--file_type FILE_TYPE]  
                                                        [--libsize_perc LIBSIZE_PERC]  
                                                        [--mit_cutoff MIT_CUTOFF]  
                                                        [--min_cells MIN_CELLS]  
  
optional arguments:  
  -h, --help            show this help message and exit  
  --data_dir DATA_DIR   path to direcory with 10X outputs: matrix.mtx.gz,  
                        features.tsv.gz, barcodes.tsv.gz)  
  --file_type FILE_TYPE  
                        Note if the input count matrix is in format mtx, csv,  
                        or tsv (default: mtx)  
  --libsize_perc LIBSIZE_PERC  
                        int or tuple of ints, above or below which to retain a  
                        cell. Must be an integer between 0 and 100 (default:  
                        (5,99))  
  --mit_cutoff MIT_CUTOFF  
                        Remove cells with total expression of a mitochondrial  
                        genes above or below a threshold (default: 200)  
  --min_cells MIN_CELLS  
                        Filter all genes with negligible counts in all but a  
                        few cells (default: 5)  
  
example: python 00.01.filter_normalize_scale_single_cell_data.py --data_dir ./ --mit_cutoff -1 --min_cells 2 --file_type pickle  
```

#### Combine preprocessed UMI pickle tables from multiple samples  
See the section "Combine pickle files" for instructions.  
Use the --not_raw option.  
  
#### Get a list of all the identified cells (after filtering)   
See the section "A subset of cells with high expression of viral transcripts is selected." for instructions.  
Use low thresholds (1) to capture all cells.  
example: choose_cells.py --data_dir combined_data/. --sample_table sample_table_fastq.tsv --output_folder /. --sum 1 --count 1 --exp 1   
  
#### Single-cell reads are extracted from each selected cell, trimmed, and aligned  
See the section "Single-cell reads are extracted from each selected cell, trimmed, and assembled" for example.  
Use the output of the previous section.  
  
#### Transcripts and read files are named after the cell barcode  
See the section "Transcripts and read files are named after the cell barcode" for instructions.  
  
#### Transcript files from each cell are concatenated   
See the section "Transcript files from each cell are concatanated" for instructions.  
  
#### Transcripts are filtered using sortmerna against the pr2 database  
See the section "Transcripts are filtered using sortmerna against the pr2 database" for instructions.  
  
#### Reads are aligned to PR2 and metaPR2 database using blast  
See the section "Reads are aligned to PR2 and metaPR2 database using blast" for instructions.  
  
#### Blast results are summarized  
See the section "Blast results are summarized" for instructions.  
  
#### Viral transcripts are determent from transcript reference  
(frank)  
  
#### Combine raw UMI pickle tables from multiple samples  
See the section "Combine pickle files" for instructions.  
Use the --raw option.  
  
#### Implement dimensionality reduction on processed UMI tables  
```
usage: 00.02.dimentionality_reduction_single_cell_data.py [-h] --data_dir DATA_DIR  
                                                          [--pca_components PCA_COMPONENTS]  
                                                          [--tsne_perplexity TSNE_PERPLEXITY]  
                                                          [--data_file DATA]   
                                                          [--pca_data PCA]   
                                                          [--metadata_file METADATA]   
                                                          [--metadata_output METADATA_OUT]   
                                                          [--data_pca_output PCA_OUT]   
                                                          [--min_dist MINDIST]   
                                                          [--spread SPREAD]  
  
optional arguments:  
  -h, --help            show this help message and exit  
  --data_dir    DATA_DIR   
                        path to direcory with preprocessed single-cell data,  
                        including file names: data.pickle.gz and  
                        metadata.pickle.gz  
  --pca_components  PCA_COMPONENTS  
                        Number of components to calculate for the PCA  
                        (default: 50)  
  --tsne_perplexity TSNE_PERPLEXITY  
                        float. The perplexity is related to the number of  
                        nearest neighbors that is used in other manifold  
                        learning algorithms. Larger datasets usually require a  
                        larger perplexity. Consider selecting a value between  
                        5 and 50. Different values can result in significantly  
                        different results (default: 30)  
  --data_file   DATA     
                    name of data file as appear in data_dir (default: data.pickle.gz  
  --pca_data	PCA	     
                    name of the pca data to use for dimensionality reduction.   
                    If not mentioned, it will produce a new dataset (default = -1)  
  --metadata_file	METADATA  
                    name of metadata file as appear in data_dir   
                    (default: metadata.pickle.gz  
  --metadata_output	METADATA_OUT  
                    name of output metadata file   
                    (default: metadata_dimentionality_reduction.pickle.gz)  
  --data_pca_output	PCA_OUT  
                    name of output pca data file (default: data_pca.pickle.gz)  
  --min_dist    MINDIST  
                    value of the min_dist variable in UMAP (default: 0.1)  
  --spread  SPREAD  
                    value of the spread variable in UMAP  (default: 1)  
  
  
example: python 00.02.dimentionality_reduction_single_cell_data.py --data_dir ./ --min_dist 0.15 --spread 0.75 --pca_components 10   
```  
#### Annotating cells, creating UMAP projection plots oh host-virus coexpression and extracting the barcodes of infected Katablepharidacea cells.   
This is done in a jupyter notebook. explanations and input requirements are in the notebook:  
Coexpression_wrapper.ipynb  
  
#### Single-cell reads are extracted from each selected kata cell, trimmed, and poly-A removed  
```  
pull_trim_clean.sh  
   
  
usage:  
    pull_trim_clean.sh BARCODES FASTQ  
  
arguments:  
BARCODES    A file containing the cell barcodes   
FASTQ       Path to raw fastq dir of a specific sample.

```  
#### Concatenate trimmed fastq files from all cell barcodes and assemble reads altogether  

 ```
cat */*trim*trim*.fq.gz All_cells.fq.gz  
  
rnaspades.py --s 1 All_kata_cells.fq.gz -k 21 -t 8 -o $OUTDIR  
```
  
#### Transcripts are aligned to host 18s rRNA databases using blast  
Note that the identity is lower because we aim for a lower taxonomic level, therefore we except lower identity  
```  
blastn -outfmt 6 -evalue 1e-10 -perc_identity 90 -query transcripts.fasta -subject 099.metapr2.fasta -out transcripts_kata.metaPR2.tsv  
blastn -outfmt 6 -evalue 1e-10 -perc_identity 90 -query transcripts.fasta -subject 099.pr2.fasta -out transcripts_kata.PR2.tsv  
```  
#### Transcripts are aligned to host viral marker gene databases using blast  
```  
blastn -outfmt 6 -evalue 1e-10 -query transcripts.fasta -subject GVDB.markergenes.fna -out transcripts.GVDB.marker.tsv  
```  
#### After we find the virus with high certainity, we align the reads to to virus transcriptome to measure gene expression  
Prepare rsem reference:  
```
rsem-prepare-reference \  
--bowtie2 GVMAG-M-3300020187-27.fltr.fna.genes.fna \  
GVMAG-M-3300020187-27.fltr.fna.genes.fna  
```
Calculate gene expression
```  
rsem-calculate-expression -p 10 --bowtie2 --fragment-length-mean 58 All_cells.fq.gz GVMAG-M-3300020187-27.fltr.fna.genes.fna Kata_mapping  
``` 
#### Annotate the virus genes  
(Frank)
  
#### Find relative abundance of different taxa in the mesocosm experiment and the relative abundance of Katablepharidles  
This is done in a jupyter notebook:  
relative_abundance_all.ipynb  
The input data is the ASV analysis from Vincent et al. 2023 (Nature microbiology).  
Explanations can be found in the notebook  

