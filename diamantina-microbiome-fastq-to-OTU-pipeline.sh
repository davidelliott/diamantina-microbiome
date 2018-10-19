#!/bin/bash

# these commands are not designed to be run at once as a script, although that might work it will depend on file paths etc.
# the guide from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline was used heavily to start with - accessed 2018-09-07    
# The perl script map.pl was also provided there
# fastqdump tool obtained from sra-toolkit - Compiled binaries of July 24, 2018, version 2.9.2 from https://github.com/ncbi/sra-tools/wiki/Downloads
# QIIME environment needs to be set up first, e.g. as follows:
# exec bash
# conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

fasterqdump=/home/david/bin/sratoolkit.2.9.2-ubuntu64/bin/fasterq-dump
THREADS=3
PERL=$(which perl)
VSEARCH=/home/david/bin/vsearch-2.8.4-linux-x86_64/bin/vsearch
DATA=/home/david/git/diamantina-microbiome/data
SCRIPTS=/home/david/git/diamantina-microbiome
RESULTS=/home/david/git/diamantina-microbiome/results
REF=$SCRIPTS/gold.fasta

date

echo
echo Obtaining Gold reference database for chimera detection
if [ ! -e gold.fasta ]; then

   if [ ! -e Silva.gold.bacteria.zip ]; then
       wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
   fi

   echo Decompressing and reformatting...
   unzip -p Silva.gold.bacteria.zip silva.gold.align | \
       sed -e "s/[.-]//g" > gold.fasta

fi

cd $DATA

# download fastq files from SRA
# all the files will be in the $DATA folder together 
# named with suffix _1 and _2 signifying the two paired reads for each sample
# echo
# echo 'Downloading sequence data in fastq from SRA'
# 
# $fasterqdump -e $THREADS -p -x SRX4216609 -o pan5
# $fasterqdump -e $THREADS -p -x SRX4216608 -o pan4
# $fasterqdump -e $THREADS -p -x SRX4216607 -o pan1
# $fasterqdump -e $THREADS -p -x SRX4216606 -o nebkha5
# $fasterqdump -e $THREADS -p -x SRX4216605 -o pan3
# $fasterqdump -e $THREADS -p -x SRX4216604 -o pan2
# $fasterqdump -e $THREADS -p -x SRX4216603 -o nebkha2
# $fasterqdump -e $THREADS -p -x SRX4216602 -o nebkha1
# $fasterqdump -e $THREADS -p -x SRX4216601 -o nebkha4
# $fasterqdump -e $THREADS -p -x SRX4216600 -o nebkha3
# $fasterqdump -e $THREADS -p -x SRX4216599 -o filter1
# $fasterqdump -e $THREADS -p -x SRX4216598 -o filter2
# $fasterqdump -e $THREADS -p -x SRX4216597 -o filter3
# $fasterqdump -e $THREADS -p -x SRX4216596 -o filter4
# $fasterqdump -e $THREADS -p -x SRX4216595 -o filter5
# $fasterqdump -e $THREADS -p -x SRX4216594 -o channel1
# $fasterqdump -e $THREADS -p -x SRX4216593 -o channel3
# $fasterqdump -e $THREADS -p -x SRX4216592 -o channel4
# $fasterqdump -e $THREADS -p -x SRX4216591 -o channel5
# $fasterqdump -e $THREADS -p -x SRX4216590 -o channel6
# $fasterqdump -e $THREADS -p -x SRX4216589 -o duneflank3
# $fasterqdump -e $THREADS -p -x SRX4216588 -o duneflank2
# $fasterqdump -e $THREADS -p -x SRX4216587 -o duneflank1
# $fasterqdump -e $THREADS -p -x SRX4216586 -o dune5
# $fasterqdump -e $THREADS -p -x SRX4216585 -o dune4
# $fasterqdump -e $THREADS -p -x SRX4216584 -o dune3
# $fasterqdump -e $THREADS -p -x SRX4216583 -o dune2
# $fasterqdump -e $THREADS -p -x SRX4216582 -o dune1
# $fasterqdump -e $THREADS -p -x SRX4216581 -o duneflank5
# $fasterqdump -e $THREADS -p -x SRX4216580 -o duneflank4

# using the SRA download appeared to work but the following sample processing would not work.
# It appears that the SRA change the files in some way
# So have to download direct from the sequencing provider instead:

cd $DATA
wget -r --cut-dirs=2 -np -nH -R "index.html*" http://cgr.liv.ac.uk/illum/LIMS10796_978f1f02f3c84a18/Trimmed/
cd Trimmed
gzip -d -r ./
cd ..

# the files end up in a folder called "Trimmed" and in individual folders per sample

# remove unpaired reads
rm Trimmed/*/*_R0_*.fastq

# move files out of sample named directories, put them all together into $DATA folder
mv Trimmed/**/*.* .
rm -R Trimmed

# re-name the files (e.g. so 1-dune1_TAAGGCGA-TAGATCGC_L001_R1_001.fastq becomes "dune1_R1.fq"
# etc.
rename 's/[0-9]*-([a-z]*)([_a-z]*)([0-9]*).*(R[0-9]).*/$1$2$3_$4.fq/' *

# re-name the dune_flank samples to duneflank (without underscore)
rename 's/dune_flank/duneflank/' *

# start processing with vsearch

echo
echo 'Checking FASTQ format version for one file'

$VSEARCH --threads $THREADS \
    --fastq_chars $(ls -1 *.fq | head -1)


# Process samples                                          

for f in *_R1.fq; do

    r=$(sed -e "s/_R1/_R2/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================

    $VSEARCH --threads $THREADS \
        --fastq_mergepairs $f \
        --reverse $r \
        --fastq_minovlen 200 \
        --fastq_maxdiffs 15 \
        --fastqout $s.merged.fastq \
        --fastq_eeout

    # Commands to demultiplex and remove tags and primers                       
    # using e.g. cutadapt may be added here.                                    

    echo
    echo Calculate quality statistics

    $VSEARCH --threads $THREADS \
        --fastq_eestats $s.merged.fastq \
        --output $s.stats
    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
        --fastq_filter $s.merged.fastq \
        --fastq_maxee 0.25 \
        --fastq_minlen 225 \
        --fastq_maxns 0 \
        --fastaout $s.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $s.filtered.fasta \
        --strand plus \
        --output $s.derep.fasta \
        --sizeout \
        --uc $s.derep.uc \
        --relabel $s. \
        --fasta_width 0

done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

rm -f all.derep.fasta all.nonchimeras.derep.fasta
cat *.derep.fasta > all.fasta

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --threads $THREADS \
    --derep_fulllength all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)

echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size all.derep.fasta \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.preclustered.uc \
    --centroids all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)

echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
    --uchime_denovo all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection

$VSEARCH --threads $THREADS \
    --uchime_ref all.denovo.nonchimeras.fasta \
    --db $REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

$PERL $SCRIPTS/map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL $SCRIPTS/map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with sha1 hash, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel_sha1 \
    --centroids $RESULTS/all.otus.fasta \
    --otutabout $RESULTS/all.otutab.txt \
    --biomout $RESULTS/all.otutab.biom

echo
echo Number of OTUs: $(grep -c "^>" $RESULTS/all.otus.fasta)
# Number of OTUs: 19729

echo

# use QIIME to assign taxonomy

source activate qiime1

## assign taxonomy with SILVA/uclust
assign_taxonomy.py -m uclust -r /home/david/silva132/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -t /home/david/silva132/taxonomy/16S_only/99/majority_taxonomy_7_levels.txt -i $RESULTS/all.otus.fasta -o $RESULTS/

source deactivate

echo Done

date

# the next steps are done in R

# the following files produced by this script will be used in R:
# - $RESULTS/all.otutab.biom
# - $RESULTS/all.otus_tax_assignments.txt

# In addition:
# - $RESULTS/all.otus.fasta should be retained as it contains all the OTU sequences


