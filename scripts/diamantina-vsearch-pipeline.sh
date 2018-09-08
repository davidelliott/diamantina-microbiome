#!/bin/sh                                                                       

# This is an example of a pipeline using vsearch to process data in the         
# Mothur 16S rRNA MiSeq SOP tutorial dataset to perform initial paired-end      
# read merging, quality filtering, chimera removal and OTU clustering.      

# As starting point using guide from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline - accessed 2018-09-07    

THREADS=3
REF=../gold.fasta
PERL=$(which perl)
VSEARCH=/home/david/bin/vsearch-2.8.4-linux-x86_64/bin/vsearch

date

# echo
# echo Obtaining Gold reference database for chimera detection
#
# if [ ! -e gold.fasta ]; then
# 
#     if [ ! -e Silva.gold.bacteria.zip ]; then
#         wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
#     fi
# 
#     echo Decompressing and reformatting...
#     unzip -p Silva.gold.bacteria.zip silva.gold.align | \
#         sed -e "s/[.-]//g" > gold.fasta
# 
# fi
                                                         

# before the next step I took all the fastq sequences and put them in the "data" folder directly (originally they were each in folders per sample)

cd ../data

# re-naming the files (e.g. so 1-dune1_TAAGGCGA-TAGATCGC_L001_R1_001.fastq becomes "dune1_R1.fq"
rename 's/[0-9]*-([a-z]*)([_a-z]*)([0-9]*).*(R[0-9]).*/$1$2$3_$4.fq/' *

# the dune_flank samples were then re-named to duneflank (without underscore)

echo
echo Checking FASTQ format version for one file

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

$PERL ../map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL ../map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus.fasta \
    --otutabout all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" all.otus.fasta)

echo

### xxx

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel_sha1 \
    --centroids all.otus.fasta \
    --otutabout all.otutab.txt \
    --biomout all.otutab.biom


source activate qiime1
assign_taxonomy.py -i all.otus.fasta -r /home/david/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -t /home/david/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt

parallel_assign_taxonomy_rdp.py -r /home/david/gg_13_8_otus/rep_set/99_otus.fasta -t /home/david/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt -i all.otus.fasta -o rdp-assigned-tax-99/ -O 3

biom add-metadata -i vsearch-derep.ms2.biom -o vsearch-derep.ms2.rdp-tax-99.biom --observation-metadata-fp rdp-assigned-tax-99/rep-set.ms2_tax_assignments.txt --observation-header ID,taxonomy --sc-separated taxonomy

biom head -i all.otutab.biom
head all.otus.fasta

biom add-metadata -i all.otutab.biom -o all.otutab.w.taxa.biom --observation-metadata-fp uclust_assigned_taxonomy/all.otus_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy

biom head -i all.otutab.w.taxa.biom


#biom convert -i all.otutab.txt -o all.otutab.biom --to-hdf5 --table-type "OTU table"

#biom add-metadata -i all.otutab.biom -o all.otutab.w.taxa.biom --observation-metadata-fp uclust_assigned_taxonomy/all.otus_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy

#biom add-metadata -i all.otutab.biom -o all.otutab.w.taxa.biom --observation-metadata-fp uclust_assigned_taxonomy/all.otus_tax_assignments.txt --observation-header OTUID,taxonomy


echo Done

date


