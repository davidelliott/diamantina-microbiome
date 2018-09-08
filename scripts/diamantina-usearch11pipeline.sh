#!/bin/bash
# commands based on https://drive5.com/usearch/manual/ex_miseq.html
# memory limit was exceeded so trying vsearch instead (commands were adjusted as needed)
# some parts based on https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline

export usearch=/home/david/bin/usearch11.0.667_i86linux32
export vsearch=/home/david/bin/vsearch-2.8.4-linux-x86_64/bin/vsearch

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

rm -rf ../out
mkdir -p ../out
cd ../out

# before the next step I took all the fastq sequences and put them in the "data" folder directly (originally they were each in folders per sample)

cd ../data

# re-naming the files (e.g. so 1-dune1_TAAGGCGA-TAGATCGC_L001_R1_001.fastq becomes "dune1_R1.fq"
rename 's/[0-9]*-([a-z]*)([_a-z]*)([0-9]*).*(R[0-9]).*/$1$2$3_$4.fq/' *

cd ../out

# Merge paired reads
# Add sample name to read label (-relabel option)
# Pool samples together in raw.fq (Linux cat command)
for Sample in channel1 channel2 channel3 channel4 channel5 channel6 dune1 dune2 dune3 dune4 dune5 dune_flank1 dune_flank2 dune_flank3 dune_flank4 dune_flank5 filter1 filter2 filter3 filter4 filter5 nebka1 nebka2 nebka3 nebka4 nebka5 pan1 pan2 pan3 pan4 pan5
do
	$vsearch -fastq_mergepairs ../data/${Sample}*_R1.fq -reverse ../data/${Sample}*_R2.fq -fastqout $Sample.merged.fq -relabel $Sample.
	cat $Sample.merged.fq >> all.merged.fq
done


# 	cat $Sample.merged.fq >> all.merged.fq

# Strip primers (V4F is 19, V4R is 20)
# $usearch -fastx_truncate all.merged.fq -stripleft 19 -stripright 20 \
#   -fastqout stripped.fq
# - this is not needed because the sequencing provider did it already

# Calculate quality statistics
$vsearch --threads 3 --fastq_eestats all.merged.fq --output all.merged.fq.stats

# Quality filter
$vsearch -fastq_filter all.merged.fq -fastq_maxee 1.0 -fastaout filtered.fa -relabel Filt
# 9195657 sequences kept (of which 0 truncated), 177854 sequences discarded.


# Dereplicate at sample level and relabel with sample_n
$vsearch --threads 3 --derep_fulllength filtered.fa --strand plus --output $s.derep.fasta \
        --sizeout \
        --uc $s.derep.uc \
        --relabel $s. \
        --fasta_width 0



# Find unique read sequences and abundances
$usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

# Make 97% OTUs and filter chimeras
$usearch -cluster_otus uniques.fa -otus otus.fa -relabel Otu

# Denoise: predict biological sequences and filter chimeras
$usearch -unoise3 uniques.fa -zotus zotus.fa

##################################################
# Downstream analysis of OTU sequences & OTU table
# Can do this for both OTUs and ZOTUs, here do
# just OTUs to keep it simple.
##################################################

# Make OTU table
$usearch -otutab all.merged.fq -otus otus.fa -otutabout otutab_raw.txt

# Normalize to 5k reads / sample
$usearch -otutab_norm otutab_raw.txt -sample_size 5000 -output otutab.txt

# Alpha diversity
$usearch -alpha_div otutab.txt -output alpha.txt

# Make OTU tree
$usearch -cluster_agg otus.fa -treeout otus.tree

# Beta diversity
mkdir beta/
$usearch -beta_div otutab.txt -tree otus.tree -filename_prefix beta/

# Rarefaction
$usearch -alpha_div_rare otutab.txt -output rare.txt

# Predict taxonomy
$usearch -sintax otus.fa -db ../data/rdp_16s_v16.fa -strand both \
  -tabbedout sintax.txt -sintax_cutoff 0.8

# Taxonomy summary reports
$usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank g -output genus_summary.txt
$usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank p -output phylum_summary.txt

# Find OTUs that match mock sequences
$usearch -uparse_ref otus.fa -db ../data/mock_refseqs.fa -strand plus \
  -uparseout uparse_ref.txt -threads 1
