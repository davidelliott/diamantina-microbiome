# Switch to bash:

exec bash
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1

mkdir diamantina
cd diamantina
wget -r --cut-dirs=2 -np -nH -R "index.html*" http://cgr.liv.ac.uk/illum/LIMS10796_978f1f02f3c84a18/Trimmed/
cd Trimmed
gzip -d -r ./

###
validate_mapping_file.py -m map.txt -o map_check/ -p -b

# remove unpaired reads
rm */*_R0_*.fastq

cd ..

multiple_join_paired_ends.py -i ./Trimmed/ -o joined

# remove unjoined reads (https://github.com/biocore/qiime/issues/1971)
rm ./joined.test/*/*.un*
# rm ./joined/*/*.un*

# rename the folders because they will later be used as sample identifiers.
# 1-dune1_TAAGGCGA-TAGATCGC_L001_R1_001/
rename 's/[0-9]*-([a-z]*)_*([a-z]*)([0-9]*).*/$1$2$3/' joined/*

# split libraries
multiple_split_libraries_fastq.py -i joined -o ./ --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

# pick OTUs
pick_open_reference_otus.py -i /home/manager/git/diamantina/seqs.fna -o /home/manager/git/diamantina/otus

# pick rep set
pick_rep_set.py -i otus/seqs_otus.txt -f seqs.fna -o rep_set.fna

# assign taxonomy - done by pick rep set
# make otu table - done by pick rep set

# Plot taxa summaries for all samples:
summarize_taxa_through_plots.py -o taxa_summary -i otus/otu_table_mc2_w_tax.biom -m map.txt
beta_diversity_through_plots.py -o beta -i otus/otu_table_mc2_w_tax.biom -m map.txt -t otus/rep_set.tre


pick_open_reference_otus.py -i ~/git/diamantina/Trimmed/16S-merged.fq -o ~/git/diamantina/cluster
