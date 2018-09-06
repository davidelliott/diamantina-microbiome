# Switch to bash, set up QIIME environment:

exec bash
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1

# obtain raw data
mkdir data
cd data
wget -r --cut-dirs=2 -np -nH -R "index.html*" http://cgr.liv.ac.uk/illum/LIMS10796_978f1f02f3c84a18/Trimmed/
cd Trimmed
gzip -d -r ./

cd ..

# remove unpaired reads
rm Trimmed/*/*_R0_*.fastq

# join paired reads (needs fastq-join which can be installed using sudo apt install ea-utils)
multiple_join_paired_ends.py -i ./Trimmed/ -o joined

# remove unjoined reads (https://github.com/biocore/qiime/issues/1971)
rm ./joined/*/*.un*

# rename the folders because they will later be used as sample identifiers.
# 1-dune1_TAAGGCGA-TAGATCGC_L001_R1_001/
rename 's/[0-9]*-([a-z]*)_*([a-z]*)([0-9]*).*/$1$2$3/' joined/*

# check the mapping file
validate_mapping_file.py -m map.txt -o map_check/ -p -b

# split libraries
multiple_split_libraries_fastq.py -i joined -o ./ --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

### new section here 5 September 2018 ####

# chimera check using usearch 6.1 (see http://qiime.org/tutorials/chimera_checking.html)
identify_chimeric_seqs.py -i /home/david/git/diamantina-microbiome/data/seqs.fna -m usearch61 -o usearch_checked_chimeras/ -r /home/david/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta

# filter chimeras
filter_fasta.py -f /home/david/git/diamantina-microbiome/data/seqs.fna -o /home/david/git/diamantina-microbiome/data/seqs_chimeras_filtered.fna -s /home/david/git/diamantina-microbiome/usearch_checked_chimeras/chimeras.txt -n

# pick otus, assign taxonomy, make otu table
pick_open_reference_otus.py -i /home/david/git/diamantina-microbiome/data/seqs_chimeras_filtered.fna -o /home/david/git/diamantina-microbiome/data/usearch61_picked_otus -m usearch61 -a -O 3

# copy the biom file which we will be using in downstream analyses and providing as supplementary data:
cp data/usearch61_picked_otus/otu_table_mc2_w_tax_no_pynast_failures.biom data/otus/Elliott_Diamantina_Microbiome_Data_S1[otu-tax].biom

## note
# at wernerlab http://www.wernerlab.org/software/macqiime/citations it says "Another new chimera screening option is to use the UCHIME chimera detection during the OTU-picking step (pick_otus.py -m usearch61 ...), which means you never have to use identify_chimeric_seqs.py.". 
# at QIIME documentation http://qiime.org/scripts/pick_otus.html it says "Chimera checking with usearch 6.X is implemented in identify_chimeric_seqs.py. Chimera checking should be done first with usearch 6.X, and the filtered resulting fasta file can then be clustered. ".

# These 2 seem in contradiction. I have done it the way stated in QIIME docs for now but this should be checked by also doing it like this (without prior chimera check):
### pick_open_reference_otus.py -i /home/david/git/diamantina-microbiome/data/seqs.fna -o /home/david/git/diamantina-microbiome/data/otus -m usearch61




### end of new section - stuff below here was in the original version ###

# pick OTUs, assign taxonomy, make otu table
pick_open_reference_otus.py -i /home/david/git/diamantina-microbiome/data/seqs.fna -o /home/david/git/diamantina-microbiome/data/otus

# copy the biom file which we will be using in downstream analyses and providing as supplementary data:
cp data/otus/otu_table_mc2_w_tax.biom data/otus/Elliott_Diamantina_Microbiome_Data_S1[otu-tax].biom

# prepare a folder for results from R script
mkdir results
