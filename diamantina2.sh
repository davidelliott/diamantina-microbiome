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

# pick OTUs, assign taxonomy, make otu table
pick_open_reference_otus.py -i /home/david/git/diamantina-microbiome/data/seqs.fna -o /home/david/git/diamantina-microbiome/data/otus

# copy the biom file which we will be using in downstream analyses and providing as supplementary data:
cp data/otus/otu_table_mc2_w_tax.biom data/otus/Elliott_Diamantina_Microbiome_Data_S1[otu-tax].biom

# prepare a folder for results from R script
mkdir results