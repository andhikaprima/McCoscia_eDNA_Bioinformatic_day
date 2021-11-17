####Original tutorial
# Wolves’ diet based on DNA metabarcoding
# https://pythonhosted.org/OBITools/wolves.html
# Ref: Shehzad W, Riaz T, Nawaz MA, Miquel C, Poillot C, Shah SA, Pompanon F, Coissac E, Taberlet P (2012) Carnivore diet analysis based on next generation sequencing: application to the leopard cat (Prionailurus bengalensis) in Pakistan. Molecular Ecology, 21, 1951-1965.


wget --no-check-certificate https://pythonhosted.org/OBITools/_downloads/wolf_tutorial.zip
unzip wolf_tutorial.zip
cd wolf_tutorial

# run obitools
obitools

#assign SWARM
export PATH=$PATH:~/extrastorage/SOFTWARE/swarm/bin

# convert *.tab to *.csv in R
writeLines(gsub("\t", ",", readLines("input.tab")), "output.csv")


# 1. Paired-end alignment. Keep reads with quality > 30
illuminapairedend -r wolf_F.fastq wolf_R.fastq | obiannotate -S goodali:'"Good_WOLF" if score>30.00 else "Bad_WOLF"' | obisplit -t goodali

mkdir Output
cp Good_WOLF.fastq
cd Output

# 2. Demultiplex (need ngsfilter formatted .tsv file)
ngsfilter -t ../wolf_diet_ngsfilter.txt --fasta-output -u unidentified_WOLF.fasta Good_WOLF.fastq > WOLF.filtered.fasta

# 3. Filter sequences with lengths between 80 and 120 bp and with only 'ACGT'
obigrep -p 'seq_length>80' -p 'seq_length<120' -s '^[ACGT]+$' WOLF.filtered.fasta > WOLF.filtered_length.fasta

# 4. Group unique seqs
obiuniq -m sample WOLF.filtered_length.fasta >  WOLF.unique.fasta

# 5. Exchange the identifier to a short index _(obiannotate)_
obiannotate --seq-rank WOLF.unique.fasta | obiannotate --set-identifier '"'WOLF'_%09d" % seq_rank' > WOLF.fasta

# 6. Change formats to vsearch (in Owenia)
Rscript /home/naiara/SOFTWARE/scripts/owi_obifasta2vsearch -i WOLF.fasta -o WOLF.vsearch.fasta

# 6.1. Editing vsearch.fasta
vim WOLF.vsearch.fasta

# delete space between Motus name and ;
:%s/ ;size=/;size=/g
:wq

# 7. Desnoising. Run UCHIME de novo in VSEARCH for removing chimeras
vsearch --uchime_denovo WOLF.vsearch.fasta --sizeout --minh 0.90 --nonchimeras WOLF.nonchimeras.fasta --chimeras WOLF.chimeras.fasta --uchimeout WOLF.uchimeout.txt


#############################
## Clustering using SWARM
#############################

# 8. Clustering using SWARM
# -d 1
swarm -d 1 -z -t 10 -o WOLF_d1_SWARM3nc_output -s WOLF_d1_SWARM3nc_stats -w WOLF_d1_SWARM3nc_seeds.fasta WOLF.nonchimeras.fasta

# 9. Create the tab file
obitab -o WOLF.fasta >  WOLF.new.tab

# 10.1 Recount after SWARM
Rscript /home/naiara/SOFTWARE/scripts/owi_recount_swarm WOLF_d1_SWARM3nc_output WOLF.new.tab

# running locally
Rscript owi_recount_swarm WOLF_d1_SWARM3nc_output WOLF.new.tab

Reading swarm database...
Cluster database read including 1394 total clusters.
Calculating number of reads in each cluster
Kept only 63 clusters of size greater than or equal to 2 reads.
Reading tabulated database. This could take a while...
Database read including 3947 total different sequences and 4 samples.
Kept only 2614 sequences for calculations.
File WOLF_d1_SWARM3nc_output.counts.csv written

# 10.2. Select only non singleton MOTUs
# Open  YOURFILE_SWARM3nc_seeds.fasta in vim. Add a space in every header by changing “;size=” by “; size=”. 
vim  WOLF_d1_SWARM3nc_seeds.fasta
:%s/;size=/; size=/g
:wq

# 11. Select only non singletons
obigrep -p 'size>1' WOLF_d1_SWARM3nc_seeds.fasta > WOLF_d1_seeds_nonsingletons.fasta

# 12. Annotate with ecotag
ecotag -d ../embl_r117 -R ../db_v05_r117.fasta WOLF_d1_seeds_nonsingletons.fasta > WOLF_d1_ecotag.fasta 

47.792% of the alignments was cached

# 13. Add high level taxa 
Rscript /home/andhika/extrastorage/eDNA_scripts/owi_add_taxonomy WOLF_d1_ecotag.fasta

Reading ecotagged fasta file
Read 63 records
Output file WOLF_d1_ecotag.fasta.annotated.csv written with 63 sequences

# 14. Combine Ecotag and abundance files
Rscript /home/naiara/SOFTWARE/scripts/owi_combine -i WOLF_d1_ecotag.fasta.annotated.csv -a WOLF_d1_SWARM3nc_output.counts.csv -o WOLF_swarm-d1_MOTUs.csv

Reading ecotag database...
Ecotag database read including 63 total MOTUs.
Reading abundance database...
Abundances database read including 63 total MOTUs and 4 samples.
File WOLF_swarm-d1_MOTUs.csv written, including 63 MOTUs with 39786 total reads in 4 samples.
(63 non-singletons MOTUs).

#############################


#############################
## Clustering using SUMACLUST
#############################

# 8. Exchange the identifier to a short index for Sumaclust
obiannotate --seq-rank WOLF.unique.fasta | obiannotate --set-identifier '"'WOLF'_%09d" % seq_rank' > WOLF.new.fasta

# 9. Cluster at 99% with sumaclust
sumaclust -t 0.99 WOLF.new.fasta > WOLF.sumaclust99.fasta

===========================================================
 SUMACLUST version 1.0.31
 Alignment using SSE2 instructions.
===========================================================
Reading dataset...
3947 sequences
Indexing dataset... : Done
Sorting sequences by count...
Maximum ratio between the counts of two sequences to connect them: 1.000000
Clustering sequences when similarity >= 0.990000
Aligning and clustering... 
Done : 100 %       3093 clusters created.  

# 10. Create the tab file
obitab -o WOLF.sumaclust99.fasta >  WOLF.sumaclust99.tab

# 11. Recount abundances by sample
Rscript /home/naiara/SOFTWARE/scripts/owi_recount_sumaclust -i WOLF.sumaclust99.tab -o WOLF.sumaclust99.counts.csv

Reading obitab database...
Cluster database read including 3947 total sequences.
Database includes 3093 different clusters.
With 552 non-singleton clusters.
Calculating number of reads in non-singleton clusters
File WOLF.sumaclust99.counts.csv written, including 41119 reads in 3093 clusters.
(552 non-singletons clusters).

# 12. Get cluster centers
obigrep -p 'cluster_center' WOLF.sumaclust99.fasta >  WOLF.sumaclust99.centers.fasta

# 13. Taxonomic assignment using ecotag
ecotag -d ../embl_r117 -R ../db_v05_r117.fasta WOLF.sumaclust99.centers.fasta > WOLF.suma99.ecotag.fasta

96.060% of the alignments was cached


# 14. Add taxa above order level
Rscript /home/naiara/SOFTWARE/scripts/owi_add_taxonomy WOLF.suma99.ecotag.fasta

Reading ecotagged fasta file
Read 3093 records
Output file WOLF.suma99.ecotag.fasta.annotated.csv written with 3093 sequences

# 14. Combine ecotag and Sumaclust abundance files
Rscript /home/naiara/SOFTWARE/scripts/owi_combine -i WOLF.suma99.ecotag.fasta.annotated.csv -a WOLF.sumaclust99.counts.csv -o WOLF_suma99_MOTUs.csv

Reading ecotag database...
Ecotag database read including 3093 total MOTUs.
Reading abundance database...
Abundances database read including 3093 total MOTUs and 4 samples.
File WOLF_suma99_MOTUs.csv written, including 3093 MOTUs with 41119 total reads in 4 samples.
(552 non-singletons MOTUs).

#############################


## EXTRA. Collapse MOTUs (merging similar species)
## start column for sample (-s) to end column of sample (-e)

# Swarm
Rscript Rscript /home/andhika/extrastorage/eDNA_scripts/owi_collapse_ap -i WOLF_swarm-d1_MOTUs.csv -t 0.70 -s 28 -e 31 

Reading WOLF_swarm-d1_MOTUs.csv database
Database WOLF_swarm-d1_MOTUs.csv read including 63 sequences.
Database WOLF_swarm-d1_MOTUs.csv includes 63 sequences with species name,
Belonging to 17 unique species.
Now collapsing species with best_identity higher than 0.7.
55 sequences collapsed to 9 species.
Written WOLF_swarm-d1_MOTUs_collapsed.csv with 17 MOTUs.


# Sumaclust
Rscript Rscript /home/andhika/extrastorage/eDNA_scripts/owi_collapse_ap -i WOLF_suma99_MOTUs.csv -t 0.70 -s 17 -e 20

Reading WOLF_suma99_MOTUs.csv database
Database WOLF_suma99_MOTUs.csv read including 3093 sequences.
Database WOLF_suma99_MOTUs.csv includes 3093 sequences with species name,
Belonging to 29 unique species.
Now collapsing species with best_identity higher than 0.7.
3090 sequences collapsed to 26 species.
Written WOLF_suma99_MOTUs_collapsed.csv with 29 MOTUs.


