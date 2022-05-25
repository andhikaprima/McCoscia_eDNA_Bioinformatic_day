# Blasting MOTUs
# Original Script: Samuel Browett
# 5/22/2020
# 1. Set-up
# File 1. MOTU counts for each sample
# File 2. Sample details table
# 2. Set working directory and load libraries
# 3. Load Data
# 4. Blast
# 5. Filter Blast Search Results
# 6. Assign Taxonomy
# 7. Restrict Taxa Assignment
# 8. Create New Phyloseq Object
# 1. Set-up
# I reccommend you create a directory called 'BlastMOTUs_projname' (replace projname with whatever project name you want).
# 
# Inside this directory, create another directory called 'input_data'. This is the directory you want to keep your input files.
# 
# You will need two input .csv files:

# Modified Script: Andhika Prasetyo
# 11/17/2021
# 1. Create directory for output
# 2. phyloseq object to data frame
#

# 2. Set working directory and load libraries
# setwd("~/.../BlastMOTUs_projname/")

# Install required package
# install.packages("dplyr")

# install.packages("BiocManager")
# BiocManager::install(version = "3.9")
# BiocManager::install("phyloseq")



# install.packages("knitr")
#library(devtools)
#install_github("helixcn/seqRFLP")# install.packages("taxize")
# install.packages("stringr")

# Load the library
library(dplyr)
library(phyloseq)
library(knitr)
library(seqRFLP) # the seqRFLP package is required for turning a table into a fasta file
library(taxize)
library(stringr)

# Create new directory for output
dir.create("output")

# 3. Load Data
# Change to your MOTU count file names
MOTUcounts <- read.csv("./input_data/MOTU_counts.csv") # Read in Otu abundance table
otu_mat <- MOTUcounts %>% select (-c("total_reads", "sequence"))

# Each row is an individual MOTU The first column is the MOTU ID. The second column has the sequence.
tax_mat <- MOTUcounts %>% select (c("id", "sequence"))

# Change to your sample sheet file name
samples_df <- read.csv("./input_data/sample_details_sheet.csv")

# Phyloseq requires row names for each table
row.names(otu_mat) <- otu_mat$id        # Name rows by OTU ID
otu_mat <- otu_mat %>% select (-id)     # Since we have rownames, We no longer need the 'id' column

row.names(tax_mat) <- tax_mat$id        # Name rows by OTU ID
tax_mat <- tax_mat %>% select (-id)     # Since we have rownames, We no longer need the 'id' column

row.names(samples_df) <- samples_df$Sample   # Name rows by sample ID
samples_df <- samples_df %>% select (-Sample)  # Since we have rownames, We no longer need the 'sample' column

# Phyloseq requires the taxa and OTU tables to be in a Matrix format
otu_mat <- as.matrix(otu_mat)  # Change data frame format to Matrix
tax_mat <- as.matrix(tax_mat)

# Create your Phyloseq Object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE) # Formats otu_table for phyloseq
TAX = tax_table(tax_mat)                       # Formats taxa_table for phyloseq
samples = sample_data(samples_df)              # # Formats sample_sheet for phyloseq

phylo <- phyloseq(OTU, TAX, samples)

phylo # just entering the phyloseq object will give its details
save(phylo, file="./output/phyloseq_raw.RData")

# Tidy up workspace
# Removes all objects except 'phylo'
rm(list =setdiff(ls(), "phylo"))

# Remove singletons
# remove MOTUs represented by only a single read
ps <- prune_taxa(taxa_sums(phylo) > 1, phylo)
ps

## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 841 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 3 sample variables ]
## tax_table()   Taxonomy Table:    [ 841 taxa by 1 taxonomic ranks ]


# 4. Blast
# Take list of sequences and create a fasta file for blast search

# isolate the sequences from the phyloseq object
taxadf <- as.data.frame(as.matrix(tax_table(ps)[,"sequence"]))
# Keep only the MOTU identifier column and the sequence column - for blasting
motudf <- as.data.frame(cbind(MOTU = row.names(taxadf), sequence = as.character(taxadf[,"sequence"])))

head(motudf)

# convert table into a fasta file
seqRFLP::dataframe2fas(motudf, file="./output/MOTU_list_to_blast.fasta")

# Initiate blast search
# Take note of parameter settings

###############
## Next section is done outside of R
## we use blastn
## possibly activate from within R

query_file <- "output/MOTU_list_to_blast.fasta"
#blastn <- "blastn"
blastn <- "/extrastorage/home/andhika/SOFTWARE/ncbi-blast-2.9.0+/bin/blastn"

#check the system
#system2(blastn, 
#        c("-h"))

system2(blastn, 
        c("-db /extrastorage/home/refdbs/blastdb/nt/nt", # path to blastn database
          "-num_threads 23", # how many cores to use
          "-max_target_seqs 25", # maximum sequence targets found to be reported
          "-outfmt '6 std qlen ssciname staxid'", # what information to extract
          "-out ./output/motus.blasthits.txt", # name of output file
          "-qcov_hsp_perc 90", # minimum coverage %
          "-perc_identity 80", # minimum identity %
          "-query", query_file)) # the input/query file

# 5. Filter Blast Search Results
# So far, this is when the blast results are already available

path <- "./"
IDtable_name <- file.path(path,"./output/motus.blasthits.txt")

IDtable=read.csv(IDtable_name,
                 sep='\t',
                 header=F,
                 as.is=TRUE)

names(IDtable) <- c("qseqid","sseqid","pident",
                    "length","mismatch","gapopen",
                    "qstart","qend","sstart",
                    "send","evalue","bitscore",
                    "qlen","ssciname","staxid")

# Filter list of hits so it contains the tops hits for each MOTU (top hits defined as the best hit and ~0.5% down - set by 'margin')

margin <- 0.51 # 
new_IDtable <- IDtable[0,] # prepare filtered matchlist
ids <- names(table(IDtable$qseqid))
i=1
o=length(ids)
for (name in ids){
  print(paste0("progress: ", round(((i/o) * 100),0) ,"%")) # make a progressline
  test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
  max <- max(test$pident)
  test <- test[which(test$pident > (max-margin)),] # select all lines for a query
  #These lines can be included if analysing a taxonomic group with a lot of
  #"unassigned" sequences in GenBank, to exclude those from further evaluation.
  #test2 <- test[!grepl("uncultured eukaryote",
  #          test$truncated_ssciname,ignore.case = TRUE),]
  #if (nrow(test2) > 1) {test <- test2}
  #test <- test[!grepl("Environmental",
  #          test$truncated_ssciname,ignore.case = TRUE),]
  if (nrow(test) > 0 ) { test$string <- toString(names(table(test$ssciname))) }
  new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
  i=i+1
}

# Determine the most common taxid amongst hits
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


new_IDtable$majority_taxid <-  with(new_IDtable, ave(staxid, qseqid , FUN=Mode))
IDtable2 = new_IDtable[!duplicated(new_IDtable[c(1,17)]),]


# 6. Assign Taxonomy
# Assign taxonomy to MOTUs using the taxids This is where you will need the API key to increase speed

all_staxids <- names(table(IDtable2$staxid)) # get all taxids for table
all_classifications <- list() # prepare list for taxize output
o=length(all_staxids) # number of taxids

Start_from <- 1 # change if loop needs to be restarted due to time-out

for (cl in Start_from:o){ # the taxize command "classification" can be run on
  #the all_staxids vector in one line, but often there is
  #a timeout command, therefor this loop workaround.
  
  #make a progressline (indicating the index the loops needs to be
  #restarted from if it quits)
  print(paste0("processing: ", cl , " of ", o , " taxids"))
  all_classifications[cl] <- classification(all_staxids[cl], db = "ncbi")
}

# Reformat Taxonomy Details
#####################
output <- data.frame(staxid=character(),taxpath=character(),
                     stringsAsFactors=FALSE)
totalnames <- length(all_staxids)
for (curpart in seq(1:totalnames)){
  print(paste0("progress: ", round(((curpart/totalnames)
                                    * 100),0) ,"%")) # make a progressline
  currenttaxon <- all_classifications[curpart][[1]]
  if ( !is.na(currenttaxon)) {
    spec <- all_staxids[curpart]
    gen <- currenttaxon[which(currenttaxon$rank == "genus"),"name"]
    fam <- currenttaxon[which(currenttaxon$rank == "family"),"name"]
    ord <- currenttaxon[which(currenttaxon$rank == "order"),"name"]
    cla <- currenttaxon[which(currenttaxon$rank == "class"),"name"]
    phy <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"]
    kin <- currenttaxon[which(currenttaxon$rank == "kingdom"),"name"]
    spe <- currenttaxon[which(currenttaxon$rank == "species"),"name"]
    currentpath <- gsub(" ", "_",
                        #paste0("k__",kin,";p__",phy,";c__",cla,";o__",ord,";f__",fam,";g__",gen,";s__",spe))
                        paste0(kin,";",phy,";",cla,";",ord,";",fam,";",gen,";",spe))
    output[curpart,"staxid"] <-  spec # add row to the filtered IDtable
    output[curpart,"taxpath"] <-  currentpath # add row to the filtered IDtable
  }
}

# You may get the following warning:
# the condition has length > 1 and only the first element will be used
# This is expected and fine

# Save Taxonomy Details
# Next we save the file with all the taxonomic information
taxonomic_info <- merge(IDtable2,output,by = "staxid", all=TRUE)
tbname <- file.path(path,"./output/Table_otu_taxonomy1.txt")
{write.table(taxonomic_info, tbname, sep="\t",quote=FALSE, col.names = NA)}

# Another Tidy Workspace
rm(list =setdiff(ls(), c("phylo", "ps", "path")))


# Restructure the taxonomy information file: Part 2
#Split taxonomic string into levels for the OTU data
tab_name <- file.path(path,"./output/Table_otu_taxonomy1.txt")
otutaxonomy <- read.table(tab_name, sep="\t", header=TRUE, as.is=TRUE)
otulevels <- str_split_fixed(otutaxonomy$taxpath, ";", 7)
otulevels <- gsub(".__","",otulevels)
otulevels <- as.data.frame(otulevels)
names(otulevels) <- c("kingdom","phylum","class","order","family","genus",
                      "species")
otutaxlevels <- cbind(otutaxonomy,otulevels)
tab_name <- file.path(path,"./output/Table_otu_taxonomy_levels1.txt")
{write.table(otutaxlevels, tab_name, sep="\t",quote=FALSE, col.names = NA)}

# Everything has been taxonomically assigned
# All information is in the file "Table_otu_taxonomy_levels1.txt"


# Restructure the taxonomy information file: Part 3
tax.ass.tab <- read.csv("./output/Table_otu_taxonomy_levels1.txt", sep="\t", header = TRUE, as.is=TRUE)

taxonomy_table <- tax.ass.tab[,c("qseqid",  # the query sequence ID - i.e. the MOTU value/ID
                                 #"pident",  # percentage identity
                                 #"evalue",  # the e value
                                 "kingdom",  #
                                 "phylum",
                                 "class",
                                 "order",
                                 "family",
                                 "genus",
                                 "species",
                                 "pident",
                                 "evalue")]


newtaxatable <- as.data.frame(tax_table(ps))
newtaxatable$qseqid <- rownames(newtaxatable)
newtaxatable <- merge(taxonomy_table, newtaxatable, by= "qseqid", all=TRUE)
rownames(newtaxatable) <- newtaxatable$qseqid
newtaxatable <- newtaxatable %>% select (-qseqid)


# Clean up taxonomy table
tax.clean <- as.data.frame(newtaxatable)
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- "none" # if the values in this dataframe is NA, replace with blank i.e. ""

for (i in 1:nrow(tax.clean)){
  
  #Fill in missing taxonomy
  
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}


# 7. Restrict Taxa Assignment
# Restrict taxonomic level according to percentage identity
# 
# Any blast hits between 95% and 98% identity are restricted to genus level identification
# 
# Any blast hits between 93% and 95% identity are restricted to family level identification
# 
# Any blast hits between 90% and 93% identity are restricted to order level identification
# 
# Any blast hits less than 90% identity are restricted to class level identification
# 
# You can change these restrictions by altering the values in the following script

tax.clean$pident <- as.numeric(tax.clean$pident)
tax.clean[is.na(tax.clean)] <- 0

for (i in 1:nrow(tax.clean)){
  
  # identity percentage filter
  if (tax.clean[i,"pident"] == 0){
    tax.clean[i, 1:7] <- "none"
  } else if (tax.clean[i,"pident"] < 90){
    class <- paste("Class_", tax.clean[i,"class"], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,"pident"] < 93){
    order <- paste("Order_", tax.clean[i,"order"], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,"pident"] < 95){
    family <- paste("Family_", tax.clean[i,"family"], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,"pident"] < 98){
    genus <- paste("Genus_", tax.clean[i,"genus"], sep = "")
    tax.clean[i, "species"] <- genus
  } 
}


# 8. Create New Phyloseq Object
# Now to create a new phyloseq object with MOTUs wih taxonomy assigned using blast and NCBI. Singletons have already been removed.

#phylo.new <- phyloseq()
new.tax.mat <- as.matrix(tax.clean)
tax.assign.phylo <- phyloseq(tax_table(new.tax.mat), 
                             sample_data(sample_data(ps)), 
                             otu_table(otu_table(ps), 
                                       taxa_are_rows = TRUE))

tax.assign.phylo
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 100 taxa and 51 samples ]
## sample_data() Sample Data:       [ 51 samples by 3 sample variables ]
## tax_table()   Taxonomy Table:    [ 100 taxa by 10 taxonomic ranks ]


# Save final phyloseq object
save(tax.assign.phylo, file="./output/phyloseq_taxa_assigned_no_singletons.RData")


# Create CSV file

library("readxl")    # To read Excel files into R
#library("dplyr")     # To manipulate dataframes
#library("ggplot2")   # graphics

# install.packages("BiocManager")
# BiocManager::install("phyloseq")
#library("phyloseq") 

# Load Rdata
load("./output/phyloseq_taxa_assigned_no_singletons.RData")

# call the function to convert phyloseq object to data frame
source("phyloseq_to_df.R")


# ############################ Extract tax.assign.phylo data to study it
# 1st option
# Convert phyloseq object to data frame (for exporting)
# https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_to_df.html
# run phyloseq_to_df.R FIRST
direct_blast <- phyloseq_to_df(tax.assign.phylo, addtax = T, addtot = T, addmaxrank = F,
                             sorting = "abundance")

# Create cvs
write.csv(direct_blast,"./output/All_MOTUs_Direct_blast.csv", row.names = FALSE)


