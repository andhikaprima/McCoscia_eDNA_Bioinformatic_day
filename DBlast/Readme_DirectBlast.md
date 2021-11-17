# Initial Set-up
# Create directory for Direct Blast then go to the directory
mkdir DBlast
cd Dblast
mkdir input_data

# Create 2 files and placed into "input_data" folder
# Work in excel then upload
# File 1. MOTU counts for each sample which consist of id, total_reads, individual sample's reads and sequence
# File 2. Sample details table which consistf of number, sample name, other info (eg. location, sample tipy, extraction, etc)

# Copy two file into directory Blasting_MOTUs_ap.r and phyloseq_to_df.R

# Run the Rscript
Rscript Blasting_MOTUs_ap.r

# wait the process
# your output file *.csv will be in folder "output" with and named as "All_MOTUs_Direct_blast.csv"
