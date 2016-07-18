1. merged

1.0 main directory

Contains the merged dataset. For each chunk:

- .info
A list of snps, with the following columns for each wave: info, freq, pass.

- .map

- .fam

- .gz
The dosage file

- .ngt
For each SNP, has value 1 if the SNP is genotyped in all waves, 0 otherwise.

- .ngt.matrix
A matrix with SNPs on one axis, waves on the other axis, and values 1 if the
SNP is genotyped in the wave, 0 otherwise.

------------

2. check

Contains checks of the original dataset and the merging.


2.0. Main directory

- check.txt
Contains the names of all chunks that failed one or more checks. (Should be empty.)

- covar.txt
Three columns: FamID, indID and wave of that individual.

- missing_stats.txt
Merging of all missing_stats.txt files in the 'chunks' directory

- missing.txt
Merging of all missing.txt files in the 'chunks' directory

- missing_with_indco.txt
Merging of all missing_with_indco.txt files in 'chunks' directory.

- snp_totals.txt
Lists the number of SNPs in each wave, for each chunk.


2.1. 'chunk' folder:

- chunks/[chunk]_check.txt
Contains the result of all checks run on the chunk

- chunks/[chunk]_missing.txt
A matrix of SNPs * waves indicating whether the SNPs are present in the waves.
Only SNPs that are not present in all waves are listed.

- chunks/[chunk]_missing_stats.txt
For numbers 1-[number of waves], indicates how many SNPs are present in exactly that number of waves.

- chunks/[chunk]_missing_with_indco.txt
Indicates how many individuals pass the inclusion criteria for each SNP.

