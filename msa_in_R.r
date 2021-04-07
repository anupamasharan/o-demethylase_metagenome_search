#installing the software and packages 

#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install(version = "3.10") #make sure to check documentation page for the correct version of BiocManager compatible with the version of R you are running

#install some packages on which biomartr is dependent

#BiocManager::install("msa")
#BiocManager::install("Biostrings")

#load all libraries

library(msa)
library(Biostrings)

CP_db <- readAAStringSet("input_to_msa/all_CP_with_template.fasta")

MT1_db <- readAAStringSet("input_to_msa/all_MT1_with_template.fasta")

MT2_db <- readAAStringSet("input_to_msa/all_MT2_with_template.fasta")

#this file is the output generated with the python code plus the three template files that was used to get these sequences from RefSeq NCBI and we want to build hmm profile using this

CP_muscle_alignment <- msaMuscle(CP_db, maxiters = 10, verbose = TRUE)
MT1_muscle_alignment <- msaMuscle(MT1_db, maxiters = 10, verbose = TRUE)
MT2_muscle_alignment <- msaMuscle(MT2_db, maxiters = 10, verbose = TRUE)

#there are also optiopns to do ClustalW and ClustalOmega depending on user's preference, in addition to varying the parameters for the Muscle 

#Biostrings package has a lot of cool reading and writing functions for a number of common biological data file formats

write.phylip(CP_muscle_alignment, "msa_output/CP_alignment.phylip")
write.phylip(MT1_muscle_alignment, "msa_output/MT1_alignment.phylip")
write.phylip(MT2_muscle_alignment, "msa_output/MT2_alignment.phylip")

#this output can now be fed to the hmm utilities on the command line or a linux terminal.

