# <i>O</i>-demethylase_metagenome_search
Instructional tool to help retrieve co-localised o-demethylase operons from metagenomic data

## Background
Corrinoid dependent <i>O</i>-demethylase enzyme systems (referred to henceforth as o-demethylases) typically consist of multiple component gene products that catalyse the important chemical modification of methylation and demethylation of aromatic ether compounds. They present a compelling biotransformation alternative to the more common S-adenosyl-methionine (SAM) dependent methyl transfer systems which catalyse the same process (as depicted in figure 1) irreversibly and therefore cannot be scaled up very well making the former systems more potent biocatalysts. Aromatic methyl ethers are degradation products of lignin, a very abundant aromatic renewable biopolymer. Demethylation of these compounds therefore is very important for development of the global bioeconomy [1].

![alt text](https://github.com/anupamasharan/o-demethylase_metagenome_search/blob/main/reaction_mechanism.png?raw=true)

Anaerobic acteogens and organohalide-respiring organisms encode this complicated 4-component system consisting of two methyltransferases [2] (referred henceforth as MT1 and MT2) a cobalamin dependent corrinoid protein(CP) which is the shuttle and an activating enzyme (AE) that keeps regenerating CP in the catalytic cycle. The first three components occur together, co-localised either in triplets or in pairs of one CP and MT1 or MT2. The AE occurs elsewhere in the genome (figure 2)[3]. The co-localised gene products have been heterologously produced and expressed over several years, but AE remains elusive. However, the information about these enzyme systems is not organized. As a researcher working with testing anaerobic metagenomes on biodegradation of aromatics, I had to face a lot of trouble getting streamlined information about biochemically validated components to look for similar enzymes in-silico. There are different names for the enzymes as present on publicly searchable databases and three different proteins make it further complicated. Given the promise of these enzymes in the area of lignin-aromatic biotransformation, there should be a streamlined method to search for gene candidates encoding these enzymes in genomes or metagenomes of interest. There are unique features about the enzyme system that can be manipulated to selectively search for robust candidates for heterologous production and assay of the <i>o</i>-demethylase activity.

![alt text](https://github.com/anupamasharan/o-demethylase_metagenome_search/blob/main/gene_operon_organization.png?raw=true)

## The pipeline 

This is a tutorial-style guide to enable users with minimal bioinformatics experience to obtain co-localised o-demethylase sequences from public or custom metagenomic datasets, following creation of a hidden markov model (hmm) generated reference dataset. The sequences are retrieved using one of the following two co-localisation criteria:

1. Strict numerical co-localisation - suitable for databases where contig/scaffold information is not present in protein file headers (such as NCBI)
(use script NCBI_Protein_hit_co_localised_sequence_headers_retrieval.py)
2. (Recommended) Contig based co-localisation - wherever data linking genes to contigs/scaffols are available this script should be used as it is more representative of how these systems would be naturally found in envrionmental genomes. (JGI_IMG_Protein_hit_co_localised_sequence_headers_retrieval.py)
Through this pipeline I suggest an approach to obtain protein targets suitable for direct biochemical characterization and assay development for functional o-demethylation. However, if the user's objective is to do explorative or comparitive metagenomeic studies for o-demethylation function, all sequence-search hits should be analysed. Also, all parts of the workflow are tunable to the user preference and level of expertise.

###### Note : The scripts and code used in the pipeline are provided as separate file uploads. The sample file formats needed for input in the different steps of the pipeline and output files generated are provided in a winrar zipped archive "sample_input_output_allsteps". The user should save all folders within the zipped archive to their working directory to test the pipeline.

### Part 1 Creating a reference hmm profile

#### 1.1: Retrieving sequences

###### Note : You can make your own reference database by linking to any online repositories similar to NCBI.  

In this part we will see how the reference database supplied on the project page was generated by mining sequences from <a href = "https://www.ncbi.nlm.nih.gov/refseq/"> NCBI RefSeq Database </a> using 3 previously characterised sequence templates. This is done using Biopython tools and tidied in python, supplied as a jupyter lab notebook <b>(ref_db_creation.ipynb)</b>.

#### 1.2: Creating co-localised reference sequence database

The CP, MT1 and MT2 sequences obtained in part 1.1 can now be mined for the co-localised operons to create the reference database. Use the python script attached <b>(python_header_function.py)</b> and modify as per your header files. Sample csv files are attached as reference generated from step 1.1 in the <b>python header function input and output generated in the output folder</b>. Following extraction of co-localised headers, the csv output file can be fed as input to the python script <b>"bio_efetch.py"</b> to get the fasta sequences of the IDs from the header files. Example output from this script is provided in the <b>"efetch_output"</b> folder. This script uses the "efetch" module from the <a href = "https://biopython.org/docs/1.75/api/Bio.Entrez.html"> Bio.Entrez </a> package. Outside biopython this can also be done outside python using NCBI's <a href = "https://www.ncbi.nlm.nih.gov/sites/batchentrez"> batch enterez tool </a> to obtain all the fasta sequences.

###### Note : The NCBI_Protein_hit_co_localised_sequence_headers_retrieval.py script only searches for genes/operons that are immediately in the vicinity of each other (positions differ by +/- 1). However, demethylase enzyme systems not specific to aryl-methyl ethers can have different orientations/organizations as well. Read more in this review by Schilhabel et al., 2009. <a href = "https://jb.asm.org/content/191/2/588">[3] </a>. The script can be easily modified to accommodate combinations of numerical position based organization (lines 59, 67). The bio_efetch script will need the input headers formatted as per the example files supplied here to produce the desired output. Also, the output includes pairs of CP with MT1, MT2 in addition to triplets where CP, MT1 and MT2 are all co-localised. For the most stringent searches, such as obtaining candidates for heterologous production, only the triplets should be considered.

#### 1.3: Build hidden markov model (using hmmer)

From this step onwards we will be working in mostly R and a little in Linux (for hidden markov model generation). The bioconductor software in R has packages that allow for constructing and working with multiple sequence alignments (msa) directly using well-known algorithms such as ClustalW and MUSCLE. To keep use biopython, you would need to have a msa tool installed in command line and then run a <a href = "http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec92"> wrapper </a> from within Biopython. The user should keep in mind that running these algorithms from within R is computationally expensive. For bigger fasta file inputs, consider working on the command line directly. Working in R has another advantage of direct interfacing with some visualization tools.

##### 1.3.1 Multiple Sequence Alignment using msa package in R

There is a <a href = "https://bioconductor.org/packages/release/bioc/html/msa.html"> msa package </a> that can be installed with BioConductor which let’s us do msa very easily in R. Make sure <a href = "https://bioconductor.org/install/"> BioConductor </a> is installed beforehand compatible with your version of RStudio. Refer to the R script <b>msa_in_R.r</b> and generate the phylip file compatible with hmmer command line tools. The algorithm used here is MUSCLE with basic parameters. For more details on the two packages used in the script, refer to the online manuals, <a href = "https://bioconductor.org/packages/release/bioc/manuals/msa/man/msa.pdf"> msa </a> and <a href= "https://bioconductor.org/packages/release/bioc/vignettes/Biostrings/inst/doc/MultipleAlignments.pdf"> Biostrings </a>. Sample input and output files for this step has been provided in msa_input and msa_output folders.

###### Note : The hmm database will contain all the sequences we obtained from RefSeq database plus the 3 template sequences we used to retrieve those sequences as we want the hmm profile to contain all these sequences. The input fasta to the msa script is provided in the folder with the same name. The output from msa step is the input to the hmm building step.

##### 1.3.3 hmmer tools using Linux command line/terminal

If you already have a Linux or OS system, you can install the hmmmer binaries in your path and directly access them through R. With a windows system it is much more difficult and needs cgywin installation. Also, a bigger motivation to do this outside R was because the hmmsearch function within R cannot search using custom hmm-profiles very easily without interfacing with some other functions/packages, is not very intuitive. 

In python on the other hand, it seems possible to run hmm commands using the <a href = "https://pyhmmer.readthedocs.io/en/stable/index.html"> pyhmmer </a> package. However, when installing in a windows system, it might give some incompatibility with Microsoft Visual C++ and might require installing of Microsoft visual studio. 

On terminal it can be done very easily. The following command was used (assuming <a href = "http://hmmer.org/"> hmmer is </a> already installed following instructions on website).

Building hmm profile

$ hmmbuild o/p_filename.hmm input_alignment.phylip

Sample hmms for CP, MT1 and MT2 sequences are provided in the <b>hmm_output</b> folder.

### Part 2: Retrieve o-demethylase sequences from custom metagenomic database

This part of the tool demonstrates how to download metagenomic datasets from online databases if the user wishes to do an exploratory data analysis from publicly available dataset. If the user has their own dataset, you can skip to step 2.2. 

#### 2.1: Download publicly available metagenomic datasets

##### 2.1.1 NCBI

This is demonstrated using <a href = "https://docs.ropensci.org/biomartr/"> biomartr </a> package in R. This package interfaces with the NCBI portal through R. One shortcoming of this package is that it downloads whole datasets with no way of selecting specific metagenomic datasets and/or files within those datasets. The <b>metegenome_dataset_download</b> R markdown file provided here also includes other packages to overcome these shortcomings and allow the user to be specific with dataset selection. This file can be easily modified to get any metagenome dataset of interest on NCBI and can also be converted into an executable script. The file downloaded for the anaerobic digester metagenome can be found in the <b>metagenomic_data</b> folder.

###### Note: Usually very big datasets can be more efficiently downloaded using online NCBI FTP site or the JGI-IMG download portal. The R script included in this tutorial is to enable easy interfacing with the other steps, implemented in R.

##### 2.1.2 JGI IMG/M

These datasets can be found at <a href = "https://img.jgi.doe.gov/"> the website home page </a> using either accession IDs or keyword search. Please note that to download data you will need a user account. There are several tutorials available on the website through which the user can get familiar with the JGI workspace and downloading data. This is the <a href = "https://img.jgi.doe.gov/docs/IMGWorkspaceUserGuide.pdf"> guide </a> to using the workspace. The dataset used to generate sample files attached has IMG Genome ID: 3300050645 (<a href = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=3300050645"> link </a> to dataset).

#### 2.2: Use hmmsearch to get hits

We will use the hmm profile generated in part 1.3.3 to get o-demethylase putative hits from our dataset. The following is a sample code that cane be modified as per the users requirement. Other forms of output can also be generated. In this pipeline, we will specifically use the optional table format of output (--tblout option) to get the desired headers.

$ hmmsearch -o path/filename.txt -A path/msa_of_all_hits.sto --tblout path/filename.txt desired_hmm_profile_name.hmm sequence_database_to_search.faa

Alternatively, hmmscan can also be used but hmmsearch is faster

Following this the table output is imported in python and using the SearchIO package within biopython (refer to script <b>parsing_hmm_output.py</b>), the file is parsed to retrieve the headers (<b>parsing_hmm_input/output</b> folder). Next, depending on you co-localisation crietria (explained in "The pipeline" section previously) use the output as input to either of the co-localised header retrieval python scripts (NCBI or JGI_IMG_Protein_hit_co_localised_sequence_headers_retrieval.py). Sample input/output file formats for both scripts based on  datasets used as examples are included in the <b>python_header_metageome_output</b> folders.

#### 2.3: Retrieving co-localised sequence files from metagenome dataset

Once you have the co-localised header list, you can go back to the metagenomic data sets to retrieve the sequences. This is a good tool which can be run on Linux terminal which parses big metagenome datasets using a list of headers supplied as a text file (convert the headers in the csv file in the <b>python_function_metagenome_output</b> folder to separate text files and input to an executable tool called <a href = "https://github.com/santiagosnchez/faSomeRecords"> faSomeRecords </a> that accepts the list of headers as a simple text file. Any similar tool can be used to get the sequences. A more simplistic version based on "grep" command in Linux/bash can also be used using a -f flag for the header file input for patterns to be matched usually with a linearized fasta file (with all sequences in one line under each header). As an example, the co-localised fasta files retrieved from the NCBI metagenome dataset for co-localised MT1 hits are included (<b>fasta_sequence_retrieval_output_fasome</b>).

### Part 3: Visualisation in trees or network diagrams

Once we have the fasta sequences, it will be interesting to see the diversity of the functional genes within a metagenome or in comparison to the template and RefSeq reference sequences. This diversity in sequences can be visualised either using a phylogenetic tree or using a sequence similarity network (SSN) diagram. This will be helpful in understanding if this is a niche function limited to few taxonomic groups or if it is spread over a broad range of taxa. There are several ways to visualise relationships between sequences phylogenetically, two examples are given below. Having more than one diagram is also good for comparison.

###### Note : These steps are demonstrated with the combined fasta file for the MT1 proteins from the NCBI AD metagenome and the template and reference sequences (sample input and outputs to the different steps in part 3 are included in the MT1_AD_fasta_example including the fasta file, phyolgenetic tree output images, SSNPipe output, the network diagram input (igraph input) and the network diagram image.

#### 3.1 : Tree diagrams using phangorn and ape packages in R 

Phylogenetic trees are a biologists essential tool, and they were the first visualization method for understanding evolutionary relationships. There are several algorithms for tree construction now and they depend on the user research questions. Most of these algorithms analyse the substitutions in a given set of sequence alignment and infer phylogenetic or evolutionary distance accordingly. The <a href ="https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html#introduction"> resource </a> referenced in this example discusses uses the phangorn package which offers three algorithms for tree construction - neighbour-joining (nj), parsimony and maximum-likelihood (ml). The latter two are dependent on the construction of nj distance matrix first. Nj by itself is quite fast but will not be used because it relies on defining closest nodes as neighbours and keeps joining them together until tree construction is complete, but this distance joining is not guaranteed to produce the most accurate tree since changing the distance calculation algorithm also changes the tree. 

Earlier in the tutorial (section 1.3.1) we saw that we can use the msa package within R to generate the alignment. And since the tree construction is also within R, the alignment object can be passed directly, making this quite seamless. The ml algorithm is used for tree construction. Refer to the "phylogenetic_tree_construction_with_R" markdown file. The tree diagrams are also provided as image files in the same folder as the fasta file.

### 3.2 : Network diagrams using open-source <a href = "https://github.com/ahvdk/SSNpipe"> SSNpipe </a> tool and igraph package in R

SSNs enable the visualization of very large protein datasets without the need for multiple sequence alignment by doing individual pairwise comparisons. This is very fast computationally and allows inference of orthogonal information, or functional trends, in very large sets of evolutionarily related proteins. They give more control on the user end and are very useful when it is not practical to align very large dataset of sequences (such as those from metagenomes) [4]. The output from the tools used to create SSNs allow for easy visualization and can be streamlined with genomic neighbourhood maps or pfam families to enable the inferring of potential function of proteins with unknown annotations [5]. 

To generate the files needed for visualisation of the network diagrams, an open-source software called SSNPipe encoded in python was used. The combined fasta files with metagenomic and reference sequences, same as the one used in section 3.1 was used as the input to this software. The run parameters together with the instructions for using igraph package in R to visualise the network is included in R markdown file "network_diagrams_in_R".

If the user compares the ml_tree_color image with the MT1_network diagram image, one can appreciate the similarity in both the diagrams. There are at least 4 MT1 gene sequences from the AD metagenome which are quite distinct from the reference sequences and form exciting targets for heterologous production and functional testing.

###### Note : You can also use the open-source tools such as <a href = "https://cytoscape.org/"> Cytoscape </a> to make more complex network diagrams using the output from SSNPipe.

## Conclusion

This tutorial provides an easily-adaptable approach to search for co-localised o-demethylase sequences in the user's metagenomic dataset. The detailed instructions make it easy for beginners to follow along while more advanced users can adapt this workflow into a seamless tool for customised search for co-localised genes or genes related in function that occur in tandem in biological genomes using either publicly available dataset reference sequences or their own datasets.

## References
[1] J. Farnberger et al., "Biocatalytic methylation and demethylation via a shuttle catalysis concept involving corrinoid proteins", Communications Chemistry, vol. 1, no. 1, 2018. Available: 10.1038/s42004-018-0083-2.

[2] F. Mingo, S. Studenik and G. Diekert, "Conversion of phenyl methyl ethers byDesulfitobacteriumspp. and screening for the genes involved", FEMS Microbiology Ecology, vol. 90, no. 3, pp. 783-790, 2014. Available: 10.1111/1574-6941.12433.

[3] A. Schilhabel et al., "The Ether-Cleaving Methyltransferase System of the Strict Anaerobe Acetobacterium dehalogenans: Analysis and Expression of the Encoding Genes", Journal of Bacteriology, vol. 191, no. 2, pp. 588-599, 2008. Available: 10.1128/jb.01104-08.

[4] H. J. Atkinson, J. H. Morris, T. E. Ferrin, and P. C. Babbitt, “Using Sequence Similarity Networks for Visualization of Relationships Across Diverse Protein Superfamilies,” PLoS ONE, vol. 4, no. 2, 2009. 

[5] R. Zallot, N. Oberg, and J. A. Gerlt, “The EFI Web Resource for Genomic Enzymology Tools: Leveraging Protein, Genome, and Metagenome Databases to Discover Novel Enzymes and Metabolic Pathways,” Biochemistry, vol. 58, no. 41, pp. 4169–4182, 2019.
