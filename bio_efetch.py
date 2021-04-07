def retrieving_fasta_files_online():

    all_headers = input('Enter_path_to_the_co-localised_header_file: ')
   
    all_headers_in = open(all_headers)
    
    import pandas as pd
    
    df_header_final = pd.read_csv(all_headers_in)
    
    #we need to tidy the dataframe
    
    df_header_final = df_header_final.rename(index={0: 'CP_paired_MT2', 1: 'MT2_paired_CP', 2: 'MT1_paired_CP', 3: 'CP_paired_MT1', 4: 'CP_triplicates', 5: 'CP_pairs_only'})
    df_header_final = df_header_final.drop(df_header_final.columns[[0]], axis=1)
    
    #separating the relevant IDs out
    
    df_MT2 = df_header_final.loc['MT2_paired_CP']
    df_MT2 = df_MT2.dropna()
    df_MT1 = df_header_final.loc['MT1_paired_CP']
    df_MT1 = df_MT1.dropna()
    df_CP_triplets = df_header_final.loc['CP_triplicates']
    df_CP_triplets = df_CP_triplets.dropna()
    df_CP_pairs = df_header_final.loc['CP_pairs_only']
    df_CP_pairs = df_CP_pairs.dropna()
    
    #making list for input to the efetch function
    
    MT2_IDs = df_MT2.values.tolist()
    MT1_IDs = df_MT1.values.tolist()
    CP_IDs_triplets = df_CP_triplets.values.tolist()
    CP_IDs_pairs = df_CP_pairs.values.tolist()
    CP_IDs_all = CP_IDs_triplets + CP_IDs_pairs #since the co-localised operon search function was anchored on CP, the final list is re-adjusted to make sure there are no duplicate values
   
    from Bio import Entrez

    #Entrez search with efetch

    Entrez.email = "anupama.sharan@mail.utoronto.ca" #update with your email address

    handle_CP = Entrez.efetch(db="protein", id=",".join(CP_IDs_all), rettype="fasta", retmode="text")
    handle_MT1 = Entrez.efetch(db="protein", id=",".join(MT1_IDs), rettype="fasta", retmode="text")
    handle_MT2 = Entrez.efetch(db="protein", id=",".join(MT2_IDs), rettype="fasta", retmode="text")
    
    CP_ID_fasta = open("fasta_files_output/all_CP.fasta", "w") 
    CP_ID_fasta.write(handle_CP.read())
    CP_ID_fasta.close()

    MT1_ID_fasta = open("fasta_files_output/all_MT1.fasta", "w") 
    MT1_ID_fasta.write(handle_MT1.read())
    MT1_ID_fasta.close()

    MT2_ID_fasta = open("fasta_files_output/all_MT2.fasta", "w") 
    MT2_ID_fasta.write(handle_MT2.read())
    MT2_ID_fasta.close()
    
    return(print("The reference database fasta files are ready in the fasta_files_output directory "))

retrieving_fasta_files_online()
#and we're done....for now!