def co_localised_operon_hunt():
    CP_headers = input('Enter path to CP header_file: ') #use only file name if it is already in working directory
    MT1_headers = input('Enter path to MT1_header_file:')
    MT2_headers = input('Enter path to MT2 header_file:')
    
    CP_headers_in = open(CP_headers)
    MT1_headers_in = open(MT1_headers)
    MT2_headers_in = open(MT2_headers)
    
    import pandas as pd
    
    df_CP_headers_all = pd.read_csv(CP_headers_in)
    df_CP_headers_all.columns = ['Sequence_header_CP']
    
    #separating the header using the unique '_' delimiter found in JGI protein .faa file headers, update to your specific sequence header to get contig info
    df_CP_headers_all[['CP_JGI_Project_ID','CP_Contig_ID','CP_Gene_Start','CP_Gene_End']] = df_CP_headers_all.Sequence_header_CP.str.split("_", expand = True)
    
    df_MT1_headers_all = pd.read_csv(MT1_headers_in)
    df_MT1_headers_all.columns = ['Sequence_header_MT1']
    df_MT1_headers_all[['MT1_JGI_Project_ID','MT1_Contig_ID','MT1_Gene_Start','MT1_Gene_End']] = df_MT1_headers_all.Sequence_header_MT1.str.split("_", expand = True)
    
    df_MT2_headers_all = pd.read_csv(MT2_headers_in)
    df_MT2_headers_all.columns = ['Sequence_header_MT2']
    df_MT2_headers_all[['MT2_JGI_Project_ID','MT2_Contig_ID','MT2_Gene_Start','MT2_Gene_End']] = df_MT2_headers_all.Sequence_header_MT2.str.split("_", expand = True)
    
    #time to find sequences co-localised on the same contig

    ctr_MT1 = 0 #it is necessary to start with the protein that has the highest number of hits so that co-localised hits are not missed, update based on your search results
    MT1_index_IDs = []
    MT1_1_index_IDs = []
    ctr_CP = 0
    CP_index_IDs= []
    ctr_MT2 = 0
    MT2_index_IDs = []

    for x in df_MT1_headers_all['MT1_Contig_ID']:
        ctr_CP = 0
        for y in df_CP_headers_all['CP_Contig_ID']:
            if(x == y): 
                MT1_index_IDs.append(ctr_MT1)
                CP_index_IDs.append(ctr_CP) 
                
            ctr_CP = ctr_CP+1
        ctr_MT2 = 0
    
        for z in df_MT2_headers_all['MT2_Contig_ID']:
            if(x == z):
                MT1_1_index_IDs.append(ctr_MT1)
                MT2_index_IDs.append(ctr_MT2) 
                
            ctr_MT2 = ctr_MT2+1
        ctr_MT1 = ctr_MT1+1

#we run an additional loop to find CP and MT2 pairs
                                
    ctr_CP_1 = 0
    CP_1_index_IDs = []
    ctr_MT2_1 = 0
    MT2_1_index_IDs = []

    for i in df_CP_headers_all['CP_Contig_ID']:
        ctr_MT2_1 = 0
        for j in df_MT2_headers_all['MT2_Contig_ID']:
            if(i == j): 
                CP_1_index_IDs.append(ctr_CP_1)
                MT2_1_index_IDs.append(ctr_MT2_1)   
            
            ctr_MT2_1 = ctr_MT2_1+1
        ctr_CP_1 = ctr_CP_1+1
                                
#remove duplicate index instances from first round of contig ID matching
    
    MT1_unique_index = list(set(MT1_index_IDs))
    MT1_1_unique_index = list(set(MT1_1_index_IDs))
    CP_unique_index = list(set(CP_index_IDs))
    CP_1_unique_index = list(set(CP_1_index_IDs))
    MT2_unique_index = list(set(MT2_index_IDs))
    MT2_1_unique_index = list(set(MT2_1_index_IDs))

    MT1_all_index = MT1_unique_index + MT1_1_unique_index
    CP_all_index = CP_unique_index + CP_1_unique_index
    MT2_all_index = MT2_unique_index + MT2_1_unique_index

#subsetting list to get indices to retrieve headers from dataframe
#we are focusing on triplicates or contings that have at least one copy each of MT1, CP and MT2 
#remove comment from pair sentences if you also need paired headers 

    ls_MT1_triplicates = list(set([x for x in MT1_all_index if MT1_all_index.count(x) > 1]))
    #ls_MT1_pairs = list(set([x for x in MT1_all_index if MT1_all_index.count(x) == 1]))
    ls_CP_triplicates = list(set([x for x in CP_all_index if CP_all_index.count(x) > 1]))
    #ls_CP_pairs = list(set([x for x in CP_all_index if CP_all_index.count(x) == 1]))
    ls_MT2_triplicates = list(set([x for x in MT2_all_index if MT2_all_index.count(x) > 1]))
    #ls_MT2_pairs = list(set([x for x in MT2_all_index if MT2_all_index.count(x) == 1]))
    
#retreiving desired headers as list

    MT1_triplicate_headers = pd.DataFrame(list(df_MT1_headers_all.iloc[ls_MT1_triplicates,0]), columns = ['MT1_triplicate_headers'])
    MT1_triplicate_headers.to_csv("MT1_metagenome_triplicate_headers.csv", index = False)
    MT1_triplicate_headers[['MT1_JGI_Project_ID','MT1_Contig_ID','MT1_Gene_Start','MT1_Gene_End']] = MT1_triplicate_headers.MT1_triplicate_headers.str.split("_", expand = True)
    MT1_repeat_and_unique_ls = list(MT1_triplicate_headers['MT1_Contig_ID'])
    MT1_repeat_contig_IDs = list(set([x for x in MT1_repeat_and_unique_ls if MT1_repeat_and_unique_ls.count(x) > 1]))
    MT1_multicopy_contig_IDs = pd.DataFrame(MT1_repeat_contig_IDs, columns = ['MT1_multicopy_contig_IDs'])
    MT1_unique_contig_IDs = list(set([x for x in MT1_repeat_and_unique_ls if MT1_repeat_and_unique_ls.count(x) == 1]))
    MT1_single_copy_contig_IDs = pd.DataFrame(MT1_unique_contig_IDs, columns = ['MT1_single_copy_contig_IDs'])
    MT1_multicopy_contig_IDs.to_csv("MT1_metagenome_triplicates_multicopy_contig_IDs.csv", index = False)
    MT1_single_copy_contig_IDs.to_csv("MT1_metagenome_triplicates_single_copy_contig_IDs.csv", index = False)
    
                           
    #MT1_pair_headers = list(df_MT1_headers_all.iloc[ls_MT1_pairs,0])
    #MT1_pair_headers.to_csv("MT1_metagenome_pair_headers.csv", index = False)
    
    CP_triplicate_headers = pd.DataFrame(list(df_CP_headers_all.iloc[ls_CP_triplicates,0]), columns = ['CP_triplicate_headers'])
    CP_triplicate_headers.to_csv("CP_metagenome_triplicate_headers.csv", index = False)
    CP_triplicate_headers[['CP_JGI_Project_ID','CP_Contig_ID','CP_Gene_Start','CP_Gene_End']] = CP_triplicate_headers.CP_triplicate_headers.str.split("_", expand = True)
    CP_repeat_and_unique_ls = list(CP_triplicate_headers['CP_Contig_ID'])
    CP_repeat_contig_IDs = list(set([x for x in CP_repeat_and_unique_ls if CP_repeat_and_unique_ls.count(x) > 1]))
    CP_multicopy_contig_IDs = pd.DataFrame(CP_repeat_contig_IDs, columns = ['CP_multicopy_contig_IDs'])
    CP_unique_contig_IDs = list(set([x for x in CP_repeat_and_unique_ls if CP_repeat_and_unique_ls.count(x) == 1]))
    CP_single_copy_contig_IDs = pd.DataFrame(CP_unique_contig_IDs, columns = ['CP_single_copy_contig_IDs'])
    CP_multicopy_contig_IDs.to_csv("CP_metagenome_triplicates_multicopy_contig_IDs.csv", index = False)
    CP_single_copy_contig_IDs.to_csv("CP_metagenome_triplicates_single_copy_contig_IDs.csv", index = False)
                           
    #CP_pair_headers = list(df_CP_headers_all.iloc[ls_MT1_pairs,0])
    #CP_pair_headers.to_csv("CP_metagenome_pair_headers.csv", index = False)
                           
    MT2_triplicate_headers = pd.DataFrame(list(df_MT2_headers_all.iloc[ls_MT2_triplicates,0]), columns = ['MT2_triplicate_headers'])
    MT2_triplicate_headers.to_csv("MT2_metagenome_triplicate_headers.csv", index = False)
    MT2_triplicate_headers[['MT2_JGI_Project_ID','MT2_Contig_ID','MT2_Gene_Start','MT2_Gene_End']] = MT2_triplicate_headers.MT2_triplicate_headers.str.split("_", expand = True)
    MT2_repeat_and_unique_ls = list(MT2_triplicate_headers['MT2_Contig_ID'])
    MT2_repeat_contig_IDs = list(set([x for x in MT2_repeat_and_unique_ls if MT2_repeat_and_unique_ls.count(x) > 1]))
    MT2_multicopy_contig_IDs = pd.DataFrame(MT2_repeat_contig_IDs, columns = ['MT2_multicopy_contig_IDs'])
    MT2_unique_contig_IDs = list(set([x for x in MT2_repeat_and_unique_ls if MT2_repeat_and_unique_ls.count(x) == 1]))
    MT2_single_copy_contig_IDs = pd.DataFrame(MT2_unique_contig_IDs, columns = ['MT2_single_copy_contig_IDs'])
    MT2_multicopy_contig_IDs.to_csv("MT2_metagenome_triplicates_multicopy_contig_IDs.csv", index = False)
    MT2_single_copy_contig_IDs.to_csv("MT2_metagenome_triplicates_single_copy_contig_IDs.csv", index = False)
    
    #MT2_pair_headers = list(df_MT2_headers_all.iloc[ls_MT2_pairs,0])
    #MT2_pair_headers.to_csv("MT2_metagenome_pair_headers.csv", index = False)
    
    return(print("The co-localised header files and contig IDs with repeat gene copies are ready in your working directory"))
    
co_localised_operon_hunt()