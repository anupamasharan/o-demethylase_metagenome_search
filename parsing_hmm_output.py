def hmmer_input_parse():
    CP_hmmsearch = input('Enter path to CP hmmsearch table text file:')
    MT1_hmmsearch = input('Enter path to MT1 hmmsearch table text file:')
    MT2_hmmsearch = input('Enter path to MT2 hmmsearch table text file:')
    
    CP_hmmsearch_in = open(CP_hmmsearch)
    MT1_hmmsearch_in = open(MT1_hmmsearch)
    MT2_hmmsearch_in = open(MT2_hmmsearch)
    
    from Bio import SearchIO
    import pandas as pd
    
    CP_metagenome_headers = []
    MT1_metagenome_headers = []
    MT2_metagenome_headers = []

    with CP_hmmsearch_in as handle: 
        for record in SearchIO.parse(handle, 'hmmer3-tab'):

            query_id = record.id #seqID from fasta 
            hits = record.hits 
            num_hits = len(hits) # number of hits per query

            if num_hits > 0: # if there are more than 0 hits per query then we need to extract the info 
                for i in range(0,num_hits): 
                    CP_metagenome_headers.append(hits[i].id) # hit name
                #for more information you could use the same handle and extract the other information
                #hmm_description = hits[i].description # hit decription 
                #current_evalue = hits[i].evalue # evalue of hit
                
        handle.close()
    
    with MT1_hmmsearch_in as handle: 
        for record in SearchIO.parse(handle, 'hmmer3-tab'):

            query_id = record.id #seqID from fasta 
            hits = record.hits 
            num_hits = len(hits) # number of hits per query

            if num_hits > 0: # if there are more than 0 hits per query then we need to extract the info 
                for i in range(0,num_hits): 
                    MT1_metagenome_headers.append(hits[i].id) # hit name
                #for more information you could use the same handle and extract the other information
                #hmm_description = hits[i].description # hit decription 
                #current_evalue = hits[i].evalue # evalue of hit
                
        handle.close()
    
    with MT2_hmmsearch_in as handle: 
        for record in SearchIO.parse(handle, 'hmmer3-tab'):

            query_id = record.id #seqID from fasta 
            hits = record.hits 
            num_hits = len(hits) # number of hits per query

            if num_hits > 0: # if there are more than 0 hits per query then we need to extract the info 
                for i in range(0,num_hits): 
                    MT2_metagenome_headers.append(hits[i].id) # hit name
                #for more information you could use the same handle and extract the other information
                #hmm_description = hits[i].description # hit decription 
                #current_evalue = hits[i].evalue # evalue of hit
                
        handle.close()
        
    import os
    os.makedirs('co_localised_sequence_search_input')
    os.chdir('co_localised_sequence_search_input') 
    #WARNING: The working directory is changed at this step! Make sure to change it back to the working directory.
    CP_metagenome_headers_df = pd.DataFrame(CP_metagenome_headers, columns = ['CP_ID'])
    CP_metagenome_headers_df.to_csv("CP_metagenome_headers.csv", index = False)
    
    MT1_metagenome_headers_df = pd.DataFrame(MT1_metagenome_headers, columns = ['MT1_ID'])
    MT1_metagenome_headers_df.to_csv("MT1_metagenome_headers.csv", index = False)
    
    MT2_metagenome_headers_df = pd.DataFrame( MT2_metagenome_headers, columns = ['MT2_ID'])
    MT2_metagenome_headers_df.to_csv("MT2_metagenome_headers.csv", index = False)
    
    return(print("The files are ready in the co_localised_sequence_search_input"))

hmmer_input_parse()