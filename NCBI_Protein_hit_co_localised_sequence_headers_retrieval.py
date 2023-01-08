def co_localised_operon_hunt():
    CP_headers = input('Enter_path_to_CP_header_file: ')
    MT1_headers = input('Enter_path_to_MT1header_file:')
    MT2_headers = input('Enter_path_to_MT2header_file:')
    
    CP_headers_in = open(CP_headers)
    MT1_headers_in = open(MT1_headers)
    MT2_headers_in = open(MT2_headers)
    
    df_CP_headers_all = pd.read_csv(CP_headers_in)
    df_CP_headers_all.columns = ['Sequence_ID_CP']
    
    df_MT1_headers_all = pd.read_csv(MT1_headers_in)
    df_MT1_headers_all.columns = ['Sequence_ID_MT1']
    
    df_MT2_headers_all = pd.read_csv(MT2_headers_in)
    df_MT2_headers_all.columns = ['Sequence_ID_MT2']
    
    df_operon_merged = [df_CP_headers_all['Sequence_ID_CP'], df_MT1_headers_all['Sequence_ID_MT1'], df_MT2_headers_all['Sequence_ID_MT2']]
    df_operon_hunt = pd.concat(df_operon_merged, axis = 1)
    
    TBA = 'KAF|\.1' #defining regex for replacement
    df_operon_hunt["Numeric_ID_CP"] = df_operon_hunt["Sequence_ID_CP"].str.replace(TBA, "")
    df_operon_hunt["Numeric_ID_CP"] = df_operon_hunt["Numeric_ID_CP"].astype(int)
    
    #both MT1 and MT2 empty rows will be filled with 0s
    
    df_operon_hunt['Sequence_ID_MT1'] = df_operon_hunt['Sequence_ID_MT1'].fillna(0)
    df_operon_hunt["Numeric_ID_MT1"] = df_operon_hunt["Sequence_ID_MT1"].str.replace(TBA, "")
    df_operon_hunt['Numeric_ID_MT1'] = df_operon_hunt['Numeric_ID_MT1'].fillna(0)
    df_operon_hunt["Numeric_ID_MT1"] = df_operon_hunt["Numeric_ID_MT1"].astype(int)
    
    df_operon_hunt['Sequence_ID_MT2'] = df_operon_hunt['Sequence_ID_MT2'].fillna(0)
    df_operon_hunt["Numeric_ID_MT2"] = df_operon_hunt["Sequence_ID_MT2"].str.replace(TBA, "")
    df_operon_hunt['Numeric_ID_MT2'] = df_operon_hunt['Numeric_ID_MT2'].fillna(0)
    df_operon_hunt["Numeric_ID_MT2"] = df_operon_hunt["Numeric_ID_MT2"].astype(int)
    
    #working up until this step
    
    #Now time to find those co-localised indices
    #remeber to use the protein list with highest number of hits as basis for loop, in this case it was CP, update as per your hits result
    #adjust numerical range to suit your data, the code is looks for extremely tight co-localization +/- 1
    
    ctr_CP = 0
    CP_numeric_IDs = []
    ctr_MT2 = 0
    MT2_numeric_IDs= []
    MT2_IDs  = []
    CP_1_numeric_IDs = []
    ctr_MT1 = 0
    MT1_numeric_IDs = []

    for CP in df_operon_hunt['Numeric_ID_CP']:
        ctr_MT2 = 0
        for MT2 in df_operon_hunt['Numeric_ID_MT2']:
            if((CP==(MT2+1))|(CP==(MT2-1))): 
                CP_numeric_IDs.append(CP)
                MT2_numeric_IDs.append(MT2) 
                
            ctr_MT2 = ctr_MT2+1
        ctr_MT1 = 0
    
        for MT1 in df_operon_hunt['Numeric_ID_MT1']:
            if((CP==(MT1+1))|(CP==(MT1-1))):
                CP_1_numeric_IDs.append(CP)
                MT1_numeric_IDs.append(MT1)
                
            ctr_MT1 = ctr_MT1+1
        ctr_CP = ctr_CP+1
        
#our golden list is ready, now we need to convert back to string form to get complete IDs

    ls_CP = map(str, CP_numeric_IDs)
    ls_MT2 = map(str, MT2_numeric_IDs)
    ls_MT1 = map(str, MT1_numeric_IDs)
    ls_CP_1 = map(str, CP_1_numeric_IDs)
    
    prefix = 'KAF'
    suffix = '.1'
    
    ls_CP_ID = [prefix + s + suffix for s in ls_CP]
    ls_MT2_ID = [prefix + s + suffix for s in ls_MT2]
    ls_MT1_ID = [prefix + s + suffix for s in ls_MT1]
    ls_CP_1_ID = [prefix + s + suffix for s in ls_CP_1]
    ls_CP_all = ls_CP_ID + ls_CP_1_ID
    ls_CP_triplicates = list(set([x for x in ls_CP_all if ls_CP_all.count(x) > 1]))
    ls_all_headers = [ls_CP_ID,ls_MT2_ID,ls_MT1_ID,ls_CP_1_ID,ls_CP_triplicates]
    
    big_header_list_df = pd.DataFrame(ls_all_headers, index = ['CP_paired_MT2', 'MT2_paired_CP', 'MT1_paired_CP', 'CP_paired_MT1', 'CP_triplicates'])
    big_header_list_df.to_csv("all_co_localised_headers.csv")
    return("The file called all_co_localised_headers.csv is ready in your working directory")