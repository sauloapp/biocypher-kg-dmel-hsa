import pandas as pd
import glob
import requests
import os
import csv
import psycopg2
from flybase_tsv_reader import FlybasePrecomputedTable
import gzip

out_dir = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/fca2/gene_data/"
out_dir = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2"
symbols_to_fbgn_file = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/flybase/fbgn_fbtr_fbpp_expanded_fb_2024_03.tsv.gz"
symbols_to_fbgn_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/full/flybase/fbgn_fbtr_fbpp_expanded_fb_2024_03.tsv.gz"

def download_rnaseq_data(fbgn_id, out_dir="./"):
    if not (fbgn_id.startswith("FBgn") or fbgn_id.startswith("FBtr")):
        raise ValueError("FBgn_id must start with 'FBgn' or 'FBtr'")

    for gene_type in ['gene', 'transcriptGene', 'mir', 'transcriptMir']:
        url = f"https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn={fbgn_id}&tableOut={gene_type}"

        # Define the file name
        file_name = f"{out_dir}{fbgn_id}_{gene_type}.tsv"

        # Check if the file already exists
        if os.path.exists(file_name):
            print(f"File {file_name} already exists. Skipping download.")
            continue

        # Make a request to download the data
        response = requests.get(url)

        if response.status_code == 200:
            if "An error has occurred" in response.text:
                print(f'No data for {fbgn_id}, type {gene_type}')
                continue
            # Save the content to a file if it doesn't already exist
            with open(file_name, 'w') as file:
                file.write(response.text)
            print(f"Data saved to {file_name}")
        else:
            print(f"Failed to download data: {response.status_code}")

# download_rnaseq_data("FBgn0262455", out_dir)  # mirna different output format
# exit(9 )


def download_all_fbgn_data():
    symbols_to_fbgn_table = FlybasePrecomputedTable(symbols_to_fbgn_file)
    dframe = symbols_to_fbgn_table.to_pandas_dataframe()
    fbgn_list = dframe['gene_ID'].tolist()
    print(len(fbgn_list))
    fbgn_list = dframe['gene_ID'].unique().tolist()
    print(len(fbgn_list))
    print(f'{fbgn_list[0]} {fbgn_list[-1]}')
    #exit(9)

    #for gene in net_act_expanded_genes_list:
    for gene in fbgn_list:
        if gene.startswith("FBgn"):
            download_rnaseq_data(gene)


def convert_transcriptGene_files(file_list, output_file, transcript_type='regularRNA'):
    # Initialize an empty list to store dataframes
    dfs = []

    for file in file_list:
        # Read the file assuming tab-separated values (TSV)
        df = pd.read_csv(file, sep='\t', skiprows=6)

        with open(file) as f:
            metadata = [next(f).strip().split('\t')[1] for _ in range(4)]

        fb_gene_id, annotation_symbol, symbol, name, = metadata
        with open(file) as f:
            [next(f) for _ in range(5)]
            metadata = next(f).strip().split('\t')
        transcripts = [id for id in metadata if id != '']
        #print(transcripts)

        new_columns = ['Tissue stage and sex', 'Tissue']
        for i in range(2, len(df.columns)):
            new_columns.append(df.columns[i])
        print(df.columns)
        print(new_columns)
        df.columns = new_columns

        for t in range( len(transcripts) ):
            for index, row in df.iterrows():
                tissue_stage_and_sex = row['Tissue stage and sex']
                tissue = row['Tissue']
                #for i in range(2, (len(df.columns) - 1)):
                    # FlyBase ID	Transcript ID	Tissue stage sex	Tissue	FPKM	SD
                #print(f'Row: {row}')
                print(f'transcript: {transcripts[t]}\tt: {t}\tcolumn: : {df.columns[2 * t + 2]}')
                print(f'transcript: {transcripts[t]}\tt: {t}\tcolumn: : {df.columns[2 * t + 3]}')
                if transcript_type == 'microRNA':
                    metrics = 'TPM'
                elif transcript_type == 'regularRNA':
                    metrics = 'FPKM'
                dfs.append({
                    '#FBgene ID': fb_gene_id,                    
                    'Tissue stage and sex': tissue_stage_and_sex,
                    'Tissue': tissue,
                    'FBtranscript ID': transcripts[t],
                    metrics: row[ df.columns[2 * t + 2] ],
                    'SD': row[ df.columns[2 * t + 3] ]
                })
                #print(dfs)

    final_df = pd.DataFrame(dfs)
    # Write the final dataframe to a gzipped TSV file
    with gzip.open(output_file, 'wt') as f:
        final_df.to_csv(f, sep='\t', index=False)

file_list = glob.glob("/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/transcripts/*_transcriptGene.tsv")
output_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_transcriptGene_output.tsv.gz"
convert_transcriptGene_files(file_list, output_file, transcript_type='regularRNA')

file_list = glob.glob("/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/micro_rna_transcripts/*_transcriptMir.tsv")
output_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_Mir_transcript_output.tsv.gz"
convert_transcriptGene_files(file_list, output_file, transcript_type='microRNA')



def get_library_names(uniquenames):
    # Establish a connection to the FlyBase database
    conn = psycopg2.connect(
        host="chado.flybase.org",
        database="flybase",
        user="flybase"
    )
    
    # Create a cursor object
    cur = conn.cursor()
    
    # Initialize a dictionary to store results
    results = {}
    
    # Loop through each uniquename in the provided list
    for uniquename in uniquenames:
        # Execute the query to retrieve the name of the library
        cur.execute("SELECT name FROM library WHERE library.uniquename LIKE %s;", (uniquename,))
        # Fetch the result
        result = cur.fetchone()
        # If a result is found, add it to the dictionary
        if result:
            results[uniquename] = result[0]
        else:
            results[uniquename] = None
    
    # Close the cursor and connection
    cur.close()
    conn.close()
    
    return results


uniquenames = ['FBlc0003619', 'FBlc0003620', 'FBlc0003621', 'FBlc0003622', 'FBlc0003623', 
               'FBlc0003624', 'FBlc0003625', 'FBlc0003626', 'FBlc0003627', 'FBlc0003628', 
               'FBlc0003629', 'FBlc0003630', 'FBlc0003631', 'FBlc0003632', 'FBlc0003633', 
               'FBlc0003634', 'FBlc0003724', 'FBlc0003650', 'FBlc0003651', 'FBlc0003652', 
               'FBlc0003653', 'FBlc0003654', 'FBlc0003655', 'FBlc0003656', 'FBlc0003657', 
               'FBlc0003658', 'FBlc0003685', 'FBlc0003635', 'FBlc0003636', 'FBlc0003637', 
               'FBlc0003638', 'FBlc0003639', 'FBlc0003640', 'FBlc0003641', 'FBlc0003642', 
               'FBlc0003643', 'FBlc0003644', 'FBlc0003645', 'FBlc0003646', 'FBlc0003647', 
               'FBlc0003648', 'FBlc0003649', 'FBlc0003726']
# print(get_library_names(uniquenames))
# exit(9)


def build_fca2_fb_tissues_libraries_ids_dicts(file_path):
    gene_tissue_library_dict = {}
    transcript_tissue_library_dict = {}

    # Establish a connection to the FlyBase database
    conn = psycopg2.connect(
        host="chado.flybase.org",
        database="flybase",
        user="flybase"
    )    
    cur = conn.cursor()
    #uniquename = row['fb_library_id']
    # Execute the query to retrieve the (uniqueame, name) of the library
    #cur.execute("SELECT name FROM library WHERE library.uniquename LIKE %s;", (uniquename,))
    cur.execute("SELECT uniquename, name FROM library WHERE library.name LIKE 'RNA-Seq_FlyAtlas2_%';")
    results = cur.fetchall() 
    mir_cur = conn.cursor()
    mir_cur.execute("SELECT uniquename, name FROM library WHERE library.name LIKE 'microRNA-Seq_TPM_FlyAtlas2_%';")
    #mir_cur.execute("SELECT uniquename, name FROM library WHERE library.name LIKE 'microRNA-Seq%';")
    mir_results = mir_cur.fetchall()
    print(mir_results)
    print(len(mir_results))
    
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        # row:
        # Tissue stage and sex	    _gene file tissue	    _transcriptGene file tissue	    FB_tissue
        for row in reader:
            for result in results:
                if result[-1].endswith(row['FB_tissue'].replace(' ', '_')):
                    library_name = result[-1] if result else None
                    library_uniquename = result[0] if result else None
                    gene_key = f"{row['Tissue stage and sex']}_{row['_gene file tissue']}"
                    gene_value = (row['FB_tissue'], library_uniquename, library_name)

                    # Add the value to the tissue_library_dict
                    if gene_key not in gene_tissue_library_dict:
                        gene_tissue_library_dict[gene_key] = gene_value

                    transcript_key = f"{row['Tissue stage and sex']}_{row['_transcriptGene file tissue']}"
                    transcript_value = (row['FB_tissue'], library_uniquename, library_name)

                    # Add the value to the transcript_dict
                    if transcript_key not in transcript_tissue_library_dict:
                        transcript_tissue_library_dict[transcript_key] = transcript_value
                    break
            for result in mir_results:
                if result[-1].endswith(row['FB_tissue'].replace(' ', '_')):
                    library_name = result[-1] if result else None
                    library_uniquename = result[0] if result else None
                    gene_key = f"microRNA_{row['Tissue stage and sex']}_{row['_gene file tissue']}"
                    gene_value = (row['FB_tissue'], library_uniquename, library_name)

                    # Add the value to the tissue_library_dict
                    if gene_key not in gene_tissue_library_dict:
                        gene_tissue_library_dict[gene_key] = gene_value

                    transcript_key = f"microRNA_{row['Tissue stage and sex']}_{row['_transcriptGene file tissue']}"
                    transcript_value = (row['FB_tissue'], library_uniquename, library_name)

                    # Add the value to the transcript_dict
                    if transcript_key not in transcript_tissue_library_dict:
                        transcript_tissue_library_dict[transcript_key] = transcript_value
                    break
            
    cur.close()
    mir_cur.close()
    conn.close()

    return gene_tissue_library_dict, transcript_tissue_library_dict


#file_path = '/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/bizon/biocypher-kg-dmel-hsa/aux_files/fca2_to_fb_tissues.tsv'
file_path = 'aux_files/fca2_to_fb_tissues.tsv'
gene_dict, transcript_dict = build_fca2_fb_tissues_libraries_ids_dicts(file_path)
#fca2_lib_to_fb_tissue_lib_id
# Display the resulting dictionaries
print("Gene Dictionary:")
print(gene_dict)
print("\nTranscript Dictionary:")
print(transcript_dict)
exit(9)

def convert_gene_files(file_list, output_file, gene_type='regular_gene'):
    # List to store individual DataFrames
    data_frames = []

    # Iterate over each file in the list
    for file in file_list:
        # Read the TSV file, ignoring the first lines of metadata
        df = pd.read_csv(file, sep='\t', skiprows=6)

        # Name the columns correctly
        columns = ['Tissue', 'FPKM_Male_Adult', 'SD_Male_Adult', 'Enrichment_Male_Adult',
                   'FPKM_Female_Adult', 'SD_Female_Adult', 'Enrichment_Female_Adult',
                   'M/F', 'p_value', 'FPKM_Larval', 'SD_Larval', 'Enrichment_Larval']
        df.columns = columns

        # Extract metadata information from the ignored lines
        with open(file) as f:
            metadata = [next(f).strip().split('\t')[1] for _ in range(4)]

        flybase_id, annotation_symbol, symbol, name = metadata

        # Process each row of the DataFrame
        for index, row in df.iterrows():
            if row['Tissue'] != '-':
                if gene_type == 'regular_gene':
                    metrics = 'FPKM'
                elif gene_type == 'microRNA':
                    metrics = 'TPM'
                # Add data for "Male Adult"
                if row['FPKM_Male_Adult'] != '-':
                    data_frames.append({
                        '#FBgene ID': flybase_id,
                        #'Annotation Symbol': annotation_symbol,
                        #'Symbol': symbol,                        
                        'Tissue stage and sex': 'Adult Male',
                        'Tissue': row['Tissue'],
                        metrics: row['FPKM_Male_Adult'],
                        'SD': row['SD_Male_Adult'],
                        'Enrichment': row['Enrichment_Male_Adult']
                    })

                # Add data for "Female Adult"
                if row['FPKM_Female_Adult'] != '-':
                    data_frames.append({
                        '#FBgene ID': flybase_id,
                        #'Annotation Symbol': annotation_symbol,
                        #'Symbol': symbol,                        
                        'Tissue stage and sex': 'Adult Female',                        
                        'Tissue': row['Tissue'],
                        metrics: row['FPKM_Female_Adult'],
                        'SD': row['SD_Female_Adult'],
                        'Enrichment': row['Enrichment_Female_Adult']
                    })

                # Add data for "Larval"
                if row['FPKM_Larval'] != '-':
                    data_frames.append({
                        '#FBgene ID': flybase_id,
                        #'Annotation Symbol': annotation_symbol,
                        #'Symbol': symbol,                        
                        'Tissue stage and sex': 'Larval',                        
                        'Tissue': row['Tissue'],
                        metrics: row['FPKM_Larval'],
                        'SD': row['SD_Larval'],
                        'Enrichment': row['Enrichment_Larval']
                    })

    # Create a final DataFrame with all consolidated data
    final_df = pd.DataFrame(data_frames)

    # Save the final DataFrame to a gzipped TSV file
    with gzip.open(output_file, 'wt') as f:
        final_df.to_csv(f, sep='\t', index=False)


file_list = glob.glob("/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/genes/*.tsv")
output_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_gene_output.tsv.gz"
convert_gene_files(file_list, output_file, gene_type='regular_gene')

file_list = glob.glob("/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/micro_rna_genes/*.tsv")
output_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_Mir_gene_output.tsv.gz"
convert_gene_files(file_list, output_file, gene_type='microRNA')

# download_rnaseq_data("FBgn0262177", "gene")  # mirna different output format
