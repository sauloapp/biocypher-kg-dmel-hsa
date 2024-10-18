# -*- coding: utf-8 -*-
import pandas as pd
import glob
import requests
import os
import csv
import psycopg2
#from flybase_tsv_reader import FlybasePrecomputedTable
import gzip
import re

out_dir = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/fca2/gene_data/"
out_dir = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2"
symbols_to_fbgn_file = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/flybase/fbgn_fbtr_fbpp_expanded_fb_2024_03.tsv.gz"
symbols_to_fbgn_file = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/full/flybase/fbgn_fbtr_fbpp_expanded_fb_2024_03.tsv.gz"

# def download_rnaseq_data(fbgn_id, out_dir="./"):
#     if not (fbgn_id.startswith("FBgn") or fbgn_id.startswith("FBtr")):
#         raise ValueError("FBgn_id must start with 'FBgn' or 'FBtr'")

#     for gene_type in ['gene', 'transcriptGene', 'mir', 'transcriptMir']:
#         url = f"https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn={fbgn_id}&tableOut={gene_type}"

#         # Define the file name
#         file_name = f"{out_dir}{fbgn_id}_{gene_type}.tsv"

#         # Check if the file already exists
#         if os.path.exists(file_name):
#             print(f"File {file_name} already exists. Skipping download.")
#             continue

#         # Make a request to download the data
#         response = requests.get(url)

#         if response.status_code == 200:
#             if "An error has occurred" in response.text:
#                 print(f'No data for {fbgn_id}, type {gene_type}')
#                 continue
#             # Save the content to a file if it doesn't already exist
#             with open(file_name, 'w') as file:
#                 file.write(response.text)
#             print(f"Data saved to {file_name}")
#         else:
#             print(f"Failed to download data: {response.status_code}")            


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
                
                # Append the fbgn_id and gene_type to the log file
                with open(log_file_path, 'a') as log_file:
                    log_file.write(f"{fbgn_id}\t{gene_type}\n")
                    
                continue
            # Save the content to a file if it doesn't already exist
            with open(file_name, 'w') as file:
                file.write(response.text)
            print(f"Data saved to {file_name}")
        else:
            print(f"Failed to download data: {response.status_code}")


# Path to the log file for no data entries: same as this script.
log_file_path = os.path.join(os.path.dirname(__file__), 'no_data_from_fca2.txt')


def missing_gene_files(directory, gene_ids):
    files_in_dir = os.listdir(directory)
    present_ids = set()
    
    for file_name in files_in_dir:
        if file_name.startswith("FBgn") and file_name.endswith(".tsv"):
            gene_id = file_name.split('_')[0]  # Extract the identifier
            present_ids.add(gene_id)
    
    missing_ids = [gene_id for gene_id in gene_ids if gene_id not in present_ids]
    
    return missing_ids



def filter_fbgn_list_by_log(fbgn_list, log_file):
    """
    Remove FBgns from fbgn_list if they occur in the log_file for all gene_type.
    """
    # Initialize a dictionary to count the occurrences of each fbgn for each gene_type
    fbgn_gene_type_count = {fbgn: set() for fbgn in fbgn_list}
    
    # Open the log_file and check the occurrences of each fbgn for each gene_type
    with open(log_file, 'r') as f:
        for line in f:
            # Each line of the log_file should have the format 'fbgn_id gene_type'
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            fbgn_id, gene_type = parts
            if fbgn_id in fbgn_gene_type_count:
                fbgn_gene_type_count[fbgn_id].add(gene_type)

    # Keep only the fbgns that do not have all gene types ('gene', 'transcriptGene', 'mir', 'transcriptMir')
    fbgn_list_filtered = [
        fbgn for fbgn in fbgn_list
        if fbgn_gene_type_count[fbgn] != {'gene', 'transcriptGene', 'mir', 'transcriptMir'}
    ]
    

    return fbgn_list_filtered

'''
    USE THIS ONLY ONCE!
'''
def download_all_fbgn_data():
    query_reg = "SELECT uniquename FROM feature WHERE uniquename ~ '^FBgn[0-9]{7}$' AND is_obsolete = 'f';"
    query_irreg = "SELECT uniquename FROM feature WHERE uniquename ~ '^FBgn[0-9]{7}:1$' AND is_obsolete = 'f';"

    # Establish a connection to the FlyBase database
    conn = psycopg2.connect(
        host="chado.flybase.org",
        database="flybase",
        user="flybase"
    )    
    cur = conn.cursor()
    cur.execute(query_reg)
    results = cur.fetchall() 
    print(len(results))
    fbgn_set = set()
    for gene_data in results:
        fbgn_set.add(gene_data[0])

    # cur.execute(query_irreg)
    # results = cur.fetchall() 
    # print(len(results))
    # print(len(fbgn_set))
    # for gene_data in results:
    #     fbgn_set.add(gene_data[0].split(':')[0])  # gets the FBgn# without ':1'
    # exit(9)
    cur.close()
    conn.close()
    #print(f'Trying to get {len(fbgn_set)} files from FCA2\n{fbgn_set}\n{len(fbgn_set)}')
    # exit(9)

            
    symbols_to_fbgn_table = FlybasePrecomputedTable(symbols_to_fbgn_file)
    dframe = symbols_to_fbgn_table.to_pandas_dataframe()
    fbgn_list = dframe['gene_ID'].tolist()
    print(len(fbgn_list))
    fbgn_list = dframe['gene_ID'].unique().tolist()
    print(len(fbgn_list))
    fbgn_set.update(fbgn_list)
    print(len(fbgn_set))
    fbgn_list = list(fbgn_set)

    # Call the function that checks the log_file and filters the fbgn_list
    fbgn_list = filter_fbgn_list_by_log(fbgn_list, log_file)

    fbgn_list = missing_gene_files(out_dir, fbgn_list)
    last = len(fbgn_list)
    i = last - 1
    for _ in fbgn_list:
        download_rnaseq_data(fbgn_list[i])
        print(f'{i} files to download...')
        i -= 1

    # #for gene in net_act_expanded_genes_list:
    # for gene in fbgn_list:
    #     if gene.startswith("FBgn"):
    #         download_rnaseq_data(gene)




# download_all_fbgn_data()
# exit(9)




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




def parse_file(file_path, micro_rna = False):
    # Lista que armazenará os dados para o dataframe final
    data = []

    # Abre o arquivo e lê o conteúdo linha por linha
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Extrai o FBgn ID da linha correspondente
    fbgn_id = lines[0].split("\t")[1].strip()


    # Início da tabela de dados (após a linha do cabeçalho)
    start_index = 7

    # Faz um loop pelas linhas de dados da tabela
    for line in lines[start_index:]:
        # Separa os valores da linha atual
        values = line.strip().split("\t")
        
        # Ignora linhas que não têm dados válidos
        if len(values) < 3:
            continue

        # Pega os dados por condição de "Adult Male", "Adult Female" e, se houver, "Larval"
        tissue = values[0]
        #adult_male_data = values[1:4]
        adult_male_data = [ value for value in values[1:4] if value != '-' ]
        #adult_female_data = values[4:7]
        adult_female_data = [ value for value in values[4:7] if value != '-' ]
        larval_data = [ value for value in values[9:12] if value != '-' ]

        # Adiciona as linhas formatadas para Adult Male
        #if all(adult_male_data):
        if adult_male_data != []:
            data.append([fbgn_id, 'Adult Male', tissue] + adult_male_data)

        # Adiciona as linhas formatadas para Adult Female
        #if all(adult_female_data):
        if adult_female_data != []:
            data.append([fbgn_id, 'Adult Female', tissue] + adult_female_data)

        # Adiciona as linhas formatadas para Larval, se houver
        #if all(larval_data):
        if larval_data != []:
            data.append([fbgn_id, 'Larval', tissue] + larval_data)

    # Retorna um dataframe com os dados processados
    if micro_rna:
        return pd.DataFrame(data, columns=['#FBgene ID', 'Tissue stage and sex', 'Tissue', 'TPM', 'SD', 'Enrichment'])
    else:
        return pd.DataFrame(data, columns=['#FBgene ID', 'Tissue stage and sex', 'Tissue', 'FPKM', 'SD', 'Enrichment'])


def generate_fca2_fbgn_output(directory, output_file="fca2_fbgn_gene_output.tsv", micro_rna = False):
    if micro_rna:
        files = [f for f in os.listdir(directory) if re.match(r"^FBgn\d{7}_mir\.tsv$", f)]    
    else:
        # Usa regex para garantir que o nome siga o padrão FBgnXXXXXXXX_gene.tsv
        files = [f for f in os.listdir(directory) if re.match(r"^FBgn\d{7}_gene\.tsv$", f)]

    # Lista para armazenar todos os dataframes
    all_data = []

    # Lê e processa cada arquivo
    for file in files:
        file_path = os.path.join(directory, file)
        df = parse_file(file_path, micro_rna=micro_rna)  # Usa a função parse_file para processar cada arquivo
        all_data.append(df)

    # Combina todos os dataframes em um único
    combined_df = pd.concat(all_data, ignore_index=True)

    # Salva o dataframe combinado no arquivo de saída
    combined_df.to_csv(output_file, sep="\t", index=False)
    # Write the final dataframe to a gzipped TSV file
    with gzip.open(output_file, 'wt') as f:
        combined_df.to_csv(f, sep='\t', index=False)

    print(f"Arquivo {output_file} gerado com sucesso.")



# def generate_fca2_fbgn_output(directory, output_file="fca2_fbgn_gene_output.tsv"):
#     # Use a regex pattern to ensure that the file name follows the format FBgnXXXXXXXX_gene.tsv
#     files = [f for f in os.listdir(directory) if re.match(r"^FBgn\d{7}_gene\.tsv$", f)]
    
#     # List to store all the dataframes
#     all_data = []
    
#     for file in files:
#         file_path = os.path.join(directory, file)
#         # Read each file and append it to the dataframe list
#         df = pd.read_csv(file_path, sep="\t")
#         all_data.append(df)
    
#     # Combine all the dataframes into a single one
#     combined_df = pd.concat(all_data, ignore_index=True)
    
#     # Set the column names to the desired format
#     combined_df.columns = ['#FBgene ID', 'Tissue stage and sex', 'Tissue', 'FPKM', 'SD', 'Enrichment']
    
#     # Save the combined dataframe to the output file
#     combined_df.to_csv(output_file, sep="\t", index=False)

#     print(f"File {output_file} successfully generated.")

directory = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/fca2/genes_data"
out_file = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/fca2/fca2_fbgn_gene_output.tsv.gz"
generate_fca2_fbgn_output(directory, out_file=out_file)
out_file = "/mnt/hdd_2/saulo/snet/rejuve.bio/das/shared_rep/data/input/full/fca2/fca2_fbgn_mir_gene_output.tsv.gz"
generate_fca2_fbgn_output(directory, output_file = out_file, micro_rna = True)
exit(9)

#directory = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/genes"
#generate_fca2_fbgn_output(directory, output_file="/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_gene_output.tsv.gz")
#generate_fca2_fbgn_output(directory, output_file="/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_gene_output.tsv.gz")
directory = "/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/micro_rna_genes"
#generate_fca2_fbgn_output(directory, output_file="/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_genes_v2.tsv.GZ")
generate_fca2_fbgn_output(directory, output_file="/home/saulo/snet/hyperon/github/das-pk/shared_hsa_dmel2metta/data/toy/fca2/fca2_fbgn_mir_gene_output.tsv.gz", micro_rna=True)
exit(9)



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
#fca2_lib