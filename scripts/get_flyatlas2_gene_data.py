import requests


def download_rnaseq_data(fbgn_id, gene_type):
    if not (fbgn_id.startswith("FBgn") or fbgn_id.startswith("FBtr")):
        raise ValueError("FBgn_id must start with 'FBgn' or 'FBtr'")
    if gene_type not in ['gene', 'transcriptGene', 'mir', 'transcriptMir']:
        raise ValueError("GENE_TYPE must be one of 'gene', 'transcriptGene', 'mir', or 'transcriptMir'")

    url = f"https://motif.mvls.gla.ac.uk/FA2Direct/index.html?fbgn={fbgn_id}&tableOut={gene_type}"

    # Make a request to download the data
    response = requests.get(url)

    if response.status_code == 200:
        # Save the content to a file
        file_name = f"{fbgn_id}_{gene_type}.tsv"
        with open(file_name, 'w') as file:
            file.write(response.text)
        print(f"Data saved to {file_name}")
    else:
        print(f"Failed to download data: {response.status_code}")

# Example usage:
#download_rnaseq_data("FBgn0016075", "gene")
download_rnaseq_data("FBgn9999999", "gene")
