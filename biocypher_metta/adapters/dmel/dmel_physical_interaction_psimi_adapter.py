'''
1. Protein-protein: use the "representative protein" from:
    https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Unique_protein_isoforms_.28dmel_unique_protein_isoforms_fb_.2A.tsv.gz.29
    --> get its symbol nad get its FBpp from: 
    https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#FBgn_.3C.3D.3E_FBtr_.3C.3D.3E_FBpp_IDs_.28fbgn_fbtr_fbpp_.2A.tsv.29

2. Protein-RNA: same as above, but get FBtr from the same record as the representative protein ("representative transcript").
3. RNA-RNA: same as above

To speedup the processing build a dictionary mapping FBgn --> {('protein', FBpp), ()'rna', FBtr)}


ID(s) Interactor A	    0
ID(s) Interactor B	    1
Alt ID(s) Interactor A	Alt ID(s) Interactor B	Alias(es) Interactor A	Alias(es) Interactor B	
Interaction Detection Method(s)     6
Publication 1st Author(s)	
Publication ID(s)	8
Taxid Interactor A	9      # these could be of another organism (not  D. melanogaster)
Taxid Interactor B	10
Interaction Type(s) 11
Source Database(s)	Interaction Identifier(s)	Confidence Value(s)	Expansion Method(s)	Biological Role(s) Interactor A	Biological Role(s) Interactor B	Experimental Role(s) Interactor A	Experimental Role(s) Interactor B	

Type(s) Interactor A	20  # protein, ribonucleic acid
Type(s) Interactor B	21
Xref(s) Interactor A	Xref(s) Interactor B	Interaction Xref(s)	Annotation(s) Interactor A	Annotation(s) Interactor B	
Interaction Annotation(s)   27	
Host Organism(s)	Interaction Parameters	Creation Date	Update Date	Checksum Interactor A	Checksum Interactor B	Interaction Checksum	Negative	Feature(s) Interactor A	Feature(s) Interactor B	Stoichiometry Interactor A	Stoichiometry Interactor B	Identification Method(s) Participant A	Identification Method(s) Participant B			,0\\

Schemas

post translational interaction:
  is_a: expression
  inherit_properties: true
  represented_as: edge
  input_label: interacts_with
  source: protein
  target: protein
  properties:
    score: float    # NOT USED HERE!
    source_gene: gene       # A
    target_gene: gene       # B
    detection_method: mi    # Molecular Interaction Ontology
    pubmed_ids: str[]
    source_taxon_id: int
    target_taxon_id: int
    interaction_type: mi    # Molecular Interaction Ontology
    interaction_comments: str


# Schema for protein-RNA and RNA-RNA interactions
physical interaction:
  is_a: expression
  inherit_properties: true
  represented_as: edge
  input_label: physically_interacts_with
  source: [protein, transcript]
  target: [protein, transcript]
  properties:
    score: float    # NOT USED HERE!
    source_gene: gene       # A
    target_gene: gene       # B
    detection_method: mi    # Molecular Interaction Ontology
    pubmed_ids: str[]
    source_taxon_id: int
    target_taxon_id: int
    interaction_type: mi    # Molecular Interaction Ontology
    interaction_comments: str
'''

from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter
import re
import pandas


class PhysicalInteractionAdapter(Adapter):

    def __init__(self, write_properties, add_provenance, dmel_filepaths: list[str] = None):
        self.dmel_filepaths = dmel_filepaths
        self.label = 'physically_interacts_with'
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'

        super(PhysicalInteractionAdapter, self).__init__(write_properties, add_provenance)


    def get_edges(self):
        fb_psimi_table = FlybasePrecomputedTable(self.dmel_filepaths[0])
        self.version = fb_psimi_table.extract_date_string(self.dmel_filepaths[0])
        fb_psimi_table = self.pre_process(fb_psimi_table)
        print(fb_psimi_table)        
        exit(9)
        #header:
        #ID(s) Interactor A	ID(s) Interactor B	Alt ID(s) Interactor A	Alt ID(s) Interactor B	Alias(es) Interactor A	Alias(es) Interactor B	Interaction Detection Method(s)	Publication 1st Author(s)	Publication ID(s)	Taxid Interactor A	Taxid Interactor B	Interaction Type(s)	Source Database(s)	Interaction Identifier(s)	Confidence Value(s)	Expansion Method(s)	Biological Role(s) Interactor A	Biological Role(s) Interactor B	Experimental Role(s) Interactor A	Experimental Role(s) Interactor B	Type(s) Interactor A	Type(s) Interactor B	Xref(s) Interactor A	Xref(s) Interactor B	Interaction Xref(s)	Annotation(s) Interactor A	Annotation(s) Interactor B	Interaction Annotation(s)	Host Organism(s)	Interaction Parameters	Creation Date	Update Date	Checksum Interactor A	Checksum Interactor B	Interaction Checksum	Negative	Feature(s) Interactor A	Feature(s) Interactor B	Stoichiometry Interactor A	Stoichiometry Interactor B	Identification Method(s) Participant A	Identification Method(s) Participant B			,0\\
        rows = fb_psimi_table.get_rows()
        source_id = -1
        target_id = -1

        for row in rows:
            print(row)
            # exit(9)
            if "Gene symbol" in row[1]:     # to skip header (columns' names)
                continue
            id += 1
            props = {}            
            props['gene'] = row[0]
            #if row[2] != '':
            props['gene_hgnc_id'] = row[2]


    def pre_process(self, fb_psimi_table: FlybasePrecomputedTable) -> FlybasePrecomputedTable:
        df = fb_psimi_table.to_pandas_dataframe()
        '''
        ids = df['Interaction Identifier(s)'].tolist()
        print('PRIMEIRAAAAAAAAAAAAAAAAAA'+str(len(ids)))
        print(ids)

        #exit(9)
        '''
        # replace cells equal to "-" by ""
        df = df.applymap(lambda x: '' if str(x) == '-' else x)

        df["ID(s) Interactor A"] = df["ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["ID(s) Interactor B"] = df["ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        # df["Alt ID(s) Interactor A"] = df["Alt ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        # df["Alt ID(s) Interactor B"] = df["Alt ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        # df["Alias(es) Interactor A"] = df["Alias(es) Interactor A"].apply(lambda x: self._get_symbol_str(x) if x else x)
        # df["Alias(es) Interactor B"] = df["Alias(es) Interactor B"].apply(lambda x: self._get_symbol_str(x) if x else x)
        df["Publication ID(s)"] = df["Publication ID(s)"].apply(lambda x: self._get_pubmed_id(x) if x else x)
        # df["Interaction Xref(s)"] = df["Interaction Xref(s)"].apply(lambda x: self._extract_ids(x) if x else x)
        # df["Annotation(s) Interactor A"] = df["Annotation(s) Interactor A"].apply(lambda x: self._strip_comments(x) if x else x)
        # df["Annotation(s) Interactor B"] = df["Annotation(s) Interactor B"].apply(lambda x: self._strip_comments(x) if x else x)
        df["Interaction Annotation(s)"] = df["Interaction Annotation(s)"].apply(lambda x: self._strip_comments(x) if x else x)

        # These columns will be linked to MI Ontology
        mi_columns = [
            "Interaction Detection Method(s)",
            "Interaction Type(s)",
            "Source Database(s)",
            "Experimental Role(s) Interactor A",
            "Experimental Role(s) Interactor B",
            "Type(s) Interactor A",
            "Type(s) Interactor B"
        ]
        for mi_column in mi_columns:
            df[mi_column] = df[mi_column].apply(lambda x: self._get_MI_code(x) if x else x)

        df = self._remove_empty_columns(df)
        fb_psimi_table = fb_psimi_table.from_pandas_dataframe(df)

        return fb_psimi_table


# helper functions:

    def _remove_flybase_str(self, fbgn_string):
        return fbgn_string.split('|')[0].split(':')[-1]

    # Examples of alias:    flybase:ATPsynO(gene name)
    #                       flybase:"l(2)gl"(gene name)
    def _get_symbol_str(self, alias_string):
        #flybase:ATPsynO(gene name)
        #flybase:"l(2)gl"(gene name)
        return alias_string.split(':')[-1].strip('"').split('"')[0].replace("(gene name)", '')

    def _get_pubmed_id(self, db_reference_string):
        return db_reference_string.split('|')[1].split(':')[-1]

    def _get_MI_code(self, mi_string):
        try:
            parts = mi_string.split('"')
            return parts[1]
        except Exception as e:
            print(f"mi_string:: {mi_string}.  {e}")
            return ''

    def _extract_ids(self, input_string):
        ids = []
        parts = input_string.split("|")
        for part in parts:
            id_part = part.split(":")[-1]
            ids.append(id_part)
        return "|".join(ids)

    def _strip_comments(self, string):
        output_string = re.sub(r'comment:', '', string)
        output_string = re.sub(r'\|', '"|"', output_string)
        return output_string.strip('"')
    
    def _remove_empty_columns(self, df):
        # Print the names of empty columns
        empty_columns = [column for column in df.columns if df[column].dropna().empty]
        #data_columns = [column for column in df.columns if not df[column].dropna().empty]
        df = df.drop(empty_columns, axis=1)

        return df