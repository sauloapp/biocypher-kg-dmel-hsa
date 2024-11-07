'''
1. Protein-protein: use the "representative protein" from:
    https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Unique_protein_isoforms_.28dmel_unique_protein_isoforms_fb_.2A.tsv.gz.29
    --> get its symbol nad get its FBpp from: 
    https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#FBgn_.3C.3D.3E_FBtr_.3C.3D.3E_FBpp_IDs_.28fbgn_fbtr_fbpp_.2A.tsv.29

2. Protein-RNA: same as above, but get FBtr from the same record as the representative protein ("representative transcript").
3. RNA-RNA: same as above

To speedup the processing build a dictionary mapping FBgn --> {('protein', FBpp), ('transcript', FBtr)}

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
Xref(s) Interactor A	Xref(s) Interactor B	Interaction Xref(s)	
Annotation(s) Interactor A	25
Annotation(s) Interactor B	26
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
    source_comments: str
    target_comments: str    


# Schema for protein-RNA and RNA-RNA interactions
physical interaction:
  is_a: expression
  inherit_properties: true
  represented_as: edge
  input_label: physically_interacts_with
  source: [protein, transcript]                                 Type(s) Interactor A	20  # protein, ribonucleic acid
  target: [protein, transcript]                                 Type(s) Interactor B	21
  properties:
    score: float    # NOT USED HERE!
    source_gene: gene       # A
    target_gene: gene       # B
    detection_method: mi    # Molecular Interaction Ontology    Interaction Detection Method(s)     6
    pubmed_ids: str[]                                           Publication ID(s)	8
    source_taxon_id: int                                        Taxid Interactor A	9      # these could be of another organism (not  D. melanogaster)
    target_taxon_id: int                                        Taxid Interactor B	10
    interaction_type: mi    # Molecular Interaction Ontology    Interaction Type(s) 11
    interaction_comments: str                                   Interaction Annotation(s)   27	
    source_comments: str                                        "Annotation(s) Interactor A"    25
    target_comments: str                                        "Annotation(s) Interactor B"    26    
'''


from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter
import re
import pandas


class PhysicalInteractionAdapter(Adapter):

    def __init__(self, write_properties, add_provenance, label, dmel_filepaths: list[str] = None):
        self.dmel_filepaths = dmel_filepaths
        self.label = label
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'

        self.fbgn_to_fbtr_fbpp_dict = self.__build_dictionary()

        super(PhysicalInteractionAdapter, self).__init__(write_properties, add_provenance)


    def get_edges(self):
        fb_psimi_table = FlybasePrecomputedTable(self.dmel_filepaths[0])
        self.version = fb_psimi_table.extract_date_string(self.dmel_filepaths[0])
        fb_psimi_table = self.pre_process(fb_psimi_table)
        
        #header:
        #ID(s) Interactor A	ID(s) Interactor B	Alt ID(s) Interactor A	Alt ID(s) Interactor B	Alias(es) Interactor A	Alias(es) Interactor B	Interaction Detection Method(s)	Publication 1st Author(s)	Publication ID(s)	Taxid Interactor A	Taxid Interactor B	Interaction Type(s)	Source Database(s)	Interaction Identifier(s)	Confidence Value(s)	Expansion Method(s)	Biological Role(s) Interactor A	Biological Role(s) Interactor B	Experimental Role(s) Interactor A	Experimental Role(s) Interactor B	Type(s) Interactor A	Type(s) Interactor B	Xref(s) Interactor A	Xref(s) Interactor B	Interaction Xref(s)	Annotation(s) Interactor A	Annotation(s) Interactor B	Interaction Annotation(s)	Host Organism(s)	Interaction Parameters	Creation Date	Update Date	Checksum Interactor A	Checksum Interactor B	Interaction Checksum	Negative	Feature(s) Interactor A	Feature(s) Interactor B	Stoichiometry Interactor A	Stoichiometry Interactor B	Identification Method(s) Participant A	Identification Method(s) Participant B			,0\\
        rows = fb_psimi_table.get_rows()

        for row in rows:
            source_fbgn = row[0]
            target_fbgn = row[1]
            if source_fbgn not in self.fbgn_to_fbtr_fbpp_dict or target_fbgn not in self.fbgn_to_fbtr_fbpp_dict:
                print(f'Gene {source_fbgn} or/and gene {target_fbgn} has/have no transcript/protein in Flybase')
                continue
            source, target = self.get_ids(source_fbgn, target_fbgn, row)
            # ptp
            # source: [protein, transcript]
            # target: protein
            if target[0] == 'transcript' and self.label.startswith('ptp'):
                continue
            # ptt
            # source: [protein, transcript]
            # target: transcript            
            if target[0] == 'protein' and self.label.startswith('ptt'):
                continue
            props = self.get_properties(row)            
            yield source, target, self.label, props
            

    def get_ids(self, source_fbgn, target_fbgn, data_row):
        source_type = 'protein' if data_row[20] == 'MI:0326' else 'transcript'  # 'MI:0326' means 'protein'
        target_type = 'protein' if data_row[21] == 'MI:0326' else 'transcript'

        fb_source = self.fbgn_to_fbtr_fbpp_dict[source_fbgn][source_type]
        fb_target = self.fbgn_to_fbtr_fbpp_dict[target_fbgn][target_type]

        return (source_type, fb_source), (target_type, fb_target)        

    def get_properties(self, row):
        properties = {}        
        properties['source_gene'] = row[0]
        properties['target_gene'] = row[1]
        properties['detection_method'] = row[6]
        properties['pubmed_ids'] = row[8]
        properties['source_taxon_id'] = row[9]
        properties['target_taxon_id'] = row[10]
        properties['interaction_type'] = row[11]
        properties['interaction_comments'] = row[27]
        properties['source_comments'] = row[25]
        properties['target_comments'] = row[26]

        return properties

    def __build_dictionary(self):
        table_fbgn_to_fbtr_fbpp = FlybasePrecomputedTable(self.dmel_filepaths[1])
        fbgn_to_fbtr_fbpp_dict = {}
        for row in table_fbgn_to_fbtr_fbpp.get_rows():
            if row[2] not in fbgn_to_fbtr_fbpp_dict:
                fbgn_to_fbtr_fbpp_dict[row[2]] = {
                    'transcript': row[7],  # FBtr
                    'protein': row[9] if row[9] != '' else None  # FBpp ou None
                }
        return fbgn_to_fbtr_fbpp_dict
    

    def pre_process(self, fb_psimi_table: FlybasePrecomputedTable) -> FlybasePrecomputedTable:
        df = fb_psimi_table.to_pandas_dataframe()        
        # replace cells equal to "-" by ""
        df = df.map(lambda x: '' if str(x) == '-' else x)

        df["ID(s) Interactor A"] = df["ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["ID(s) Interactor B"] = df["ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        # df["Alt ID(s) Interactor A"] = df["Alt ID(s) Interactor A"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        # df["Alt ID(s) Interactor B"] = df["Alt ID(s) Interactor B"].apply(lambda x: self._remove_flybase_str(x) if x else x)
        df["Publication ID(s)"] = df["Publication ID(s)"].apply(lambda x: self._get_pubmed_id(x) if x else x)
        # df["Interaction Xref(s)"] = df["Interaction Xref(s)"].apply(lambda x: self._extract_ids(x) if x else x)
        df["Annotation(s) Interactor A"] = df["Annotation(s) Interactor A"].apply(lambda x: self._strip_comments(x) if x else x)
        df["Annotation(s) Interactor B"] = df["Annotation(s) Interactor B"].apply(lambda x: self._strip_comments(x) if x else x)
        df["Interaction Annotation(s)"] = df["Interaction Annotation(s)"].apply(lambda x: self._strip_comments(x) if x else x)

        df["Taxid Interactor A"] = df["Taxid Interactor A"].apply(lambda x: self._get_taxon_id(x) if x else x)
        df["Taxid Interactor B"] = df["Taxid Interactor B"].apply(lambda x: self._get_taxon_id(x) if x else x)

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

        print(df)
        df = self._remove_empty_columns(df)
        print(df)
        fb_psimi_table = fb_psimi_table.from_pandas_dataframe(df)
        print(fb_psimi_table)
        return fb_psimi_table


# helper functions:

    def _remove_flybase_str(self, fbgn_string):
        return fbgn_string.split('|')[0].split(':')[-1]

    def _get_taxon_id(self, taxon_id: str):
        '''
        @param: taxon_id: like 'taxid:7227("Drosophila melanogaster")'
        '''
        return taxon_id.split('(')[0].split(':')[1]


    def _get_pubmed_id(self, db_reference_string):
        # return db_reference_string.split('|')[1].split(':')[-1]
        parts = db_reference_string.split('|')
        if len(parts) > 1:
            sub_parts = parts[1].split(':')
            return sub_parts[-1] if sub_parts else None
        return None


    def _get_MI_code(self, mi_string):
        try:
            parts = mi_string.split('"')        
            return parts[3]#parts[1]
        except Exception as e:
            print(f"mi_string:::::--> {mi_string}.  {e}")
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
        output_string = re.sub(r'\|', ' AND ', output_string)
        return output_string.strip('"')
    
    def _remove_empty_columns(self, df):
        # Print the names of empty columns
        empty_columns = [column for column in df.columns if df[column].dropna().empty]
        #data_columns = [column for column in df.columns if not df[column].dropna().empty]
        df = df.drop(empty_columns, axis=1)

        return df