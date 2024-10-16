from hyperon_das import DistributedAtomSpace
import time
import psycopg2
import os
import requests
import gzip
from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable

class SimpleQueries:
    def __init__(self):
        self.annotation_ids_dict = None


    def __get_annotation_id_table(self, local_dir):
        """
        Check if the remote file is already in the local directory. 
        If not, download the remote file and save it locally.
        Return a FlybasePrecomputedTable object initialized with the downloaded file.
        
        :param local_dir: Path to the local directory where the file should be saved.
        :return: FlybasePrecomputedTable object with the local file data.
        """
        # Remote file URL
        remote_url = "http://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_annotation_ID.tsv.gz"
        
        # Local file path
        local_file_name = os.path.join(local_dir, "fbgn_annotation_ID.tsv.gz")
        
        # Check if the file already exists locally
        if not os.path.exists(local_file_name):
            # If the file doesn't exist, download it
            print(f"Downloading {remote_url} ...")
            response = requests.get(remote_url, stream=True)
            if response.status_code == 200:
                with open(local_file_name, 'wb') as f:
                    f.write(response.content)
                print(f"File downloaded and saved at {local_file_name}")
            else:
                raise Exception(f"Failed to download file from {remote_url}")
        
        # Return an object of FlybasePrecomputedTable initialized with the TSV file
        return FlybasePrecomputedTable(local_file_name)


    def __build_annotation_id_dict(self):
        annotation_id_table = self.__get_annotation_id_table("../../aux_files/")
        rows = annotation_id_table.get_rows()
        annotation_id_dict = {}
        for row in rows:
            annotation_id_dict[row[4]] = (row[0], row[2])  # (symbol, fbgn)
        return annotation_id_dict


    def _get_symbol_from_flybase(self, fbgn:str) -> str:
        """
        Only for the "net_act subset of Drosophila data" because some related gene symbols are
        not stored in DAS. So get them from Flybase.

        :param fbgn: The Flybase identifier ('uniquename') of a gene.
        :return: (gene_symbol, is_obsolete) if fbgn is the current, non-obsolete uniquename; otherwise return (symbol, current_fbgn, is_obsolete)
        :raise exception: if no current, non-obsolete data is found.
        """
        query_reg = f"SELECT name, is_obsolete FROM feature WHERE uniquename = '{fbgn}';"

        # Establish a connection to the FlyBase database
        conn = psycopg2.connect(
            host="chado.flybase.org",
            database="flybase",
            user="flybase"
        )    
        cur = conn.cursor()
        cur.execute(query_reg)
        results = cur.fetchall() 
        for result in results:
            if result[1] == False:    # current symbol found
                return result
        #return None                # no 'non-obsolete' was found, try to find the fresh one:
    
        # if no 'non-obsolete' was found, try to find a fresh one:
        symbol = results[0][0]      # symbol for non-current FBgn
        query_reg = f"SELECT uniquename, is_obsolete FROM feature WHERE name = '{symbol}';"
        cur.execute(query_reg)
        results = cur.fetchall() 
        for result in results:
            if result[1] == False and str(result[0]).startswith("FB"):    # current uniquename was found for the symbol
                return (symbol, result[0], result[1])   # 

        # No current (fbgn, symbol) was found. "Maybe" the symbol is an "annotation id" that is started with "CG"
        # in Flybase
        annot_id_dict = self.__build_annotation_id_dict()
        gene_data = annot_id_dict.get(symbol)
        current_symbol = gene_data[0]
        current_fbgn = gene_data[1]
        if current_fbgn != None:
            return (current_symbol, current_fbgn, False)
        else:
            raise Exception(f"FBgn {fbgn} is obsolete. No current Fbgn for symbol {symbol}")

    """
    Retrieves a node identifier given one of its property's value.
    Important note:
    Evidently, the property value should be unique for the node type. For example, the gene symbol
    ("gene_name") property should be used as the "property_name" parameter along with the
    "gene_name" value ("property_value" parameter) to get the gene identifier. On the other side, if
    the "gene_type" property is used as the "property_name" parameter (with "protein_code"
    value for the "property value" parameter), for instance, all gene identifiers for "protein coding'
    genes would be returned. But, sometimes this is exactly what is needed: if you wish to know 
    all "protein code genes" in the knowledge base.
    
    get_FBgn(), get_groups_for_gene(), and get_FB_gene_disease_data() functions call this function.

    Sentences that could be used to call this function:

    (gene_name (gene FBgn0038655) CG14297)                      get_node_entity_id("gene_name", "gene", "CG14297")  ->  ["FBgn0038655"]
    (transcript_name (transcript FBtr0415465) Abd-B-RJ)
    (gene (dmel_disease_model dmel_disease_model_93) FBgn0283831)
    (genes (gene_group FBgg0001425) FBgn0266100)                get_node_entity_id("genes", "gene_group", "FBgn0266100") -> ["FBgg0001425"]
    """
    def get_node_entity_id(self, property_name: str, 
                                main_entity_type: str, propert_value: str) -> list[str]:
        q = {
            "atom_type": "link",
            "type": "Expression",
            "targets": [
                {"atom_type": "node", "type": "Symbol", "name": property_name},
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": main_entity_type},
                        {"atom_type": "variable", "name": "v1"},  # id
                    ]
                },
                {"atom_type": "node", "type": "Symbol", "name": propert_value},            
            ]
        }
        answer = das.query(q)
        ids = []
        for result in answer:
            ids.append(result.subgraph['targets'][1]['targets'][1]['name'])

        return ids


    """
    Gets all property's values given its name (property_name parameter), type of the entity it is
    linked to (main_entity_type parameter) along with its id (main_entity_id parameter).

    get_gene_symbol(), get_group_genes_FBgn(), get_groups_for_gene(), and get_GO_term_name() among others use this function.

    The most common use for this function is to get the name of entities given their ids:    
    
    (gene_name (gene FBgn0038655) CG14297)              get_node_property_value("gene_name", "gene", "FBgn0038655") -> ["CG14297"]
    (term_name (go GO:2001050) negative_regulation_of_tendon_cell_differentiation)
    (group_symbol (gene_group FBgg0001425) M19P)
    (group_name (gene_group FBgg0001425) M19_METALLODIPEPTIDASES)
    (gene (transcript FBtr0415465) FBgn0000015)
    (transcript_name (transcript FBtr0415465) Abd-B-RJ)

    But it is used to get all entities belonging to another "main" one:

    (genes (gene_group FBgg0001425) FBgn0266100)        get_node_property_value("genes", "gene_group", "FBgg0001425") -> ["FBgn0266100", "FBgn0261804", "FBgn0039420", "FBgn0266100", "FBgn0261804", "FBgn0039420"]

    (genes (gene_group FBgg0001425) FBgn0261804)
    (genes (gene_group FBgg0001425) FBgn0039420)
    (genes (gene_group FBgg0001425) FBgn0266100)
    (genes (gene_group FBgg0001425) FBgn0261804)
    (genes (gene_group FBgg0001425) FBgn0039420)
    """
    def get_node_property_values(self, property_name: str, 
                                main_entity_type: str, main_entity_id: str) -> list[str]:
        q = {
            "atom_type": "link",
            "type": "Expression",
            "targets": [
                {"atom_type": "node", "type": "Symbol", "name": property_name},
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": main_entity_type},
                        {"atom_type": "node", "type": "Symbol", "name": main_entity_id} 
                    ]
                },
                {"atom_type": "variable", "name": "v1"},  # symbol
            ]
        }
        answer = das.query(q)
        query_property_value = []
        for result in answer:
            query_property_value.append(result.subgraph['targets'][2]['name'])

        return query_property_value


    """
    Gets the id for a given ***edge source***.
    An edge has the following MeTTa format:

    (edge_type (source_type source_id) (target_type target_id))

    For example:

    (belongs_to (gene FBgn0000015) (go GO:0007621))

    where 
    "belongs_to" is the edge_type;
    "(gene FBgn0000015)" is the edge SOURCE: "gene" is the source_type and "FBgn0000015" is the 
    source_id that should be returned. 
    "(go GO:0007621)" is the edge target: "go" is the target_type and "GO:0007621" i the target_id.
    
    In fact, a list of all source_ids is returned for a particular target (GO term identifier in the example). 
    So, yes, it's possible to exist many genes (sources) for a specific GO term (target).

    More examples of data that could be retrieved:

    (regulates (gene FBgn0001291) (gene FBgn0010909))
    (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485))
    (paralogs_genes (gene FBgn0023076) (gene FBgn0002723))
    (transcribed_from (transcript FBtr0083695) (gene FBgn0038654)
    (interacts_with (protein CTNA) (protein MEF2))
    """
    def get_edge_sources(self, edge_type: str, source_type: str, target_type: str, target_id: str) -> list[str]:
        q = {
            "atom_type": "link",
            "type": "Expression",
            "targets": [
                {"atom_type": "node", "type": "Symbol", "name": edge_type},
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": source_type},
                        {"atom_type": "variable", "name": "v1"},
                    ]
                },
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": target_type},          
                        {"atom_type": "node", "type": "Symbol", "name": target_id},                             
                    ]
                },
                            
            ]
        }
        answer = das.query(q)
        fbgn_list = []
        for result in answer:
            fbgn_list.append(result.subgraph['targets'][1]['targets'][1]['name'])

        return fbgn_list        


    """
    Gets all these kinds of data represented in MeTTa as specified by the parameters:

    (belongs_to (gene FBgn0000015) (go GO:0007621))

    where "(go GO:0007621)" is the edge TARGET. A list of all targets are returned for a 
    particular source (gene identifier). So, yes, it's possible to exist many targets
    (GO terms ids) for a specific gene id (source).    
    
    More examples of data that could be retrieved:

    (regulates (gene FBgn0001291) (gene FBgn0010909))
    (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485))
    (paralogs_genes (gene FBgn0023076) (gene FBgn0002723))
    (transcribed_from (transcript FBtr0083695) (gene FBgn0038654)
    (interacts_with (protein CTNA) (protein MEF2))
    """
    def get_edge_targets(self, edge_type: str, source_type: str, source_id: str, target_type: str) -> list[str]:
        q = {
            "atom_type": "link",
            "type": "Expression",
            "targets": [
                {"atom_type": "node", "type": "Symbol", "name": edge_type},
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": source_type},
                        {"atom_type": "node", "type": "Symbol", "name": source_id},
                    ]
                },
                {
                    "atom_type": "link",
                    "type": "Expression",
                    "targets": [
                        {"atom_type": "node", "type": "Symbol", "name": target_type},          
                        {"atom_type": "variable", "name": "v1"},
                    ]
                },                            
            ]
        }
        answer = das.query(q)
        target_id_list = []
        for result in answer:
            target_id_list.append(result.subgraph['targets'][2]['targets'][1]['name'])

        return target_id_list        


    """
    Gets all property's values (not all properties themselves) for an given ***edge***.

    Examples for anotation GO terms:
    Edge definition ("belongs to"):
    
    (belongs_to (gene FBgn0023076) (go GO:0032922))  

    Edge properties definition ("qualifier", "db_reference", "evidence", etc):

    (qualifier (belongs_to (gene FBgn0023076) (go GO:0032922)) involved_in)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) FB:FBrf0237417)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) PMID:29174887)
    (evidence (belongs_to (gene FBgn0023076) (go GO:0032922)) IMP)
    (taxon (belongs_to (gene FBgn0023076) (go GO:0032922)) 7227)
    (source (belongs_to (gene FBgn0023076) (go GO:0032922)) GO)
    (source_url (belongs_to (gene FBgn0023076) (go GO:0032922)) https://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz)

    Examples for orthologous genes:
    
    Edge definition ("Orthologs genes"):
    
    (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485))

    Edge properties definition ("hsa_hgnc_id", "hsa_omim_id", "hsa_hgnc_symbol",etc):

    (hsa_hgnc_id (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) HGNC:7895)
    (hsa_omim_id (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) MIM:603347)
    (hsa_hgnc_symbol (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) NPAS2)
    (DIOPT_score (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) 10)
    (source_organism (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) Drosophila_melanogaster)
    (target_organism (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) Homo_sapiens)    
    """
    def get_edge_property_values(self, property_name: str, edge_name: str, 
                                 source_name: str, source_id: str, target_name: str, target_id: str) -> list[str]:
        q = {
            "atom_type": "link",
            "type": "Expression",
            "targets": [
                {"atom_type": "node", "type": "Symbol", "name": property_name},
                {
                "atom_type": "link",
                "type": "Expression",
                "targets": [
                    {"atom_type": "node", "type": "Symbol", "name": edge_name},
                    {
                        "atom_type": "link",
                        "type": "Expression",
                        "targets": [
                            {"atom_type": "node", "type": "Symbol", "name": source_name},
                            {"atom_type": "node", "type": "Symbol", "name": source_id},
                        ]
                    },
                    {
                        "atom_type": "link",
                        "type": "Expression",
                        "targets": [
                            {"atom_type": "node", "type": "Symbol", "name": target_name},
                            {"atom_type": "node", "type": "Symbol", "name": target_id},
                        ]
                    },
                    ]
                },
                {"atom_type": "variable", "name": "v1"},     
            ]
        }
        answer = das.query(q)
        target_id_list = []
        for result in answer:
            target_id_list.append(result.subgraph['targets'][2]['name'])

        return target_id_list        

    """
    Returns *** a list *** of Flybase ids. Mostly, the list contains only one element.
    """
    def get_FBgn(self, gene_symbol: str):
        gene_symbol = gene_symbol.replace('(', r'\(').replace(')', r'\)')
        fb_id = self.get_node_entity_id('gene_name', 'gene', gene_symbol)
        #assert len(fb_id) == 1, f'No FBgn for gene {gene_symbol}'
        return fb_id[0]
    

    def get_gene_symbol(self, fbgn: str):
        return self.get_node_property_values("gene_name", "gene", fbgn)


    def get_group_genes_FBgn(self, group_FBgg: str):
        return self.get_node_property_values("genes", "gene_group", group_FBgg)


    def get_groups_for_gene(self, gene_id: str):
        return self.get_node_entity_id('genes', 'gene_group', gene_id)
        
    def get_GO_term_name(self, go_id: str):
        return self.get_node_property_values("term_name", "go", go_id)

    """
    Get all GO terms  annotated to a given gene (symbol)

    (belongs_to (gene FBgn0023076) (go GO:0032922))
    (qualifier (belongs_to (gene FBgn0023076) (go GO:0032922)) involved_in)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) FB:FBrf0237417)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) PMID:29174887)
    (evidence (belongs_to (gene FBgn0023076) (go GO:0032922)) IMP)
    (taxon (belongs_to (gene FBgn0023076) (go GO:0032922)) 7227)
    (source (belongs_to (gene FBgn0023076) (go GO:0032922)) GO)
    (source_url (belongs_to (gene FBgn0023076) (go GO:0032922)) https://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz)
    """
    def get_GO_terms_for_gene(self, fbgn: str):
        return self.get_edge_targets('belongs_to', 'gene', fbgn, 'go')


    def get_genes_for_GO_term(self, go_term_id: str):
        return self.get_edge_sources('belongs_to', 'gene', 'go', go_term_id)

    """
    Get GO term name     -----> this and queries like get_gene_name can be generalized

    (go GO:0032922)
    (term_name (go GO:0032922) circadian_regulation_of_gene_expression)
    (synonyms (go GO:0032922) circadian_regulation_of_protein_expression)
    (synonyms (go GO:0032922) diurnal_variation_of_gene_expression)
    (synonyms (go GO:0032922) diurnal_variation_of_protein_expression)
    (source (go GO:0032922) Gene_Ontology)
    (source_url (go GO:0032922) http://purl.obolibrary.org/obo/go.owl)
    (subontology (go GO:0032922) biological_process)
    """
    def get_GO_term_name(self, go_id: str):
        return self.get_node_property_values("term_name", "go", go_id)

    
    """
    Gets all paralogs of a gene.    **** IT DOESN'T CHECK DIOPT_SCORE YET...

    (paralogs_genes (gene FBgn0023076) (gene FBgn0002723))
    (source_symbol (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) Clk)
    (target_symbol (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) Met)
    (DIOPT_score (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) 3)
    (taxon_id (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) 7227)
    """
    '''
    Gets all paralogs of a gene.    **** IT DOESN'T CHECK DIOPT_SCORE YET...

    (paralogs_genes (gene FBgn0023076) (gene FBgn0002723))
    (source_symbol (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) Clk)
    (target_symbol (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) Met)
    (DIOPT_score (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) 3)
    (taxon_id (paralogs_genes (gene FBgn0023076) (gene FBgn0002723)) 7227)
    '''
    def get_FB_paralogs(self, fbgn: str):
        paralogs_ids = self.get_edge_targets('paralogs_genes', 'gene', fbgn, 'gene')
        paralogs_symbols = []
        for paralog_id in paralogs_ids:
            gene_symbol = self.get_edge_property_values("target_symbol", "paralogs_genes", "gene", fbgn, "gene", paralog_id)
            if gene_symbol == []:
                paralogs_symbols.extend([None])
            else:
                paralogs_symbols.extend(gene_symbol)
        
        return list( zip(paralogs_ids, paralogs_symbols) )


    """
    Get all genes annotated to a given GO term

    (belongs_to (gene FBgn0023076) (go GO:0032922))
    (qualifier (belongs_to (gene FBgn0023076) (go GO:0032922)) involved_in)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) FB:FBrf0237417)
    (db_reference (belongs_to (gene FBgn0023076) (go GO:0032922)) PMID:29174887)
    (evidence (belongs_to (gene FBgn0023076) (go GO:0032922)) IMP)
    (taxon (belongs_to (gene FBgn0023076) (go GO:0032922)) 7227)
    (source (belongs_to (gene FBgn0023076) (go GO:0032922)) GO)
    (source_url (belongs_to (gene FBgn0023076) (go GO:0032922)) https://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz)
    """
    def get_genes_for_GO_term(self, go_term_id: str):
        return self.get_edge_sources('belongs_to', 'gene', 'go', go_term_id)
    

    """
    Gets orthologs of a Fruit Fly gene (given by its FBgn). **** IT DOESN'T CHECK DIOPT_SCORE YET...

    (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485))
    (hsa_hgnc_id (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) HGNC:7895)
    (hsa_omim_id (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) MIM:603347)
    (hsa_hgnc_symbol (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) NPAS2)
    (DIOPT_score (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) 10)
    (source_organism (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) Drosophila_melanogaster)
    (target_organism (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) Homo_sapiens)
    """
    def get_human_FB_orthologs(self, fbgn: str):
        return self.get_edge_targets('orthologs_genes', 'gene', fbgn, 'gene')


    """
    (hsa_hgnc_id (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) HGNC:7895)
    (hsa_hgnc_symbol (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) NPAS2)
    (DIOPT_score (orthologs_genes (gene FBgn0023076) (gene ENSG00000170485)) 10)
    
    """
    def get_gene_orthologs_hgnc_data(self, gene_symbol: str) -> dict:
        fbgn = self.get_FBgn(gene_symbol)
        orthologs_ids = self.get_human_FB_orthologs(fbgn)   # list of Ensembl ids for the gene's  orthologous
        orthologs_hgnc_ids = {}
        for ortholog_id in orthologs_ids:
            orthologs_hgnc_ids[ortholog_id] = {
                "hgnc_id": self.get_edge_property_values("hsa_hgnc_id", "orthologs_genes", "gene", fbgn, "gene", ortholog_id)[0],
                "hgnc_symbol": self.get_edge_property_values("hsa_hgnc_symbol", "orthologs_genes", "gene", fbgn, "gene", ortholog_id)[0],
                "DIOPT_score": self.get_edge_property_values("DIOPT_score", "orthologs_genes", "gene", fbgn, "gene", ortholog_id)[0],
            }
        return orthologs_hgnc_ids


    """
    Retrieves Flybase disease model data.

    MeTTa representation example: in order to keep data semantics an identifiers for each record was
    created.

    (dmel_disease_model dmel_disease_model_93)
    (taxon_id (dmel_disease_model dmel_disease_model_93) 7227)
    (gene (dmel_disease_model dmel_disease_model_93) FBgn0283831)
    (gene_hgnc_id (dmel_disease_model dmel_disease_model_93) HGNC:15924)
    (do_qualifier (dmel_disease_model dmel_disease_model_93) model_of)
    (do_term_id (dmel_disease_model dmel_disease_model_93) DOID:305)
    (do_term_name (dmel_disease_model dmel_disease_model_93) carcinoma)
    (allele (dmel_disease_model dmel_disease_model_93) FBal0320822)
    (interacting_alleles (dmel_disease_model dmel_disease_model_93) FBal0093088)
    (ev_code_interact_alleles (dmel_disease_model dmel_disease_model_93) is_ameliorated_by_FLYBASE:Myc[UAS.cZa];_FB:FBal0093088)
    (reference_id (dmel_disease_model dmel_disease_model_93) FBrf0245231)
    """
    def get_FB_gene_disease_data(self, gene_symbol: str) -> dict[str, dict[str, list[str]]]:
        fbgn = self.get_FBgn(gene_symbol)
        # gets ALL disease model ids related to the gene
        model_ids = self.get_node_entity_id("gene", "dmel_disease_model", fbgn)
        props_names = ['gene_hgnc_id', 'do_qualifier', 'do_term_id', 'do_term_name',
                        'allele', 'interacting_alleles', 'ev_code_interact_alleles' ]
        return self._get_entities_node_data(props_names, "dmel_disease_model", model_ids, fbgn)
    
        # models_data = {}
        # props_names = ['gene_hgnc_id', 'do_qualifier', 'do_term_id', 'do_term_name',
        #                'allele', 'interacting_alleles', 'ev_code_interact_alleles' ]
        # for model_id in model_ids:
        #     models_data[model_id] = {}
        #     models_data[model_id]['gene'] = fbgn
        #     for prop_name in props_names:
        #         prop_value = self.get_node_property_value(prop_name, "dmel_disease_model", model_id)
        #         if prop_value != []:
        #             models_data[model_id][prop_name] = prop_value

        # return models_data

    def _get_entities_node_data(self, properties_names: list[str], entity_type: str, entities_ids: list[str]):
        entity_data = {}        
        for entity_id in entities_ids:
            entity_data[entity_id] = {}            
            for prop_name in properties_names:
                prop_value = self.get_node_property_values(prop_name, entity_type, entity_id)
                if prop_value != []:
                    entity_data[entity_id][prop_name] = prop_value

        return entity_data
    


    def _get_entities_edge_data(self, properties_names: list[str], entity_type: str, target_type, entities_target_ids: list[str], 
                                source_type:str, source_id: str):
        entity_data = {}        
        for entity_target_id in entities_target_ids:
            entity_data[entity_target_id] = {}            
            for prop_name in properties_names:
                prop_value = self.get_edge_property_values(prop_name, entity_type, target_type, entity_target_id, source_type, source_id)
                if prop_value != []:
                    entity_data[entity_target_id][prop_name] = prop_value

        return entity_data
    
    """
    Retrieves all genes regulated by a given gene refered by its  ID (FBgn, ENSG, etc)

    MeTTa data representation:

    (regulates (gene FBgn0023076) (gene FBgn0000633))
    (evidence (regulates (gene FBgn0023076) (gene FBgn0000633)) pubmed:27924024)
    (databases (regulates (gene FBgn0023076) (gene FBgn0000633)) GTRD)
    (evidence_type (regulates (gene FBgn0023076) (gene FBgn0000633)) large_scale_evidence)
    (detection_method (regulates (gene FBgn0023076) (gene FBgn0000633)) chromatin_immunoprecipitation_assay)
    (source (regulates (gene FBgn0023076) (gene FBgn0000633)) TFLink)
    (source_url (regulates (gene FBgn0023076) (gene FBgn0000633)) tflink.net)
    (taxon_id (regulates (gene FBgn0023076) (gene FBgn0000633)) 7227)
    """
    def get_regulated_gene_data(self, gene_id: str):
        regulated_ids = self.get_edge_targets("regulates", "gene", gene_id, "gene")
        properties_names = ['evidence', 'databases', 'evidence_type', 'detection_method', 'source', 'source_url', 'taxon_id']
        return self._get_entities_edge_data(properties_names, 'regulates', "gene", regulated_ids, "gene", gene_id)

    
    def _get_gene_transcripts(self, gene_id: str) -> list[str]:
        """
        Gets all transcripts of the gene 

        :param gene_id: the identifier of a gene
        :return: a list of transcripts ids
        
        MeTTa:
        (transcript FBtr0415465)
        (gene (transcript FBtr0415465) FBgn0000015)
        (transcript_name (transcript FBtr0415465) Abd-B-RJ)
        (transcript_type (transcript FBtr0415465) protein_coding)
        (chr (transcript FBtr0415465) 3R)
        (start (transcript FBtr0415465) 16927667)
        (end (transcript FBtr0415465) 16968574)
        (taxon_id (transcript FBtr0415465) 7227)
        (source (transcript FBtr0415465) GENCODE)
        (source_url (transcript FBtr0415465) https://www.gencodegenes.org/human/)
        """

    def _get_gene_proteins(self, gene_id: str) -> list[str]:
        """
        Gets all
        """

    def __get_central_dogma(gene_id): 
        """
        Retrieves the 'Ccentral Dogma of Biology': (DNA)gene -> (RNA)transcript -> protein

        1. get_gene_transcripts()
        2. get_transcript_proteins()

        
        Not a simple query...
        """

    def central_dogma(gene_symbol: str):



#######################################################################################################################
host = '127.0.0.1'
port = '8080'

print(f'Connecting to DAS...')
das = DistributedAtomSpace(query_engine='remote', host=host, port=port)
#DistributedAtomSpace().about()

print(f"Connected to DAS at {host}:{port}")
print("(nodes, links) =", das.count_atoms({'precise': True}))



querier = SimpleQueries()
print(querier._get_symbol_from_flybase('FBgn0039084'))
#exit(9)
net_act_genes_list = ["Su(var)205", "Top3beta", "Mef2", "Clk", "Dref", "TfIIB", "Myc", "AGO2", "Nipped-B", "Cp190", "TfIIA-L",
               "Trl", "ash1", "Raf", "Abd-B", "Orc2", "Rbf", "mof", "msl-1", "Hmr"]

four_tfs_symbols_list = ["Abd-B", "Clk", "Mef2", "Myc"]#["Myc"]#["Abd-B", "Clk", "Mef2", "Myc"]

gene_symbols_list = net_act_genes_list
gene_symbols_list = four_tfs_symbols_list
gene_ids_list = ['FBgn0000015', 'FBgn0023076', 'FBgn0011656', 'FBgn0262656']#['FBgn0262656'] # to speed-up dev tests#['FBgn0000015', 'FBgn0023076', 'FBgn0011656', 'FBgn0262656'] # to speed-up dev tests

print("\n\nStarted getting gene data...")

# print(f'Retrieving gene ids...')
# gene_ids_list = []
# for gene_symbol in gene_symbols_list:    
#     start = time.time()
#     print(f'Retrieving ID for gene {gene_symbol}...')
#     fbgn = querier.get_FBgn(gene_symbol)
#     gene_ids_list.append( fbgn )
#     print(f"FBgn for\t{gene_symbol}:\t{fbgn}" )
#     finish = time.time()
#     print(f'Answer time: {finish - start} sec.')
print(f'Regulatory relationships...')
for gene_symbol, gene_id in zip(gene_symbols_list, gene_ids_list):
    print(f'{gene_symbol} regulates:')
    reg_data = querier.get_regulated_gene_data(gene_id)
    for target_id in reg_data.keys():
        keys = [value for value in reg_data[target_id]]
        out_str = ""
        for key in keys:
            if isinstance(reg_data[target_id][key], str):
                out_str += reg_data[target_id][key] + '\t'
            else:                                                       # list
                for value in reg_data[target_id][key]:
                    out_str += value + '\t'        
        if out_str == "":
            gsd = querier._get_symbol_from_flybase(target_id)
            if gsd[1]:
                target_symbol = gsd     # OBSOLETE ID!
                print(f'{target_id} {str(target_symbol)}: OBSOLETE ID FOUND') 
                """
                Abd-B
                FBgn0034577 ('cpa', True): OBSOLETE ID FOUND
                FBgn0000473 ('Cyp6a2', True): OBSOLETE ID FOUND
                Myc
                FBgn0025628 ('CG4199', True): OBSOLETE ID FOUND
                FBgn0001226 ('Hsp27', True): OBSOLETE ID FOUND
                FBgn0261800 ('LanB1', True): OBSOLETE ID FOUND
                FBgn0039084 ('CG10175', True): OBSOLETE ID FOUND
                FBgn0036762 ('CG7430', True): OBSOLETE ID FOUND

                FBgn0000473 ('Cyp6a2', True): OBSOLETE ID FOUND
                FBgn0036078 ('Or67c', True): OBSOLETE ID FOUND
                """
            else:
                target_symbol = gsd[0]
                print(f'{target_id} {str(target_symbol)}: Regulatory data not in DAS') 
        else:
            target_symbol = querier.get_gene_symbol(target_id)
            if target_symbol == []:
                target_symbol = str(querier._get_symbol_from_flybase(target_id))
            print(f'{target_id} {target_symbol}: {out_str}')

#exit(9)
print(f'\nFlybase disease model data...')
#for gene_symbol in gene_symbols_list:    
for gene_symbol, gene_id in zip(gene_symbols_list, gene_ids_list):
    disease_data = querier.get_FB_gene_disease_data(gene_symbol)
    print(f'\nDisease models for {gene_symbol}:')
    for disease_model_id in disease_data.keys():
        keys = [value for value in disease_data[disease_model_id]]
        out_str = ""
        for key in keys:
            if isinstance(disease_data[disease_model_id][key], str):
                out_str += disease_data[disease_model_id][key] + '\t'
            else:                                                       # list
                for value in disease_data[disease_model_id][key]:
                    out_str += value + '\t'
        print(f'{disease_model_id}: {out_str}')


print(f'Retrieving orthologous HGNC data...')
#four_tfs_ids = []
for gene_symbol in gene_symbols_list:    
    print(f'Orthologs HGNC data for {gene_symbol}')
    start = time.time()
    hgnc_dict = querier.get_gene_orthologs_hgnc_data(gene_symbol)
    print(f'Ensembl ID\tHGNC ID\tHGNC Symbol\tDIOPT score')
    for ensembl_id in hgnc_dict.keys():
        print(f'{ensembl_id}\t{hgnc_dict[ensembl_id]["hgnc_id"]}\t{hgnc_dict[ensembl_id]["hgnc_symbol"]}\t{hgnc_dict[ensembl_id]["DIOPT_score"]}')
    #print(querier.get_gene_orthologs_hgnc_data(gene_symbol))
    finish = time.time()
    print(f'Answer time: {finish - start} sec.')


print('\n\nParalogs...\n"')
#exit(9)
for gene_id, gene_symbol in zip(gene_ids_list, gene_symbols_list):
    start = time.time()
    lps = querier.get_FB_paralogs(gene_id)
    print(len(lps))
    genes_data = set(lps)
    print(len(gene_data))
    #genes_data = set(querier.get_FB_paralogs(gene_id))
    for gene_data in genes_data:
        if gene_data[1] != None:
            print(f'Paralog for {gene_symbol}:\t{gene_data[0]}\t{gene_data[1]}')
        else:
            print(f'Paralog for {gene_symbol}:\t{gene_data[0]}')
    finish = time.time()
    print(f'Answer time: {finish - start} sec.')


print('\n\nGroups...\n"')
for fbgn, symbol in zip(gene_ids_list, gene_symbols_list):
    start = time.time()
    groups = querier.get_groups_for_gene(fbgn)
    print(f'\n\nGroups for\t{symbol}:')    
    for grp in groups:
        query_property_name = "group_symbol"            # (group_symbol (gene_group FBgg0001558) VHA-V1-CAT)
        main_entity_type = "gene_group"
        group_symbol = querier.get_node_property_value(query_property_name, main_entity_type, grp)
        query_property_name = "group_name"
        group_name = querier.get_node_property_value(query_property_name, main_entity_type, grp)[0]
        print(f"Group symbol for\t{grp}:\t{group_symbol}")
        print(f"Group name for\t{grp}:\t{group_name.replace('_', ' ')}")
        grp_fbgns = querier.get_group_genes_FBgn(grp)
        for bud_fbgn in grp_fbgns:
            bud_symbol = querier.get_gene_symbol(bud_fbgn)
            if bud_symbol == []:
                #print(f'No gene symbol for {b_fbgn} in the Net_act DAS. Retrieving it from Flybase...')
                bud_symbol = querier._get_symbol_from_flybase(bud_fbgn)
            else:
                bud_symbol = bud_symbol[0]
            print(f'{symbol} buddies in group\t{grp}:\t{bud_fbgn}\t<-->\t{str(bud_symbol).replace(",", "")}')
    finish = time.time()
    print(f'Answer time (gene {symbol}): {finish - start} sec.')


print(f'\n\nGene ontology terms...')
for gene_id, gene_symbol in zip(gene_ids_list, gene_symbols_list):
    print(f'\nStarting getting GO terms  genes related to\t{gene_symbol}...')
    # start = time.time()
    go_ids = querier.get_GO_terms_for_gene(gene_id)
    for go_id in go_ids:
        start = time.time()
        go_term_name = querier.get_GO_term_name(go_id)[0].replace("_", " ")
        print(f'\nTerm for\t{gene_symbol}:\t{go_term_name}')
        fbgns = querier.get_genes_for_GO_term(go_id)
        print(f'Genes for GO term:\t{go_id}:')
        for fbgn in fbgns:
            g_symbol = querier.get_gene_symbol(fbgn)
            if g_symbol != []:
                print(f'\t\t\t{g_symbol[0]}')
            else:
                g_symbol = querier._get_symbol_from_flybase(fbgn)
                print(f'\t\t\t{str(g_symbol).replace(",", "")}')
        finish = time.time()
        print(f'Answer time: {finish - start} sec.')



# print('\n\nOrthologous...\n')
# # print(four_tfs_ids)
# # print(four_tfs_symbols)
# for gene_id, gene_symbol in zip(four_tfs_ids, four_tfs_symbols):
#     start = time.time()
#     ensembl_ids = set( querier.get_human_FB_orthologs(gene_id))
#     for ensembl_id in ensembl_ids:
#         print(f'Human ortholog for\t{gene_symbol}:\t{ensembl_id}')
#     finish = time.time()
#     print(f'Answer time: {finish - start} sec.')


