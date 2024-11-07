from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter

class GeneToSequenceOntologyAdapter(Adapter):

    def __init__(self, write_properties, add_provenance, filepath=None):
        self.filepath = filepath
        self.label = 'gene_to_so'
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'
        super(GeneToSequenceOntologyAdapter, self).__init__(write_properties, add_provenance)


    def get_edges(self):
        gene_so_table = FlybasePrecomputedTable(self.filepath)
        self.version = gene_so_table.extract_date_string(self.filepath)
        #header:
        #AlleleID	AlleleSymbol	GeneID	GeneSymbol
        rows = gene_so_table.get_rows()
        for row in rows:
            props = {}
            source = row[0]
            target = row[3]
            props['taxon_id'] = 7227

            yield source, target, self.label, props
