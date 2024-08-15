'''
FB  data:
https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Single_Cell_RNA-Seq_Gene_Expression_.28scRNA-Seq_gene_expression_fb_.2A.tsv.gz.29

expression value:
  is_a: entity
  represented_as: edge
  inherit_properties: false
  input_label: expression_value
  _source: gene
  _target: scRNASeq library
  description: >-
    A single value for a numeric estimation of a gene expression as measured by a given methodology (represented by the library).
  properties:
    value: float[]
    value_description: str[]    # Acronym or term to clarify the expression value: "RPKM", "BIN VALUE", "SPREAD",
                                # "AFFYMETRIX CALL", "AFFYMETRIX p_VALUE", etc
                                # This property is mandatory: value_description[i] refers to value[i] for all i.
    taxon_id: int                        # 7227 for dmel / 9606 for hsa
    source: str
    source_url: str


FB cRNASeq table:
#Pub_ID	Pub_miniref	Clustering_Analysis_ID	Clustering_Analysis_Name	Source_Tissue_Sex	Source_Tissue_Stage	Source_Tissue_Anatomy	Cluster_ID	Cluster_Name	Cluster_Cell_Type_ID	Cluster_Cell_Type_Name	Gene_ID	Gene_Symbol	Mean_Expression	Spread
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0031081	Nep3	1022.4949	0.00016897600540723216
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0031088	CG15322	269.05170000000004	0.0005069280162216965
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0053217	CG33217	439.74384371428573	0.026022304832713755
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0052350	Vps11	585.499525895105	0.024163568773234202
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0024733	RpL10	3497.867660248448	0.9793849273403177
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0040372	G9a	602.1811133469388	0.05795876985468063
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0000316	cin	582.4078043088889	0.03801960121662724
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0024989	CG3777	354.646665	0.0003379520108144643

'''

from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter
from biocypher._logger import logger


class ExpressionValueAdapter(Adapter):

    def __init__(self, write_properties, add_provenance, dmel_data_filepath):
        self.dmel_data_filepath = dmel_data_filepath
        self.label = 'expression_value'
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'

        super(ExpressionValueAdapter, self).__init__(write_properties, add_provenance)


    def get_edges(self):
        scRNA_table= FlybasePrecomputedTable(self.dmel_data_filepath)
        self.version = scRNA_table.extract_date_string(self.dmel_data_filepath)
        # header:
        # Pub_ID	Pub_miniref	Clustering_Analysis_ID	Clustering_Analysis_Name	Source_Tissue_Sex	Source_Tissue_Stage
        # Source_Tissue_Anatomy	Cluster_ID	Cluster_Name	Cluster_Cell_Type_ID	Cluster_Cell_Type_Name	Gene_ID	Gene_Symbol	Mean_Expression	Spread
        rows = scRNA_table.get_rows()
        for row in rows:
            props = {}
            #source_ids = row[1].split('|')
            _source = row[11]    # gene FBgn#
            _target = row[7]     # library FBlc = Cluster ID#

            props['value_and_description'] = [
                ( float(row[13]),        # Mean_Expression
                  "Mean_Expression: is the average level of expression of the gene across all "
                                          "cells of the cluster in which the gene is detected at all"
                ),
                ( float(row[14]),    # Spread
                  "Spread: is the proportion of cells in the cluster in which the gene is detected"
                )
            ]
            # props["value"] = [float(row[13]), float(row[14])]
            # props["value_description"] = ["Mean_Expression: is the average level of expression of the gene across all "
            #                               "cells of the cluster in which the gene is detected at all",
            #                               "Spread: is the proportion of cells in the cluster in which the gene is detected"]
            props['taxon_id'] = 7227
            if self.add_provenance:
                props['source'] = self.source
                props['source_url'] = self.source_url

            yield _source, _target, self.label, props

