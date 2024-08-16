'''

expression value:
  is_a: entity
  represented_as: edge
  inherit_properties: false
  input_label: expression_value
  source: [gene, transcript, exon]
  target: RNASeq library
  description: >-
    A set of values for numeric estimation of single gene/transcript expression as measured by a given methodology (represented by the library).
  properties:
    value_and_description: list[tuple]]
    taxon_id: int                        # 7227 for dmel / 9606 for hsa
    source: str
    source_url: str


FB  data:
https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#High-Throughput_Gene_Expression_.28high-throughput_gene_expression_fb_.2A.tsv.gz.29
https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Single_Cell_RNA-Seq_Gene_Expression_.28scRNA-Seq_gene_expression_fb_.2A.tsv.gz.29
https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#RNA-Seq_RPKM_values_.28gene_rpkm_report_fb_.2A.tsv.gz.29

FB scRNASeq table:
#Pub_ID	Pub_miniref	Clustering_Analysis_ID	Clustering_Analysis_Name	Source_Tissue_Sex	Source_Tissue_Stage	Source_Tissue_Anatomy	Cluster_ID	Cluster_Name	Cluster_Cell_Type_ID	Cluster_Cell_Type_Name	Gene_ID	Gene_Symbol	Mean_Expression	Spread
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0031081	Nep3	1022.4949	0.00016897600540723216
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0031088	CG15322	269.05170000000004	0.0005069280162216965
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0053217	CG33217	439.74384371428573	0.026022304832713755
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0052350	Vps11	585.499525895105	0.024163568773234202
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0024733	RpL10	3497.867660248448	0.9793849273403177
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0040372	G9a	602.1811133469388	0.05795876985468063
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0000316	cin	582.4078043088889	0.03801960121662724
FBrf0245988	Cattenoz et al., 2020, EMBO J. 39(12): e104486	FBlc0003731	scRNAseq_2020_Cattenoz_NI_seq_clustering		larval stage	embryonic/larval hemolymph	FBlc0003732	scRNAseq_2020_Cattenoz_NI_seq_clustering_plasmatocytes	FBbt:00001685	embryonic/larval plasmatocyte	FBgn0024989	CG3777	354.646665	0.0003379520108144643


FB high-throughput table
#High_Throughput_Expression_Section	Dataset_ID	Dataset_Name	Sample_ID	Sample_Name	Gene_ID	Gene_Symbol	Expression_Unit	Expression_Value
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003619	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Brain	FBgn0000003	7SLRNA:CR32864	FPKM	202396
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003619	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Brain	FBgn0289744	EndoG	FPKM	12
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003619	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Brain	FBgn0289745	Est17	FPKM	0
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003619	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Brain	FBgn0289868	Or67c	FPKM	0
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003620	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Crop	FBgn0000003	7SLRNA:CR32864	FPKM	240220
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003620	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Crop	FBgn0000008	a	FPKM	2
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003620	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Crop	FBgn0000014	abd-A	FPKM	0
FlyAtlas2 Anatomy RNA-Seq	FBlc0003498	FlyAtlas2	FBlc0003620	RNA-Seq_Profile_FlyAtlas2_Adult_Female_Crop	FBgn0000015	Abd-B	FPKM	0

FB RPKM report table:
# Release_ID	FBgn#	GeneSymbol	Parent_library_FBlc#	Parent_library_name	RNASource_FBlc#	RNASource_name	RPKM_value	Bin_value	Unique_exon_base_count	Total_exon_base_count	Count_used
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000060	BCM_1_RNAseq	FBlc0000068	BCM_1_FA3d	85	5	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000060	BCM_1_RNAseq	FBlc0000069	BCM_1_MA3d	119	6	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000060	BCM_1_RNAseq	FBlc0000070	BCM_1_P	53	5	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000060	BCM_1_RNAseq	FBlc0000071	BCM_1_L	55	5	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000060	BCM_1_RNAseq	FBlc0000072	BCM_1_A17d	37	4	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000085	modENCODE_mRNA-Seq_development	FBlc0000086	mE_mRNA_em0-2hr	2	1	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000085	modENCODE_mRNA-Seq_development	FBlc0000087	mE_mRNA_em2-4hr	1	1	299	299	Unique
Dmel_R6.58	FBgn0000003	7SLRNA:CR32864	FBlc0000085	modENCODE_mRNA-Seq_development	FBlc0000088	mE_mRNA_em4-6hr	1	1	299	299	Unique

'''

from biocypher_metta.adapters.dmel.flybase_tsv_reader import FlybasePrecomputedTable
from biocypher_metta.adapters import Adapter
from biocypher._logger import logger


class ExpressionValueAdapter(Adapter):
    def __init__(self, write_properties, add_provenance, dmel_data_filepaths):
        self.dmel_data_filepaths = dmel_data_filepaths
        self.label = 'expression_value'
        self.source = 'FLYBASE'
        self.source_url = 'https://flybase.org/'

        super(ExpressionValueAdapter, self).__init__(write_properties, add_provenance)


    def get_edges(self):
        for dmel_data_filepath in self.dmel_data_filepaths:
            if "scRNA-Seq_gene_expression_fb" in dmel_data_filepath:
                rna_table= FlybasePrecomputedTable(dmel_data_filepath)
                self.version = rna_table.extract_date_string(dmel_data_filepath)
                # header:
                # Pub_ID	Pub_miniref	Clustering_Analysis_ID	Clustering_Analysis_Name	Source_Tissue_Sex	Source_Tissue_Stage
                # Source_Tissue_Anatomy	Cluster_ID	Cluster_Name	Cluster_Cell_Type_ID	Cluster_Cell_Type_Name	Gene_ID	Gene_Symbol	Mean_Expression	Spread
                rows = rna_table.get_rows()
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
                    props['taxon_id'] = 7227
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                    yield _source, _target, self.label, props
            elif "high-throughput_gene_expression_fb" in dmel_data_filepath:
                rna_table = FlybasePrecomputedTable(dmel_data_filepath)
                self.version = rna_table.extract_date_string(dmel_data_filepath)
                # header:
                # High_Throughput_Expression_Section	Dataset_ID	Dataset_Name	Sample_ID	Sample_Name	Gene_ID
                # Gene_Symbol	Expression_Unit	Expression_Value
                rows = rna_table.get_rows()
                for row in rows:
                    props = {}
                    # source_ids = row[1].split('|')
                    _source = row[4]  # Gene_ID   FBgn#
                    _target = row[3]  # Sample_ID  FBlc#

                    props['value_and_description'] = [
                        (float(row[8]), # Expression_Value
                         row[7]         # Expression_Unit
                         ),
                    ]
                    props['taxon_id'] = 7227
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                    yield _source, _target, self.label, props
            elif "gene_rpkm_report_fb" in dmel_data_filepath:
                rna_table = FlybasePrecomputedTable(dmel_data_filepath)
                self.version = rna_table.extract_date_string(dmel_data_filepath)
                # header:
                # Release_ID	FBgn#	GeneSymbol	Parent_library_FBlc#	Parent_library_name	RNASource_FBlc#
                # RNASource_name	RPKM_value	Bin_value	Unique_exon_base_count	Total_exon_base_count	Count_used
                rows = rna_table.get_rows()
                for row in rows:
                    props = {}
                    _source = row[1]  # FBgn#
                    _target = row[5]  # RNASource_FBlc#

                    props['value_and_description'] = [
                        (row[7], # RPKM_value
                         "RPKM",
                         row[11]        #	Count_used in the RPKM (Total or Unique)
                         ),
                        (row[8],  # Bin_value
                         "The expression bin classification of this gene in this RNA-Seq experiment, "
                         "based on RPKM value. Bins range from 1 (no/extremely low expression) to 8 (extremely high expression)"
                         ),
                        (row[9],  # Unique_exon_base_count
                         "The number of exonic bases unique to the gene (not overlapping exons of other genes). Field "
                         "will be blank for genes derived from dicistronic/polycistronic transcripts"
                         ),
                        (row[10],  # Total_exon_base_count
                         "The number of bases in all exons of this gene"
                         ),
                    ]
                    props['taxon_id'] = 7227
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url
                    yield _source, _target, self.label, props


