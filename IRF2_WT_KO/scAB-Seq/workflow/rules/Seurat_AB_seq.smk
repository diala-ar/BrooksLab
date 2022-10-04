#from snakemake.utils import R
import numpy as np

rule demultiplex_samples:
    input:
        ncounts_file="resources/seven_bridges/Combined_LUK15635_Tumor_and_draining_lymph_node_WT_vs_IRF2KO_RSEC_MolsPerCell.csv",
        tags_file="resources/seven_bridges/LUK15635_Tumor_and_draining_lymph_node_WT_vs_IRF2KO_Sample_Tag_Calls.csv",
    output:
        tagged_raw_counts="results/Seurat_AB_seq/Tum/tagged_cells_raw_counts.csv",
    params:
        samples_metadata=config['Seurat_AB_seq']['samples_metadata']
    conda:
        '../envs/Seurat_AB_seq.yaml'
    script:
        "../scripts/demultiplex_samples.R"


# perform subsequent analyses on Tumor cells only
rule seur_create_multimodal_seurat_object:
    input:
        raw_counts="results/Seurat_AB_seq/tagged_cells_raw_counts.csv"
    output:
        seur="results/Seurat_AB_seq/Seur_Tum_before_filtering.rds"
    params:
        conds=config['Seurat_AB_seq']['conds']
    conda:
        '../envs/Seurat_AB_seq.yaml'
    script:
        '../scripts/seur_create_multimodal_seurat_object.R'


rule seur_quality_control_plots:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_{before_or_after}_filtering.rds"
    output:
        hist="results/Seurat_AB_seq/QC_Tum_{before_or_after}_cell_filtering_hist.pdf",
        violon="results/Seurat_AB_seq/QC_Tum_{before_or_after}_cell_filtering_violon.pdf"
    params:
        conds=config['Seurat_AB_seq']['conds']
    conda:
        '../envs/Seurat_AB_seq.yaml'
    script:
        '../scripts/seur_quality_control_plots.R'


rule seur_filter_cells_n_genes:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_before_filtering.rds"
    output:
        seur="results/Seurat_AB_seq/Seur_Tum_after_filtering.rds"
    params:
        conds=config['Seurat_AB_seq']['conds'],
        criteria_cutoffs_file=config['Seurat_AB_seq']['criteria_cutoffs_file'],
    conda:
        '../envs/Seurat_AB_seq.yaml'
    script:
        '../scripts/seur_filter_cells_n_genes.R'


rule seur_remove_doublets_Chord:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_after_filtering.rds"
    output:
        seur="results/Seurat_AB_seq/Seur_Tum_no_doublets.rds"
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        conds=config['Seurat_AB_seq']['conds']
    script:
        '../scripts/seur_remove_doublets_Chord.R'


rule seur_regress_out_cell_cycle_n_mito_prct:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_no_doublets.rds"
    output:
        seur=protected("results/Seurat_AB_seq/Seur_Tum_preprocessed_final.rds"),
        plot="results/Seurat_AB_seq/regressing_out_cell_cycle_n_mito_prct.pdf"
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        conds=config['Seurat_AB_seq']['conds'],
        cellcycle_genes=config['Seurat_AB_seq']['cellcycle_genes']
    script:
        '../scripts/seur_regress_out_cell_cycle_n_mito_prct.R'


rule seur_elbow_plots:
    input: seur="results/Seurat_AB_seq/Seur_Tum_preprocessed_final.rds"
    output: elbow_plots="results/Seurat_AB_seq/Tum_elbow_plots.pdf"
    conda: '../envs/Seurat_AB_seq.yaml'
    script: '../scripts/seur_elbow_plots.R'

rule seur_cluster_cells_Tum:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_preprocessed_final.rds",
        markers="config/markers.csv",
    output:
        seur=protected("results/Seurat_AB_seq/Seur_Tum_clustered.rds"),
        umap_all_resol="results/Seurat_AB_seq/Tum_UMAP_all_resol.pdf",
        umap_all_resol_split_by_sample="results/Seurat_AB_seq/Tum_UMAP_all_resol_split_by_sample.pdf",
        umap_markers="results/Seurat_AB_seq/Tum_UMAP_markers_split_by_sample.pdf",
        adt_heatmap="results/Seurat_AB_seq/Tum_heatmap_all_resol_adt_markers.pdf",
        rna_top_10_heatmap="results/Seurat_AB_seq/Tum_heatmap_all_resol_top_10_up_n_down_regulated_rna_markers.pdf",
        rna_top_20_heatmap="results/Seurat_AB_seq/Tum_heatmap_all_resol_top_20_up_n_down_regulated_rna_markers.pdf",
        adt_markers="results/Seurat_AB_seq/Tum_all_resol_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_all_resol_rna_markers.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        pca_dims=config['Seurat_AB_seq']['pca_dims'],
        resolutions=np.arange(0.1, 1.6, 0.2),
        recluster=False,
    script:
        '../scripts/seur_cluster_cells.R'



rule seur_gate_on_C5_at_resol_03:
    input: seur="results/Seurat_AB_seq/Seur_Tum_clustered.rds",
    output: seur="results/Seurat_AB_seq/Seur_Tum_clustered_gated_on_resol_0.3_C5.rds",
    conda: '../envs/Seurat_AB_seq.yaml'
    params: resolution=0.3, cluster_id=5,
    script: '../scripts/seur_gate_on_clusters_at_resol.R'



rule seur_recluster_cells_v1:
    input:
        seur="results/Seurat_AB_seq/Seur_Tum_clustered_gated_on_resol_0.3_C5.rds",
        markers="config/markers.csv"
    output:
        seur=protected("results/Seurat_AB_seq/Tum_CD8_v1/Seur_Tum_CD8_clustered.rds"),
        umap_all_resol="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_UMAP_all_resol.pdf",
        umap_all_resol_split_by_sample="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_UMAP_all_resol_split_by_sample.pdf",
        umap_markers="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_UMAP_markers_split_by_sample.pdf",
        adt_heatmap="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_heatmap_all_resol_adt_markers.pdf",
        rna_top_10_heatmap="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_heatmap_all_resol_top_10_up_n_down_regulated_rna_markers.pdf",
        rna_top_20_heatmap="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_heatmap_all_resol_top_20_up_n_down_regulated_rna_markers.pdf",
        adt_markers="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_all_resol_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_CD8_v1/Tum_CD8_all_resol_rna_markers.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        pca_dims=config['Seurat_AB_seq']['pca_dims'],
        resolutions=np.arange(0.1, 1.6, 0.2),
        recluster=True,
    script:
        '../scripts/seur_cluster_cells.R'


rule seur_gate_out_C4_in_resol_03:
    input: seur="results/Seurat_AB_seq/Tum_CD8_v1/Seur_Tum_CD8_clustered.rds",
    output: seur="results/Seurat_AB_seq/Tum_CD8_v1/Seur_Tum_CD8_clustered_without_C4.rds",
    conda: '../envs/Seurat_AB_seq.yaml'
    params: resolution=0.3, cluster_id=[0,1,2,3,5],
    script: '../scripts/seur_gate_on_clusters_at_resol.R'


rule seur_recluster_cells_v2:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v1/Seur_Tum_CD8_clustered_without_C4.rds",
        markers="config/markers.csv"
    output:
        seur=protected("results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered.rds"),
        umap_all_resol="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_UMAP_all_resol.pdf",
        umap_all_resol_split_by_sample="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_UMAP_all_resol_split_by_sample.pdf",
        umap_markers="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_UMAP_markers_split_by_sample.pdf",
        adt_heatmap="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_heatmap_all_resol_adt_markers.pdf",
        rna_top_10_heatmap="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_heatmap_all_resol_top_10_up_n_down_regulated_rna_markers.pdf",
        rna_top_20_heatmap="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_heatmap_all_resol_top_20_up_n_down_regulated_rna_markers.pdf",
        adt_markers="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_all_resol_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_all_resol_rna_markers.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        pca_dims=config['Seurat_AB_seq']['pca_dims'],
        resolutions=np.arange(0.1, 1.6, 0.2),
        recluster=True,
    script:
        '../scripts/seur_cluster_cells.R'


rule seur_gate_out_C5_in_resol_05:
    input: seur="results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered.rds",
    output: seur="results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered_without_C5.rds",
    conda: '../envs/Seurat_AB_seq.yaml'
    params: resolution=0.5, cluster_id=[0,1,2,3,4],
    script: '../scripts/seur_gate_on_clusters_at_resol.R'


rule seur_recluster_cells_v3:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered_without_C5.rds",
        markers="config/markers.csv"
    output:
        seur=protected("results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds"),
        umap_all_resol="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_UMAP_all_resol.pdf",
        umap_all_resol_split_by_sample="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_UMAP_all_resol_split_by_sample.pdf",
        umap_markers="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_UMAP_markers_split_by_sample.pdf",
        adt_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_heatmap_all_resol_adt_markers.pdf",
        rna_top_10_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_heatmap_all_resol_top_10_up_n_down_regulated_rna_markers.pdf",
        rna_top_20_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_heatmap_all_resol_top_20_up_n_down_regulated_rna_markers.pdf",
        adt_markers="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_all_resol_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_all_resol_rna_markers.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        pca_dims=config['Seurat_AB_seq']['pca_dims'],
        resolutions=np.arange(0.3, 0.9, 0.1),
        recluster=True,
    script:
        '../scripts/seur_cluster_cells.R'


rule seur_visulalize_umap_n_heatmap_for_one_resol:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds",
        all_resol_adt_markers="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_all_resol_adt_markers.csv",
        all_resol_rna_markers="results/Seurat_AB_seq/Tum_CD8_v3/Tum_CD8_all_resol_rna_markers.csv",
    output:
        umap_split_by_sample="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/UMAP_Tum_CD8_split_by_sample.pdf",
        adt_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/heatmap_Tum_CD8_adt_markers.pdf",
        rna_top_10_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/heatmap_Tum_CD8_top_10_up_regulated_rna_markers.pdf",
        rna_top_20_heatmap="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/heatmap_Tum_CD8_top_20_up_regulated_rna_markers.pdf",
        adt_markers="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/Tum_CD8_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_CD8_v3/res_{res}/Tum_CD8_rna_markers.csv",
        clust_prct_barplot="results//Seurat_AB_seq/Tum_CD8_v3/res_{res}/barplot_cluster_percentages.pdf",
        clust_prct_csv="results//Seurat_AB_seq/Tum_CD8_v3/res_{res}/barplot_cluster_percentages.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        resolution="{res}",
    script:
        '../scripts/seur_visualize_umap_n_heatmap_for_one_resol.R'



rule seur_identify_DEGs:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds",
    output:
        degs_without_cutoff="results/Seurat_AB_seq/Tum_CD8_v3/res_{resolution}/Tum_CD8_DEGs_without_cutoff.csv",
        degs_with_cutoff="results/Seurat_AB_seq/Tum_CD8_v3/res_{resolution}/Tum_CD8_DEGs_with_cutoff.csv",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        org=config['Seurat_AB_seq']['org'],
        conds=config['Seurat_AB_seq']['conds'],
        resolution="{resolution}",
    script:
        '../scripts/seur_identify_DEGs.R'


rule seur_identify_DEGs_of_KO_C3_vs_WT_C0:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds",
    output:
        degs_without_cutoff="results/Seurat_AB_seq/Tum_CD8_v3/res_{resol}/Tum_CD8_DEGs_C3.KO_vs_C0_WT_without_cutoff.csv",
        degs_with_cutoff="results/Seurat_AB_seq/Tum_CD8_v3/res_{resol}/Tum_CD8_DEGs_C3.KO_vs_C0_WT_with_cutoff.csv",
        htmp_DEGs_of_interest="results/Seurat_AB_seq/Tum_CD8_v3/res_{resol}/top_40_DEGs_KO.C3_vs_WT.C0_heatmap.pdf",
        htmp_DEGs="results/Seurat_AB_seq/Tum_CD8_v3/res_{resol}/DEGs_KO.C3_vs_WT.C0_heatmap.pdf",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        clust_1='Tum_KO.3',
        clust_2='Tum_WT.0',
        cluster_by='wsnn_res.0.4',
    script:
        '../scripts/seur_KO.C3_vs_WT.C0_analyses.R'


rule export_seur_data_to_SeqGeq:
    input:
        seur="results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered.rds",
        degs="results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/Tum_CD8_DEGs_with_cutoff.csv",
        adt_markers="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_all_resol_adt_markers.csv",
        rna_markers="results/Seurat_AB_seq/Tum_CD8_v2/Tum_CD8_all_resol_rna_markers.csv",
    output:
        seur_counts="results/SeqGeq/AB_seq_Tum_CD8_v2/res_{resolution}/CD8_counts.csv",
        seur_norm_counts="results/SeqGeq/AB_seq_Tum_CD8_v2/res_{resolution}/CD8_norm_counts.csv",
        seur_metadata="results/SeqGeq/AB_seq_Tum_CD8_v2/res_{resolution}/CD8_MetaData.csv",
    params:
        out_dir="results/SeqGeq/AB_seq_Tum_CD8_v2/res_{resolution}/",
        resolution="{resolution}",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    script:
        '../scripts/export_seur_data_to_SeqGeq.R'


rule seur_visualize_DEGs:
    input:
        degs_without_cutoff="results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/Tum_CD8_DEGs_without_cutoff.csv",
        seur="results/Seurat_AB_seq/Tum_CD8_v2/Seur_Tum_CD8_clustered.rds",
        genes_of_interest="resources/from_sabelo/CD8_T_cells_genes_of_interest.csv",
        DE_ISGs="resources/interferome_DB/ISGs_GeneSearchResults.txt"
    output:
        volcano_plot="results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/DEGs_volcanoplot.pdf",
        interactive_volcano_plot="results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/DEGs_interactive_volcanoplot.html",
        selected_genes_tiles = "results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/selected_genes_tiles.pdf",
        heatmap_DE_ISGs="results/Seurat_AB_seq/Tum_CD8_v2/res_{resolution}/DE_ISGs_Interferome_DB_heatmap.pdf",
    conda:
        '../envs/Seurat_AB_seq.yaml'
    params:
        conds=config['Seurat_AB_seq']['conds'],
        resolution="{resolution}",
    script:
        '../scripts/seur_visualize_DEGs.R'
