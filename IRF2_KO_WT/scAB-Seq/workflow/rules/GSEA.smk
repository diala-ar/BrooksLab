rule GSEA_KO_C3_vs_WT_C0_of_custom_gs:
    input:
        tox_degs="resources/gene_sets/41586_2019_1325_MOESM3_ESM.csv",
        degs="results/Seurat_AB_seq/Tum_CD8_v3/res_0.4/Tum_CD8_DEGs_C3.KO_vs_C0_WT_without_cutoff.csv",
        seur="results/Seurat_AB_seq/Tum_CD8_v3/Seur_Tum_CD8_clustered.rds",
    output:
        gsea_plots="results/Seurat_AB_seq/Tum_CD8_v3/res_{resol}/GSEA/Tum_CD8_KO.C3_vs_WT.C0_custom_genesets_GSEA_plots.pdf",
    params:
        cluster_by='wsnn_res.0.4',
        clust_1='Tum_KO.3',
        clust_2='Tum_WT.0',
    conda: "../envs/Seurat_AB_seq.yaml",
    script: '../scripts/GSEA_of_custom_gs.R'
