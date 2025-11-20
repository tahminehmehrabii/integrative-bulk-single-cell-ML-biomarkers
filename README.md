# Integrative Bulk Single Cell ML Biomarkers
### Single-cell LSCC pipeline 

This script implements a stepwise Seurat-based pipeline for multi-sample LSCC scRNA-seq (GSE150321).  

#### Inputs & basic setup
- **Inputs:** Raw per-sample count matrices in `.csv.gz` format under `data_dir`.
- **Arguments:**  
  - `data_dir` – folder containing the raw `.csv.gz` files (default: `D:/Single Cell/GSE150321_RAW`)  
  - `project_id` – project name used as prefix for outputs (default: `LSCC_GSE150321`)

---

### Stage-by-stage overview

1. **Stage 1 – Load raw counts (per sample)**  
   - Reads all `.csv.gz` files in `data_dir`.  
   - Builds one Seurat object per sample using raw counts.  
   - Stores `sample_id` in metadata.  
   - Exports per-sample cell counts before QC (`*_Stage1_raw_cell_counts.csv`).

2. **Stage 2 – Per-sample QC (nFeature, nCount, percent.mt)**  
   - Computes `percent.mt` per cell.  
   - Filters cells by fixed thresholds on `nFeature_RNA`, `nCount_RNA`, and mitochondrial %.  
   - Exports per-sample cell counts after QC (`*_Stage2_postQC_cell_counts.csv`).

3. **Stage 3 – Merge all samples**  
   - Merges all QC’d Seurat objects into one integrated object.  
   - Keeps `sample_id` in the merged metadata.  
   - Uses Seurat v5 `JoinLayers` and sets `RNA` as default assay.  
   - Saves merged metadata (`*_Stage3_meta_merged.csv`) and merged Seurat object (`*_Stage3_merged_seurat.rds`).

4. **Stage 4 – Normalization, HVG selection, scaling, PCA**  
   - `NormalizeData`, `FindVariableFeatures` (3,000 HVGs), `ScaleData`, `RunPCA`.  
   - Exports HVG gene list (`*_Stage4_HVG_genes.csv`).  
   - Generates PCA elbow plot (`*_Stage4_PCA_Elbow.png`).  
   - Saves post-PCA Seurat object (`*_Stage4_postPCA_seurat.rds`).

5. **Stage 5 – Graph-based clustering & UMAP**  
   - Builds kNN graph (`FindNeighbors`) and clusters (`FindClusters`, resolution 0.5).  
   - Runs UMAP on selected PCs.  
   - Initializes `celltype_main` as cluster ID.  
   - Exports UMAP colored by cluster and by sample (`*_Stage5_UMAP_clusters.png`, `*_Stage5_UMAP_samples.png`).  
   - Saves metadata with cluster labels (`*_Stage5_meta_with_clusters.csv`) and Seurat object (`*_Stage5_postUMAP_seurat.rds`).

6. **Stage 6 – Global cluster markers (up- and down-regulated)**  
   - Runs `FindAllMarkers` across clusters (no `only.pos` filter).  
   - Adds `direction` column (“up” / “down” based on `avg_log2FC`).  
   - Exports full marker table (`*_Stage6_cluster_markers_global_up_down.csv`).

7. **Stage 7 – Canonical lineage DotPlot**  
   - Defines marker panels for T cells, B/plasma, myeloid, tumor/epithelial, fibroblasts, endothelial, NK, mast, proliferation, TAM/SPP1+ macrophages.  
   - Filters markers to those present in the Seurat object.  
   - Exports the marker panel used (`*_Stage7_marker_panel_used.csv`).  
   - Generates a lineage DotPlot across clusters (`*_Stage7_lineage_DotPlot.png`).

8. **Stage 8 – Manual mapping of clusters to major cell types**  
   - Manually maps each cluster ID to a major lineage (e.g. “Tumor/Epithelial”, “T cells”, “Myeloid”, etc.).  
   - Updates `Idents` and `celltype_main` accordingly.  
   - Exports UMAP colored by main cell type (`*_Stage8_UMAP_celltypes.png`).  
   - Writes updated metadata with `celltype_main` (`*_Stage8_meta_with_celltypes.csv`).

9. **Stage 9 – Lineage-specific subsetting & reclustering**  
   - For each lineage (T cells, B/plasma cells, myeloid, tumor/epithelial, fibroblasts, endothelial):  
     - Subsets cells of that lineage.  
     - Re-runs normalization, HVG selection, scaling, PCA, neighbor graph, clustering, and UMAP.  
     - Generates UMAP by cluster and by sample (`*_Stage9_<Lineage>_UMAP_*.png`).  
     - Runs `FindAllMarkers` within each lineage and exports markers (`*_Stage9_<Lineage>_markers.csv`).  
     - Exports lineage metadata (`*_Stage9_<Lineage>_meta.csv`) and Seurat object (`*_Stage9_<Lineage>_seurat.rds`).

10. **Stage 10 – Save global annotated object**  
    - Saves the fully annotated global Seurat object with cluster and lineage labels (`*_Stage10_Global_annotated_seurat.rds`).

11. **Stage 11 – Epithelial labeling (Malignant vs Keratinocyte_like)**  
    - Works on the epithelial/tumor subset (`epith_obj`).  
    - Computes module scores per cluster for two predefined epithelial signatures:  
      - “Malignant” (proliferation + EMT/epithelial markers).  
      - “Keratinocyte_like” (squamous/keratinization markers).  
    - Automatically assigns each epithelial cluster to **Malignant** or **Keratinocyte_like** based on module-score differences.  
    - Exports:  
      - Per-cell cluster/group labels (`*_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_cells.csv`).  
      - Cluster-level module-score summary (`*_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_autoSelection.csv`).

12. **Stage 12 – Epithelial DEGs (Malignant vs Keratinocyte_like, strict)**  
    - Reads Stage 11 labels and applies them at the cell level within `epith_obj`.  
    - Ensures minimum cell counts per group, then runs `FindMarkers` (Wilcoxon) with strict thresholds (logFC, min.pct).  
    - Exports:  
      - Group sizes (`*_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs_group_counts.csv`).  
      - Full DEG table including `direction` (“up” / “down”) (`*_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs.csv`).

---
### Bulk LSCC pipeline (GSE127165 + external validation cohorts)

This script implements a multi-stage bulk RNA-seq pipeline for LSCC, using one training cohort (GSE127165) and two external validation cohorts (GSE142083, GSE130605).  

#### Inputs & basic setup
- **Base path:**  
  - `project_path = "D:/LSCC"`
- **Required inputs (under `project_path`):**
  - Raw bulk count matrices (training + test cohorts)  
    - `GSE127165_raw_counts_GRCh38.p13_NCBI.csv`  
    - `GSE142083_raw_counts_GRCh38.p13_NCBI.csv`  
    - `GSE130605_raw_counts_GRCh38.p13_NCBI.csv`
  - Gene annotation tables (with `GeneID`, `Symbol`, `GeneType`)  
    - `GSE127165_annot.csv`  
    - `GSE142083_annot.csv`  
    - `GSE130605_annot.csv`
  - Phenotype/metadata (per sample; must include `Sample`, `group`)  
    - `Pheno_Data_GSE127165.csv`  
    - `Pheno_Data_GSE142083.csv`  
    - `Pheno_Data_GSE130605.csv`
  - Single-cell DEGs for epithelial malignant vs keratinocyte-like comparison  
    - `D:/Single Cell/GSE150321_RAW/LSCC_GSE150321_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs.csv`

---

### Stage-by-stage overview (bulk pipeline)

#### Stage 1 – Preprocessing (TMM, annotation, transpose, metadata, filtering)

**Stage 1 (part 1): TMM normalization → log2CPM**  
- For each cohort (train + two tests):  
  - Reads raw count matrix (`GeneID × sample`).  
  - Creates `edgeR::DGEList`, applies TMM normalization.  
  - Converts to log2 CPM (`cpm(log=TRUE, prior.count=1)`).  
  - Writes log2CPM matrices:  
    - `GSE127165_Stage1b_log2CPM_TMM.csv`  
    - `GSE142083_Stage1b_log2CPM_TMM.csv`  
    - `GSE130605_Stage1b_log2CPM_TMM.csv`

**Stage 1 (part 2): Add gene symbols and aggregate duplicates**  
- Merges log2CPM matrices with annotation (`GeneID → Symbol`).  
- Keeps only rows with non-empty `Symbol`.  
- Aggregates duplicate symbols by mean expression.  
- Outputs symbol-level matrices:  
  - `GSE127165_Stage2_log2CPM_withSymbol_agg.csv`  
  - `GSE142083_Stage2_log2CPM_withSymbol_agg.csv`  
  - `GSE130605_Stage2_log2CPM_withSymbol_agg.csv`

**Stage 1 (part 3): Transpose (samples as rows, genes as columns)**  
- Transposes each symbol-level matrix to `Sample × gene`.  
- Adds explicit `Sample` column.  
- Writes transposed files:  
  - `GSE127165_Stage3_transposed.csv`  
  - `GSE142083_Stage3_transposed.csv`  
  - `GSE130605_Stage3_transposed.csv`

**Stage 1 (part 4): Attach phenotype / metadata**  
- Merges expression matrices with phenotype tables by `Sample`.  
- Reorders columns to put `Sample` and `group` first.  
- Outputs metadata-augmented matrices:  
  - `GSE127165_Stage4_withMetadata.csv`  
  - `GSE142083_Stage4_withMetadata.csv`  
  - `GSE130605_Stage4_withMetadata.csv`

**Stage 1 (part 5): mRNA-only and removal of RP / MT / KRT genes**  
- From annotation, selects protein-coding / mRNA genes.  
- Removes ribosomal (`RPL*`, `RPS*`), mitochondrial (`MT-*`), and keratin (`KRT*`) genes.  
- Filters each dataset to the intersection of expression features and allowed mRNA symbols.  
- Outputs final filtered expression matrices:  
  - `GSE127165_Stage5_filtered_matrix.csv`  
  - `GSE142083_Stage5_filtered_matrix.csv`  
  - `GSE130605_Stage5_filtered_matrix.csv`

---

#### Stage 2 – Train/Validation split (GSE127165)

- Uses the filtered training matrix (`GSE127165_Stage5_filtered_matrix.csv`).  
- Ensures `Sample` and `group` columns are present and non-missing.  
- Uses `caret::createDataPartition` (70/30 stratified on `group`) to split into:  
  - **Train:** `GSE127165_Stage6_Split_Train.csv`  
  - **Valid:** `GSE127165_Stage6_Split_Valid.csv`

---

#### Stage 3 – Differential expression, QC plots, and enrichment (GO/KEGG)

**Stage 3 (part 1): DEG with limma + voomWithQualityWeights**  
- Uses:  
  - Train split (`GSE127165_Stage6_Split_Train.csv`)  
  - Raw counts (`GSE127165_raw_counts_...csv`)  
  - Annotation (`GSE127165_annot.csv`)  
- Steps:  
  - Maps `GeneID → Symbol` and aggregates counts at symbol level.  
  - Restricts to genes present in Stage 6 train features.  
  - Builds `DGEList`, applies TMM normalization.  
  - Runs `voomWithQualityWeights` with design matrix for `Normal` vs `Tumor`.  
  - Applies presence filter on logCPM (`min_lcpm`, `min_prop`).  
  - Fits limma linear model, sets contrast (Tumor – Normal), applies empirical Bayes.  
  - Writes:  
    - Full DEG table: `GSE127165_Stage7_DE_topTable_full.csv`  
    - Strict filtered DEGs (FDR/logFC/AveExpr thresholds): `GSE127165_Stage7_DE_topTable_strict.csv`  
    - Up / down / combined gene lists:  
      - `GSE127165_Stage7_DE_up_strict.txt`  
      - `GSE127165_Stage7_DE_down_strict.txt`  
      - `GSE127165_Stage7_DE_updown_strict.txt`  
    - Sample quality weights: `GSE127165_Stage7_sample_quality_weights.csv`

**Stage 3 (part 2): Volcano plots**  
- Reads full DEG table from Stage 7.  
- Computes –log10(FDR) and classifies genes as Up/Down/NS.  
- Writes volcano data table: `GSE127165_Stage8_volcano_data_full.csv`.  
- Generates:  
  - Full volcano: `GSE127165_Stage8_volcano_full.png`  
  - Optional strict volcano (if strict table exists): `GSE127165_Stage8_volcano_strict.png`

**Stage 3 (part 3): Sample–sample correlation heatmaps**  
- Uses train split (`GSE127165_Stage6_Split_Train.csv`).  
- Converts expression to numeric matrix (`Sample × gene`).  
- Computes correlation between samples:  
  - Pearson (standard heatmap).  
  - Spearman with contrast enhancement (clipping to quantiles).  
- Generates sample correlation heatmap PNG (overwritten with enhanced version):  
  - `GSE127165_Stage9A_train_sample_correlation.png`

**Stage 3 (part 4): GO & KEGG enrichment (top-10 rich-factor bubble plots)**  
- Uses strict bulk DEGs (symbols) from Stage 7 or fallback list.  
- Performs enrichment:  
  - GO (BP/CC/MP; `enrichGO`, `org.Hs.eg.db`).  
  - KEGG (human; `enrichKEGG`, ENTREZ IDs).  
- Builds bubble-plot data (rich factor = Count/set size).  
- Outputs:  
  - GO results: `Stage10_GO_results.csv`  
  - GO bubble plot: `Stage10_GO_Top10_RichFactor.png`  
  - KEGG results: `Stage10_KEGG_results.csv`  
  - KEGG bubble plot: `Stage10_KEGG_Top10_RichFactor.png`

---

#### Stage 4 – Co-expression network analysis (CEMiTool)

- Uses training split (`GSE127165_Stage6_Split_Train.csv`).  
- Builds gene × sample expression matrix and phenotype (`Normal` vs `Tumor`).  
- Runs `CEMiTool` with example GMT and interaction files.  
- Outputs directly in `project_path`:  
  - CEMiTool HTML report (`generate_report`).  
  - Module and enrichment tables (`write_files`).  
  - All CEMiTool plots (`save_plots`).  
  - Key PNGs: `beta_r2.png`, `mean_k.png`, `gsea.png`, `ora_M1.png`.  
- Extracts genes from a selected co-expression module (e.g. `M1`):  
  - Reads `module.tsv` and writes module genes: `modulesGenes.txt`.

---

#### Stage 5 – Integration: Bulk DEGs × co-expression module × scRNA DEGs

- Uses:  
  - Bulk strict DEGs (Stage 7): `GSE127165_Stage7_DE_updown_strict.txt`  
  - Module genes (from CEMiTool): `modulesGenes.txt`  
  - scRNA epithelial DEGs (Stage 12, single-cell pipeline):  
    - `LSCC_GSE150321_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs.csv`
- Extracts gene-symbol sets:  
  - `Bulk_DEGs`  
  - `Module_M1`  
  - `scRNA_DEGs`
- Generates 3-set Venn diagram:  
  - `GSE127165_Stage12_Venn_All3.png`  
- Computes triple overlap across all three gene sets and writes:  
  - `GSE127165_Stage12_Overlap_All3.csv` (key integrative biomarker candidates).

---

#### Stage 6 – LASSO-based biomarker selection + single-gene ROC

**Stage 6 (part 1): LASSO logistic regression on triple-overlap genes**  
- Uses:  
  - Train split (`GSE127165_Stage6_Split_Train.csv`).  
  - Triple-overlap genes (`GSE127165_Stage12_Overlap_All3.csv`).  
- Steps:  
  - Filters expression matrix to overlapping genes.  
  - Encodes `group` as Normal (0) vs Tumor (1).  
  - Fits logistic LASSO (`glmnet`, `family="binomial"`) with 10-fold CV.  
  - Selects `lambda.min`.  
- Outputs:  
  - Lambda–deviance curve: `GSE127165_Stage13_Lambda_BinomialDeviance.png`  
  - Selected λ: `GSE127165_Stage13_lambda_min.txt`  
  - Non-zero coefficients and gene list:  
    - `GSE127165_Stage13_coefficients.csv`  
    - `GSE127165_Stage13_lassoGenes.txt`  
  - Horizontal barplot of non-zero coefficients:  
    - `GSE127165_Stage13_nonZeroCoefGenes.png`

**Stage 6 (part 2): ROC curves for selected genes (Train / Valid)**  
- Uses LASSO-selected genes from Stage 13.  
- Keeps genes present in both Train and Valid matrices.  
- For each gene, computes ROC in Train and Valid using `pROC::roc`.  
- Generates:  
  - Multi-panel ROC for all selected genes: `GSE127165_Stage14_ROC_all.png`  
  - ROC subset (removing pre-specified indices): `GSE127165_Stage14_ROC_subset.png`

---

#### Stage 7 – Machine learning models (RF / SVM / GBM) + external validation

- Uses:  
  - Train & Valid splits (Stage 6).  
  - Test cohorts (GSE142083, GSE130605) filtered at Stage 5.  
  - Biomarker genes from Stage 13 (optionally dropping some indices).  
- Harmonizes feature space across Train, Valid, Test1, Test2, and biomarker list (`common_genes`).  
- Normalizes group labels to `Normal` / `Tumor`.  

**Models trained:**
1. **Random Forest** (`randomForest`)  
   - Tunes `mtry` with `tuneRF`.  
   - Trains RF with 500 trees and variable importance.

2. **Linear SVM** (`e1071::svm`)  
   - Tunes `cost` and `gamma` on Train.  
   - Trains SVM with probability estimates.

3. **Gradient Boosting (GBM via H2O)**  
   - Initializes H2O, imports temporary CSVs for Train/Valid/Test.  
   - Trains `h2o.gbm` on Train.  

**Evaluation (per model):**
- For each dataset (Train, Valid, GSE142083, GSE130605):  
  - Predicts class labels and Tumor probability.  
  - Saves confusion matrices:  
    - e.g. `RF.train.confusion.csv`, `SVM.valid.confusion.csv`, `GBM.test_GSE130605.confusion.csv`  
  - Computes metrics (AUC, CI, accuracy, sensitivity, specificity, precision, recall, F1):  
    - e.g. `RF.valid.metrics.csv`, `SVM.test_GSE142083.metrics.csv`, `GBM.test_GSE130605.metrics.csv`
  - Builds ROC curves using `pROC::roc`.

**Combined ROC plots (model comparison):**
- Plots RF vs SVM vs GBM AUC on each dataset:  
  - `ROC_Training.png`  
  - `ROC_Validation.png`  
  - `ROC_Test_GSE142083.png`  
  - `ROC_Test_GSE130605.png`

**Model export & cleanup:**
- Saves trained models to disk:  
  - `RF_model.rds`  
  - `SVM_model.rds`  
  - H2O GBM model directory (`h2o.saveModel` in `project_path`)  
- Shuts down H2O (`h2o.shutdown(prompt = FALSE)`).

