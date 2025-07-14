# Full Annotation Meta-data on Multimodal WNN Clusters

This repository contains the detailed annotation metadata for multimodal Weighted Nearest Neighbor (WNN) clustering, integrating RNA and ATAC modalities.  
It provides transparent documentation of the full decision process used to assign final cell-type labels, along with multiple supporting metrics.

---

## File Descriptions

This CSV file contains the **detailed annotation of each WNN cluster**.  
Prefix numbers (e.g. `(1)`, `(2)`, ...) indicate the sequential strategy steps used to determine or support the final cell-type annotations.

| Column Name                  | Description |
|-------------------------------|-------------|
| **cluster**                   | WNN multiome clusters (Leiden communities) computed with `k=30` neighbors at resolution `2`. |
| **(1) ct_frequency_in2ct**     | Most frequently matched cell-type class per cluster, based on the top 50 DE genes (FDR < 5%) mapped to the reference. |
| **frequency_repeated_ct**      | All repeated cell-type matches per cluster, determined by the same top 50 DE genes. |
| **(2) ct_CRegistration_fine**  | Manual label of significant cell-types observed from fine-resolution Spatial-Registration comparing multiome vs snRNA-seq data. |
| **(3) ct_MeanRatio_support**   | Most supported cell-type based on nearest mean-ratio scores across the top 10 genes. |
| **mean_ratio_repeated_ct**     | Highest individual cell-type match based on nearest mean-ratio scores. For example, `C.34.Oligo x5` means 5 genes match to `C.34.Oligo`. |
| **mean_ratio_top2**            | The two top cell-type matches by nearest mean-ratio scores. |
| **mean_ratio_detail**          | Full breakdown of all matches from nearest mean-ratio scores among the top 10 genes. E.g. `C.34.Oligo x5`. |
| **mean_ratio_range**           | Gene mean-ratio range specificity score, calculated as:<br>`MeanRatio = mean expr in highest non-target / mean expr in target celltype`.<br>See [DeconvoBuddies Marker Finding](https://www.bioconductor.org/packages/devel/bioc/vignettes/DeconvoBuddies/inst/doc/DeconvoBuddies.html#32_Marker_Finding). |
| **(4) HClust_support**         | Pairwise support from hierarchical clustering, showing bottom leaves from clustering WNN centroids in PCA space using the `ward.D2` method. |
| **ct_final**                   | **Final cell-type assignments** after integrating:<br>• cell-type frequency<br>• clustering registration vs snRNA-seq<br>• mean-ratio scores<br>• hierarchical clustering clades<br>• 1vsALL enrichment t-stats<br>• and spatial localization (Visium-HD) for ambiguous clusters. |
| **(5) description_support**    | Technical and biological rationale summarizing why the final annotation was chosen, especially in ambiguous cases, often referencing Visium-HD support. |
| **(6) number_cells**           | Number of cells in this cluster. |
| **cluster_percentage**         | Percentage of all cells assigned to this cluster. |


---

## Purpose

This annotation table serves to:
- Provide a **transparent record** of how final cell-type labels were assigned.
- Include intermediary data (mean-ratio scores, frequency matches, HClust support) to enable critical review and validation.
- Facilitate downstream analysis by supplying ready-to-use cluster annotations.

---

## How to cite or use

If you use this annotation metadata in your analysis or publication, please cite this repository (or associated manuscript) and acknowledge the multi-step approach combining:
- differential expression statistics,
- multi-modal clustering,
- spatial registration,
- hierarchical clustering on PCA space,
- and marker-based mean-ratio specificity scoring.

---

## Questions?

If you have questions or want to request additional files (e.g., raw DE results, marker mean-ratio tables), please open an issue or contact us directly.

---

**Thank you for using this resource!**

