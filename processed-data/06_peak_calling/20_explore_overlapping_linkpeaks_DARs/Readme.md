# 📊 Processed Data Tables

This folder contains processed summary tables related to LinkPeak overlaps with differentially accessible regions (DARs) at an FDR threshold of **2.0**.

| File | Description |
|----|----|
| `overlaps_summary_linkPeak_DARs_unique_ct_FDR2.0.csv` | Summary of all LinkPeaks overlapping DARs at FDR 2.0. Each row represents a LinkPeak annotated with: total number of DARs (`n_DARs`), number of distinct cell types (`n_cell_types`), and their combined list (`cell_types`). Used to classify LinkPeaks as **Unique** (one cell type) or **Shared** (multiple). |
| `overlaps_summary_linkPeak_DARs_shared_2ct_FDR2.0.csv` | Subset of the above, containing only LinkPeaks shared by exactly **two** cell types. Used to highlight cell-type pairs that co-regulate the same accessible regions. |
| `overlaps_summary_linkPeak_ALL_shared_top3_freq_ct_2.0.csv` | Summary of the top 3 most frequent cell-type combinations for each shared level (`n_cell_types = 2–10`). Useful for identifying recurring regulatory partnerships among cell types. |
| `overlaps_linkPeaks_unique_shared_detailed_FDR2.0.csv` | Long-format table listing each `peak_id × cell_type` pair with its classification (**Unique** or **Shared**). Each LinkPeak appears once per associated cell type. Used for cell-type–specific barplots. |
| `overlaps_unique_linkPeak_DARs_detail_FDR2.0.csv` | Detailed table containing only **unique** LinkPeak–DAR overlaps. Each row corresponds to a single DAR linked to a unique LinkPeak in one cell type. Used for cell-type–specific “unique overlap” plots. |
| `overlaps_summary_linkPeak_DARs_shared_2ct_FDR2.0.csv` | Same structure as the unique file, but restricted to **shared** LinkPeaks (`n = 2`) to explore chromatin regions accessible in two cell types. |
| `overlaps_summary_linkPeak_DARs_unique_ct_FDR2.0.csv` | Summary counts of all **unique and shared** overlaps aggregated by cell type, used in stacked barplots comparing Unique vs. Shared overlap frequencies. |
