library(tidyverse)
library(here)
library(openxlsx)
library(sessioninfo)

access_type = c('open', 'restricted')[
    as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
]

raw_map_path = here(
    'processed-data', '19_data_uploads', '01_file_map', 'map.csv'
)
processed_map_path = here(
    'processed-data', '19_data_uploads', '02_processed_data', 'map.csv'
)
out_path = here(
    'processed-data', '19_data_uploads', '03_manifest',
    sprintf('manifest_%s.xlsx', access_type)
)
demo_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/14_supp_tables/donor_demographics.csv'
data_subtype_opts = c(
    'aligned_fiducials', 'aligned_tissue_image', 'detected_tissue_image',
    'scalefactors_json', 'tissue_hires_image', 'tissue_lowres_image',
    'per_barcode_metrics', 'spatial_enrichment'
)
TODO_var = '[needs_a_value]'

dir.create(dirname(out_path), showWarnings = FALSE)

raw_map = read_csv(raw_map_path, show_col_types = FALSE) |>
    mutate(library_aliquot_name = paste0(library_id, '_1'))
processed_map = read_csv(processed_map_path, show_col_types = FALSE) |>
    mutate(library_aliquot_name = NA)

file_map = raw_map |>
    select(
        file_path, technique, library_aliquot_name, open_access, md5_checksum
    ) |>
    bind_rows(processed_map) |>
    filter(open_access == (access_type == 'open'))

demo_df = read_csv(demo_path, show_col_types = FALSE) |>
    dplyr::rename(
        subject_name = donor, age_value = age, postmortem_interval = PMI
    ) |>
    mutate(
        race = case_when(
            ancestry == 'European' ~ 'White',
            ancestry == 'African' ~ 'Black_or_African_American',
            TRUE ~ NA_character_
        )
    )

subject_df = raw_map |>
    distinct(donor) |>
    dplyr::rename(subject_name = donor) |>
    mutate(
        subject_source = TODO_var, # Need to be given a value from NeMO people
        subject_source_id = subject_name,
        subject_source_catalog_number = NA,
        subject_type = 'individual',
        species = 'NCBI:txid9606',
        strain_subspecies_name = NA,
        subject_genotype = NA,
        grant_number = 'R01DA055823',
        grant_name = 'R01DA055823_maynard',  # Need to be given a value from NeMO people
        project_short_name = TODO_var,  # Need to be given a value from NeMO people
        cohort_id = TODO_var,  # Need to be given a value from NeMO people
        metadata_access = 'open',
        controlled_access_fields = NA,
        additional_assays = NA,
        age_unit = 'years',
        subject_event_name = 'postmortem',
        subject_comments = NA,
        time_of_death = 'not_known',
        autopsy_year = NA,
        transgender = 'no',
        ethnicity = 'not_Hispanic_or_Latino',
        cause_of_death = NA,
        viral_status = 'presumed_none',
        hiv_strain = NA,
        years_with_hiv = 'not_applicable',
        cd4_counts = 'not_applicable',
        cd4_nadir = 'not_known',
        cognitive_status = NA,
        brain_pathology = NA,
        spinal_cord_pathology = NA,
        antiretroviral_therapy = 'no_art',
        plasma_viral_load_measurement = 'not_applicable',
        months_since_last_plasma_viral_measurement = 'not_known',
        tox_history_amp = 'negative',
        tox_history_bar = 'negative',
        tox_history_bzo = ifelse(subject_name == 'Br9017', 'positive', 'negative'),
        tox_history_bup = 'negative',
        tox_history_thc = 'negative',
        tox_history_coc = 'negative',
        tox_history_mtd = 'negative',
        tox_history_met = 'negative',
        tox_history_opi = ifelse(subject_name == 'Br9902', 'not_known', 'negative'),
        tox_history_oxy = 'negative',
        tox_history_pcp = 'negative',
        tox_history_tca = 'negative',
        tox_history_comments = NA,
        postmortem_toxicology_nms = NA,
        postmortem_toxicology_nsu = NA,
        diagnosis_within_two_years_of_death = 'Medical record review: No Substance Use Disorders',
        lifetime_diagnosis = 'Medical record review: No Substance Use Disorders',
        overdose_death = 'no'
    ) |>
    left_join(demo_df, by = 'subject_name') |>
    mutate(sex = ifelse(sex == 'M', 'male', 'female')) |>
    select(
        subject_name,
        subject_source,
        subject_source_id,
        subject_source_catalog_number,
        subject_type,
        species,
        strain_subspecies_name,
        subject_genotype,
        grant_number,
        grant_name,
        project_short_name,
        cohort_id,
        metadata_access,
        controlled_access_fields,
        sex,
        additional_assays,
        age_value,
        age_unit,
        subject_event_name,
        subject_comments,
        postmortem_interval,
        time_of_death,
        autopsy_year,
        transgender,
        race,
        ethnicity,
        cause_of_death,
        viral_status,
        hiv_strain,
        years_with_hiv,
        cd4_counts,
        cd4_nadir,
        cognitive_status,
        brain_pathology,
        spinal_cord_pathology,
        antiretroviral_therapy,
        plasma_viral_load_measurement,
        months_since_last_plasma_viral_measurement,
        tox_history_amp,
        tox_history_bar,
        tox_history_bzo,
        tox_history_bup,
        tox_history_thc,
        tox_history_coc,
        tox_history_mtd,
        tox_history_met,
        tox_history_opi,
        tox_history_oxy,
        tox_history_pcp,
        tox_history_tca,
        tox_history_comments,
        postmortem_toxicology_nms,
        postmortem_toxicology_nsu,
        diagnosis_within_two_years_of_death,
        lifetime_diagnosis,
        overdose_death
    )

sample_df = raw_map |>
    distinct(sample_id, donor) |>
    dplyr::rename(sample_name = sample_id, subject_name = donor) |>
    mutate(
        sample_source = TODO_var, # Need to be given a value from NeMO people
        sample_source_id = sample_name,
        subject_event_name = 'postmortem',
        project_short_name = TODO_var, # Need to be given a value from NeMO people
        lab = 'Maynard',
        sample_type = 'individual',
        parent_sample_name = NA,
        anatomical_site = 'UBERON:0001904', # habenula
        sample_storage_length = NA,
        sample_comments = NA
    ) |>
    select(
        sample_name,
        sample_source,
        sample_source_id,
        subject_name,
        subject_event_name,
        project_short_name,
        lab,
        sample_type,
        parent_sample_name,
        anatomical_site,
        sample_storage_length,
        sample_comments
    )

library_df = raw_map |>
    distinct(library_id, sample_id, technique) |>
    dplyr::rename(library_name = library_id, parent_name = sample_id) |>
    mutate(
        library_aliquot_name = paste0(library_name, '_1'), # NeMO people will answer a question about this
        library_type = 'individual',
        subspecimen_type = case_when(
            str_detect(technique, 'Multiome') ~ 'nuclei',
            str_detect(technique, 'Visium') ~ 'bulk'
        ),
        project_short_name = TODO_var, # Need to be given a value from NeMO people
        lab = 'Maynard',
        parent_type = 'sample',
        library_demultiplexing = NA,
        library_batch = NA,
        library_comments = NA
    ) |>
    select(
        library_name,
        library_aliquot_name,
        library_type,
        technique,
        subspecimen_type,
        project_short_name,
        lab,
        parent_name,
        parent_type,
        library_demultiplexing,
        library_batch,
        library_comments
    )

file_df = file_map |>
    mutate(
        program = 'SCORCH',
        file_name = basename(file_path),
        summary_file = 'no',
        grant_number = 'R01DA055823',
        grant_name = 'R01DA055823_maynard', # Need to be given a value from NeMO people
        project_short_name = TODO_var, # Need to be given a value from NeMO people
        lab = 'Maynard',
        data_type = case_when(
            str_detect(file_name, '\\.fastq') ~ 'sequence_reads',
            str_detect(file_name, '\\.(tiff|jpg|png)$') ~ 'image',
            str_detect(file_name, '\\.(rds|h5|mtx)(\\.gz)?$') ~ 'counts',
            TRUE ~ 'analysis_metrics'
        ),
        file_derived_from = NA,
        species = 'NCBI:txid9606',
        file_format = str_extract(
                file_name, '\\.([A-Za-z0-9]+)(\\.gz)?$', group = 1
            ) |> str_replace('h5', 'H5'),
        data_subtype = case_when(
            str_detect(file_name, '_R1_[0-9]+\\.fastq') ~ 'r1_fastq',
            str_detect(file_name, '_R2_[0-9]+\\.fastq') ~ 'r2_fastq',
            str_detect(file_name, '_R3_[0-9]+\\.fastq') ~ 'r3_fastq',
            str_detect(file_name, '_I1_[0-9]+\\.fastq') ~ 'index1_fastq',
            str_detect(file_name, '_I2_[0-9]+\\.fastq') ~ 'index2_fastq',
            str_detect(file_name, 'tissue_positions') ~ 'tissue_positions',
            str_detect(file_name, '\\.(tsv|h5|rds)(\\.gz)?$') ~ str_extract(
                file_name, '\\.(tsv|h5|rds)(\\.gz)?$', group = 1
            ),
            str_detect(file_name, '\\.mtx(\\.gz)?$') ~ 'matrix_mex',
            str_detect(file_name, '\\.tif') ~ 'tif',
            str_detect(file_name, paste(data_subtype_opts, collapse = '|')) ~ str_extract(
                file_name, paste(data_subtype_opts, collapse = '|')
            )
        ),
        access = ifelse(open_access, 'open', 'restricted'),
        data_use_condition = 'DUO:0000004',
        data_use_specific_limit = NA,
        cohort_id = TODO_var, # Need to be given a value from NeMO people
        pipeline_name = NA,
        pipeline_rrid = NA,
        pipeline_version = NA,
        pipeline_container_url = NA,
        data_type_specific_tool = NA,
        genome_build = NA,
        gene_set_release = NA,
        sequencing_batch = NA,
        file_comments = NA
    ) |>
    select(
        program,
        file_name,
        summary_file,
        library_aliquot_name,
        technique,
        grant_number,
        grant_name,
        project_short_name,
        lab,
        data_type,
        file_derived_from,
        species,
        file_format,
        data_subtype,
        access,
        data_use_condition,
        data_use_specific_limit,
        cohort_id,
        md5_checksum,
        pipeline_name,
        pipeline_rrid,
        pipeline_version,
        pipeline_container_url,
        data_type_specific_tool,
        genome_build,
        gene_set_release,
        sequencing_batch,
        file_comments
    )

write.xlsx(
    list(
        subject = subject_df,
        sample = sample_df,
        library = library_df,
        file = file_df
    ),
    file = out_path
)

session_info()
