library(tidyverse)
library(here)
library(sessioninfo)

file_map_path = here(
    'processed-data', '19_data_uploads', '01_file_map', 'map.csv'
)
demo_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/14_supp_tables/donor_demographics.csv'
TODO_var = 'whatever'

file_map = read_csv(file_map_path, show_col_types = FALSE)

demo_df = read_csv(demo_path, show_col_types = FALSE) |>
    dplyr::rename(
        subject_name = donor, age_value = age, postmortem_interval = PMI
    ) |>
    mutate(
        race = case_when(
            race == 'European' ~ 'White',
            race == 'African' ~ 'Black_or_African_American',
            TRUE ~ NA_character_
        )
    )

subject_df = file_map |>
    distinct(donor) |>
    dplyr::rename(subject_name = donor) |>
    mutate(
        subject_source = TODO_var,
        subject_source_id = subject_name,
        subject_source_catalog_number = NA,
        subject_type = 'individual',
        species = 'NCBI:txid9606',
        strain_subspecies_name = NA,
        subject_genotype = NA,
        grant_number = 'R01DA055823',
        grant_name = TODO_var,
        project_short_name = TODO_var,
        cohort_id = TODO_var,
        metadata_access = 'open',
        controlled_access_fields = NA,
        additional_assays = NA,
        age_unit = 'years',
        subject_event_name = 'postmortem',
        subject_comments = NA,
        time_of_death = 'not_known',
        autopsy_year = NA,
        transgender = 'no',
        ethnicity = TODO_var,
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
        tox_history_bzo = TODO_var,
        tox_history_bup = 'negative',
        tox_history_thc = 'negative',
        tox_history_coc = 'negative',
        tox_history_mtd = 'negative',
        tox_history_met = 'negative',
        tox_history_opi = TODO_var,
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
    