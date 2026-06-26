#Quantifying the RNAscope data of the VGAT probes in the habenula



library(here)




plot_path = here('plots','14_vgat_rnascope', '01_halo_data')
if (!dir.exists(plot_path)) dir.create(plot_path)


halo_path = here('processed-data','14_vgat_rnascope','HALO_Trimmed')


all_files = list.files(halo_path)
csv_files = all_files[grepl('.csv',all_files)]



current_data = data.table::fread(file.path(halo_path,csv_files[1]))

hist(current_data$`VGAT+`)
colnames(current_data) 


plot(current_data$`POU4F1_690 (Opal 690) Cell Intensity`, current_data$`SLC32A1_570 (Opal 570) Cell Intensity`)

plot(current_data$`SLC17A6_620 (Opal 620) Cell Intensity`, current_data$`SLC32A1_570 (Opal 570) Cell Intensity`)

plot(current_data$`SLC17A6_620 (Opal 620) Cell Intensity`, current_data$`POU4F1_690 (Opal 690) Cell Intensity`)


