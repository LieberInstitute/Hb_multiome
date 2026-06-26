#Saving the color palette used for consistency across figures


# my_colors_mid = c(Excit.Thal = "#4d55b7",
#   LHb.4 = "#00607A",
#   Inhib.Thal = "#9a9fe7",
#   Astrocyte = "#532222",
#   MHb.1.2 = "#92007C",
#   LHb.1 = "#008092",
#   OPC = "#829454",
#   Oligo =  "#384a08",
#   Microglia = "#141b02",
#   LHb.2.7 = "#6C9FA9",
#   Endo = "#d95f02",
#   LHb.1.3.4 = "#A8B8BC",
#   MHb.1 = "#A86A9A",
#   MHb.2 = "#BCA6B6",
#   MHb.3 = "#56204eff",
#   LHb.1.3 = "#C6C6C6",
#   Inhib_LHb_4.1 = "#8B0000",
#   Inhib_LHb_4.2 = "#DC143C",
#   Ependymal = "#f5a105ff"
# ) 
new_colors = MetBrewer::met.brewer(name="Moreau", n=13, type="continuous")

my_colors_mid = c(
  Astrocyte = "#972f2f",
  OPC = "#829454",
  Oligo =  "#384a08",
  Microglia = "#141b02",
  Endo = "#f65a45",
  Ependymal = "#dbb369",
  
  Excit.Thal = "#2e6296",
  Inhib.Thal = "#8DADCA",
  LHb.4 = "#082844",
  Inhib_LHb_4.1 = "#a93905",
  Inhib_LHb_4.2 = "#aa6f16",
  LHb.1.3 = "#527BAA",
  LHb.1 = "#0C383E",
  LHb.2.7 = "#ee9630",
  LHb.1.3.4 = "#306171",
  MHb.1 = "#5e0c01",
  MHb.1.2 = "#f67104",
  MHb.2 = "#943f02",
  MHb.3 = "#f4d5ab"
) 

my_colors_class <- c(
    LHb = "#0269ae",
    MHb = "#ad1d8c",
    `Non-neurons` = "#532222", 
    Thalamus = "#5906df"
    
)
