library("rsconnect")
library("here")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appDir = here("code", "07_iSEE_app"),
    appFiles = c("app.R", "sce_Habenula_iSEE.rds"),
    appName = "Habenula_multiome",
    account = "libd",
    server = "shinyapps.io"
)
