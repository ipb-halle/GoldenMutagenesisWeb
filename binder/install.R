## General dependencies:
install.packages("seqinr")
install.packages("stringr")
install.packages("devtools")
install.packages("shinyalert")
install.packages("rlist")
install.packages("RColorBrewer")
install.packages("BiocManager")
install.packages("plotly")
install.packages("DT")
install.packages("shinyjs")
install.packages("rmarkdown")
install.packages("knitr")
## BioC dependencies

install.packages("BiocManager")
BiocManager::install(c("Biostrings", "sangerseqR"))
devtools::install_github('jbryer/DTedit')
devtools::install_github("ipb-halle/GoldenMutagenesis")
