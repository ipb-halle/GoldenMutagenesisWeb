## General dependencies:
install.packages("seqinr")
install.packages("stringr")
install.packages("devtools")
install.packages("shinyalert")
install.packages("rlist")
install.packages("RColorBrewer")
install.packages("BiocManager")
## BioC dependencies

install.packages("BiocManager")
BiocManager::install(c("Biostrings", "sangerseqR"))

devtools::install_github("ipb-halle/GoldenMutagenesis@1.1.1")
