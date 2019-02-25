#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GoldenMutagenesis)
library(seqinr)
library(stringr)
library(shinyalert)
library(rlist)
library(RColorBrewer)
library(graphics)

sequence_check<-function(input_sequence){
  input_sequence<-str_to_upper(input_sequence)
  if(nchar(input_sequence)%%3!=0) {
    stop(paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "))
  }
  codon_seq<-splitseq(s2c(input_sequence))
  met<-which(str_detect(codon_seq, "ATG"))
  if(length(met) == 0) {
    stop("No Methionine in the provided sequence. Stopping here. Please check the provided sequence.")
  }
  
  if(min(met) != 1){
    warning(paste("No Methionine at first codon found! Please check the provided sequence! Took codon #", min(met), "as start.", sep=" "))
    codon_seq<-codon_seq[min(met):length(codon_seq)]
  } #else(codon_seq<-codon_seq[-1])
  
  stop<-which(str_detect(codon_seq, "(TAA)|(TGA)|(TAG)"))
  if(length(stop) == 0) {
    stop("No stop codon in the provided sequence. Stopping here. Please check the provided sequence!")
  }
  
  if(max(stop) != length(codon_seq)) {
    warning(paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "))  
    codon_seq<-codon_seq[1:max(stop)]
  }# else {
  #codon_seq <- codon_seq[-length(codon_seq)]
  #}
  return(codon_seq)
}

print_sequence_dom<-function(sequence, mutations) {
  codon_seq<-sequence_check(sequence)
  aa_seq<-translate(s2c(sequence))
  renderlist<-list()
  positions_aa<-c()
  if(length(mutations)>0) {
    positions_aa<-c()
    names_aa<-c()
    for(i in 1:length(mutations)) {
      position_aa<-as.numeric(mutations[[i]][1])
      name_aa<-mutations[[i]][2]
      positions_aa<-c(positions_aa, position_aa)
      names_aa<-c(names_aa, name_aa)
    }
  }
  #print(positions_aa)
  renderlist[[1]]<-HTML("<table style=\"    
    border-collapse: separate;
    border-spacing: 5px 10px;
    font-family: monospace;\"><tr>")
  for (i in 1:length(aa_seq)) {
    if(i %in% positions_aa) {
      ind<-which(positions_aa == i)
      if(names_aa[ind]==aa_seq[positions_aa[ind]]){
        renderlist[[length(renderlist)+1]]<-HTML(paste("<td align=\"center\" bgcolor=\"#fcf81e\">", i, "<br>", aa_seq[i], "<br>", codon_seq[i],"</td>" , sep=""))
      } else {
        renderlist[[length(renderlist)+1]]<-HTML(paste("<td align=\"center\" bgcolor=\"#ff9999\">", i, "<br>", aa_seq[i], "<br>", codon_seq[i],"</td>" , sep=""))
      }
        
    } else {
      renderlist[[length(renderlist)+1]]<-HTML(paste("<td align=\"center\">", i, "<br>", aa_seq[i], "<br>", codon_seq[i], "</td>" , sep=""))
    }
    if((i %% 20 == 0) & (i != length(aa_seq))){
      renderlist[[length(renderlist)+1]]<-HTML(paste("</tr>","<tr>",sep=""))
    }
  }
  renderlist[[length(renderlist)+1]]<-HTML("</tr><table>")
  return(renderlist)
}

print_primer_fancy<-function(primerset) {
  renderlist<-list()
  for(i in 1:length(primerset@fragments)){
    panel<-wellPanel(style="background: #e3e7ff78",
    fluidRow(h2(paste("Fragment ", i, sep=""))),
    fluidRow(column(2, h5("Start"), primerset@fragments[[i]]@start), column(2, h5("Stop"),primerset@fragments[[i]]@stop), column(2, h5("Length"), (primerset@fragments[[i]]@stop - primerset@fragments[[i]]@start)+1)),
    fluidRow(h5("Primer Forward"), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", primerset@primers[[i]][[1]]@prefix, "</span>",
                                                         "<span style=\"background-color: #fc9191; font-size: large;\">", primerset@primers[[i]][[1]]@restriction_enzyme, "</span>",
                                                         "<span style=\"background-color: #e2d544; font-size: large;\">", primerset@primers[[i]][[1]]@suffix, "</span>",
                                                         "<span style=\"background-color: #9fff8e; font-size: large;\">", primerset@primers[[i]][[1]]@vector, "</span>",
                                                         "<span style=\"background-color: #76fcb7; font-size: large;\">", primerset@primers[[i]][[1]]@overhang, "</span>",
                                                         "<span style=\"background-color: #f9ffd1; font-size: large;\">", primerset@primers[[i]][[1]]@extra, "</span>",
                                                         "<span style=\"background-color: #a5f1ff; font-size: large;\">", primerset@primers[[i]][[1]]@binding_sequence, "</span>",
                                                         sep="")))),
    fluidRow(column(4, h5("Melting temperature binding site"), round(primerset@primers[[i]][[1]]@temperature, digits=3)), column(3, h5("Temperature difference to setting"), round(primerset@primers[[i]][[1]]@difference, digits=3))),
    fluidRow(h5("Primer Reverse"), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", primerset@primers[[i]][[2]]@prefix, "</span>",
                                                         "<span style=\"background-color: #fc9191; font-size: large;\">", primerset@primers[[i]][[2]]@restriction_enzyme, "</span>",
                                                         "<span style=\"background-color: #e2d544; font-size: large;\">", primerset@primers[[i]][[2]]@suffix, "</span>",
                                                         "<span style=\"background-color: #9fff8e; font-size: large;\">", primerset@primers[[i]][[2]]@vector, "</span>",
                                                         "<span style=\"background-color: #76fcb7; font-size: large;\">", primerset@primers[[i]][[2]]@overhang, "</span>",
                                                         "<span style=\"background-color: #f9ffd1; font-size: large;\">", primerset@primers[[i]][[2]]@extra, "</span>",
                                                         "<span style=\"background-color: #a5f1ff; font-size: large;\">", primerset@primers[[i]][[2]]@binding_sequence, "</span>",
                                                         sep="")))),
    fluidRow(column(4, h5("Melting temperature binding site"), round(primerset@primers[[i]][[2]]@temperature, digits = 3)), column(3, h5("Temperature difference to forward primer"), round(primerset@primers[[i]][[2]]@difference, digits = 3)))
    )
    renderlist<-list.append(renderlist, panel)
    
    #cat("Length ",(primerset@fragments[[i]]@stop - primerset@fragments[[i]]@start)+1, "\n", sep="")
    #cat("Forward\n")
    #print_primer(primerset@primers[[i]][[1]])
    #cat("Reverse\n")
    #print_primer(primerset@primers[[i]][[2]])
    #cat("\n")
  }
  return(renderlist)
}

base_distribution_shiny<-function(input_sequence, ab1file, replacements, trace_cutoff=80){
  plotlist<-c()
  sanger_seq<-sangerseqR::readsangerseq(ab1file) #reading in the data
  global_Align<-Biostrings::pairwiseAlignment(input_sequence, sanger_seq@primarySeq)
  global_Align_rev<-Biostrings::pairwiseAlignment(input_sequence, Biostrings::reverseComplement(sanger_seq@primarySeq))
  reverse=F
  if(global_Align_rev@score > global_Align@score) {
    reverse=T
    global_Align<-global_Align_rev
    print("Reverse sequence detected!")
  }
  mismatches<-Biostrings::mismatchTable(global_Align)
  replacements_basepairs<-as.vector(sapply(replacements, FUN<-function(x){return(c(x*3-2, x*3-1, x*3))}))
  candidates<-unlist(sapply(replacements_basepairs, FUN = function(x){which(mismatches[,"PatternStart"]==x)}, simplify = array))
  mismatches_candidates<-mismatches[candidates, ]
  mismatches_candidates$pos<-mismatches_candidates[,"PatternStart"]%%3
  mismatches_candidates[mismatches_candidates["pos"]==0, "pos"]<-3
  subject_pos<-vector()
  pattern_pos<-vector()
  for (i in 1:nrow(mismatches_candidates)) {
    subject_start<-mismatches_candidates[i, "SubjectStart"]
    pos<-mismatches_candidates[i, "pos"]
    pattern_start<-mismatches_candidates[i, "PatternStart"]
    if(pos==1) {
      subject_pos<-c(subject_pos, subject_start, subject_start+1, subject_start+2)
      pattern_pos<-c(pattern_pos, pattern_start, pattern_start+1, pattern_start+2)
      
    }
    if(pos==2) {
      subject_pos<-c(subject_pos, subject_start-1, subject_start, subject_start+1)
      pattern_pos<-c(pattern_pos, pattern_start-1, pattern_start, pattern_start+1)
      
    }
    if(pos==3) {
      subject_pos<-c(subject_pos, subject_start-2, subject_start-1, subject_start)
      pattern_pos<-c(pattern_pos, pattern_start-2, pattern_start-1, pattern_start)
      
    }
  }
  subject_pos<-unique(subject_pos)
  pattern_pos<-unique(pattern_pos)
  
  if(reverse==T) {
    subject_pos<-length(sanger_seq@primarySeq)-subject_pos+1
  }
  tracematrix_subject<-sangerseqR::traceMatrix(sanger_seq)[sangerseqR::peakPosMatrix(sanger_seq)[subject_pos],]
  sums_row<-which(rowSums(tracematrix_subject)>=trace_cutoff)
  tracematrix_subject<-as.data.frame(tracematrix_subject[sums_row,])
  for(element in sums_row) {
    # plotting as pie chart
    sliceit <- dplyr::slice (tracematrix_subject,element)
    slices <- as.numeric(sliceit)
    lbls <- c("Adenine", "Cytosine", "Guanine", "Thymine")
    if(reverse==T) {
      lbls <- c("Thymine", "Guanine", "Cytosine", "Adenine")
    }
    pct <- round(slices/sum(slices)*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    print(element)
    file<-tempfile(fileext = ".png")
    plotlist<-c(plotlist,file)
    png(file)
      pie(slices,labels = lbls, col=brewer.pal(4,"Spectral"),main = paste("Peak intensity distribution for \nPosition", pattern_pos[element], "(Template) -", subject_pos[element], "(Sequencing)", sep=" "))
    dev.off()
  }
  print(plotlist)
  return(plotlist)
}



ui <- fluidPage(
  useShinyalert(),
   # Application title
   titlePanel("GoldenMutagenesis Webtool Beta"),
   navlistPanel(id="MainNav", widths = c(4, 7),
     tabPanel("Welcome", h4("Welcome to the GoldenMutagenesis Webtool."), h4("Please select the desired task.")),
     "Pre- and Postprocessing",
        tabPanel("Domestication", 
          mainPanel(width = 15,
            tabsetPanel(id="d_s_p",type="tabs",
                tabPanel(
                "Sequence Input",
                br(), 
                p("You can paste a sequence into the textbox or upload your own fasta file."),
                tabsetPanel(type="pills",
                    tabPanel("Manual Input",
                            textAreaInput("d_seq", h4("Enter your sequence"), cols=60, rows = 10, resize = "both",
                              placeholder = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
                    ),
                    tabPanel("FASTA Upload",
                             p("If your fasta file has more than one sequence, only the first one will be used."),
                             fileInput("d_file", h4("Select .FA/.FASTA"), multiple = FALSE, accept = c(".fa", ".fasta"), width = NULL)
                    ),
                    fluidRow(column(2,actionButton(inputId = "d_n_i", 'Next')))
                )),
                tabPanel("Configuration", h1("Mutagenesis Configuration") ,  uiOutput("mut_conf") , 
                        fluidRow(h1("Domestication Configuration"), 
                                 wellPanel(style = "background: #ddffdd", p("Please select the restriction enzymes you want to domesticate."),
                                           fluidRow(column(6,
                                           fluidRow(column(3,checkboxInput("d_bsai", "BsaI",value = T)), column(3,checkboxInput("d_bbsi", "BbsI", value = T))),
                                           fluidRow(column(10,textInput("d_cu", "Custom Recognition Site"))),
                                           fluidRow(column(2,actionButton("d_add_cu", "Add")), column(2, actionButton("d_rm_cu", "Remove")))),
                                           column(6,h4("Restriction Sites"), uiOutput("d_cu_list"))
                                 ))
                                 ),
                        fluidRow(column(2,actionButton(inputId = "d_n_c", 'Next'), br())
                   ) 
                ),
                tabPanel("Preview and Selection",
                         wellPanel(uiOutput("preview")), wellPanel(style="background: #b7fff0",
                          fluidRow(
                             column(2,uiOutput("codonnum"))
                          ),
                          fluidRow(
                            column(2,p("Aminoacid:"),uiOutput("d_aa")),
                            column(2,checkboxInput("d_sm", label = "Silent Mutation?",value = F)),
                            column(2, actionButton("d_sm_apply", "Apply"))
                          )
                         ),                         
                         fluidRow(column(2,actionButton(inputId = "d_n_p", 'Next'), br()))
                         ),
                tabPanel("Results", 
                         fluidRow(column(12, h2("Legend")), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", "Prefix", "</span>",
                                                                              "<span style=\"background-color: #fc9191; font-size: large;\">", "Restriction Enzyme", "</span>",
                                                                              "<span style=\"background-color: #e2d544; font-size: large;\">", "Suffix", "</span>",
                                                                              "<span style=\"background-color: #9fff8e; font-size: large;\">", "Vector", "</span>",
                                                                              "<span style=\"background-color: #76fcb7; font-size: large;\">", "Overhang", "</span>",
                                                                              "<span style=\"background-color: #f9ffd1; font-size: large;\">", "Extra", "</span>",
                                                                              "<span style=\"background-color: #a5f1ff; font-size: large;\">", "Binding Sequence", "</span>",
                                                                              sep="")))),br(),fluidRow(column(5, actionButton("d_r_sp", "PointMutagenesis: Keep Mutations and Sequence")),column(3, actionButton("d_r_sps", "PointMutagenesis: Keep Sequence")), column(3, actionButton("d_r_ms", "SaturationMutagenesis: Keep Sequence"))),uiOutput("dom_primers"))
                )
          )),
     tabPanel("Quick-Quality-Control", 
              tabsetPanel(id="qqc", type="tabs", tabPanel("Configuration",
                          wellPanel(style="background: #c9ffe6", textAreaInput("qqc_seq", h4("Enter your sequence"), cols=60, rows = 10, resize = "both",
                                                                                            placeholder = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
                                                 ), 
                          wellPanel(style="background: #75ffbf", fileInput("qqc_file", h4("Select .ab1/.abf"), multiple = FALSE, accept = c(".ab1", ".abf"), width = NULL)),
                          wellPanel(style="background: #1df28f", h4("Mutation Positions"), fluidRow(column(3,numericInput("qqc_pos", "Sequence Position", value=1)), column(1,actionButton("qqc_add","Add")),column(1,actionButton("qqc_remove","Remove"))),
                                    fluidRow(uiOutput("qqc_mut_display"))),
                          fluidRow(actionButton("qqc_next", "Next"))
              ), tabPanel("Results", uiOutput("plotout"))
              )
            ),
     "Mutagenesis",
     tabPanel("Point-Mutagenesis", 
              mainPanel(width = 15,
                        tabsetPanel(id="s_p",type="tabs",
                                    tabPanel(
                                      "Sequence Input",
                                      br(), 
                                      p("You can paste a sequence into the textbox or upload your own fasta file."),
                                      tabsetPanel(type="pills",
                                                  tabPanel("Manual Input",
                                                           textAreaInput("sp_seq", h4("Enter your sequence"), cols=60, rows = 10, resize = "both",
                                                                         placeholder = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
                                                  ),
                                                  tabPanel("FASTA Upload",
                                                           p("If your fasta file has more than one sequence, only the first one will be used."),
                                                           fileInput("sp_file", h4("Select .FA/.FASTA"), multiple = FALSE, accept = c(".fa", ".fasta"), width = NULL)
                                                  ),
                                                  fluidRow(column(2,actionButton(inputId = "sp_n_i", 'Next')))
                                        )
                                    ),                 
                                    tabPanel("Configuration", h1("Mutagenesis Configuration") , uiOutput("sp_mut_conf"),
                                            fluidRow(column(2,actionButton(inputId = "sp_n_c", 'Next'), br()))
                                      ),
                                    tabPanel("Preview and Selection",
                                             wellPanel(uiOutput("sp_preview")), wellPanel(style="background: #b7fff0",
                                                                                       fluidRow(
                                                                                         column(2,uiOutput("sp_codonnum")),
                                                                                         column(2,uiOutput("sp_m_o"))
                                                                                       ),
                                                                                       fluidRow(
                                                                                         column(2,p("Aminoacid:"),uiOutput("sp_aa")),
                                                                                         column(2, actionButton("sp_m_apply", "Apply"))
                                                                                       )
                                             ),                         
                                             fluidRow(column(2,actionButton(inputId = "sp_n_p", 'Next'), br()))
                                    ),
                                    tabPanel("Results", 
                                             fluidRow(column(12, h2("Legend")), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", "Prefix", "</span>",
                                                                                                      "<span style=\"background-color: #fc9191; font-size: large;\">", "Restriction Enzyme", "</span>",
                                                                                                      "<span style=\"background-color: #e2d544; font-size: large;\">", "Suffix", "</span>",
                                                                                                      "<span style=\"background-color: #9fff8e; font-size: large;\">", "Vector", "</span>",
                                                                                                      "<span style=\"background-color: #76fcb7; font-size: large;\">", "Overhang", "</span>",
                                                                                                      "<span style=\"background-color: #f9ffd1; font-size: large;\">", "Extra", "</span>",
                                                                                                      "<span style=\"background-color: #a5f1ff; font-size: large;\">", "Binding Sequence", "</span>",
                                                                                                      sep="")))),br(),uiOutput("sp_primers"))
                                    
                                    )
                        )
              ),
     tabPanel("Saturation-Mutagenesis", mainPanel(width = 15,
                                                  tabsetPanel(id="ms_p",type="tabs",
                                                              tabPanel(
                                                                "Sequence Input",
                                                                br(), 
                                                                p("You can paste a sequence into the textbox or upload your own fasta file."),
                                                                tabsetPanel(type="pills",
                                                                            tabPanel("Manual Input",
                                                                                     textAreaInput("ms_seq", h4("Enter your sequence"), cols=60, rows = 10, resize = "both",
                                                                                                   placeholder = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
                                                                            ),
                                                                            tabPanel("FASTA Upload",
                                                                                     p("If your fasta file has more than one sequence, only the first one will be used."),
                                                                                     fileInput("ms_file", h4("Select .FA/.FASTA"), multiple = FALSE, accept = c(".fa", ".fasta"), width = NULL)
                                                                            ),
                                                                            fluidRow(column(2,actionButton(inputId = "ms_n_i", 'Next')))
                                                                )
                                                              ),
                                                              tabPanel("Configuration", h1("Mutagenesis Configuration") , uiOutput("ms_mut_conf"),
                                                                        fluidRow(column(2,actionButton(inputId = "ms_n_c", 'Next'), br()))
                                                              ),
                                                              tabPanel("Preview and Selection",
                                                                       wellPanel(uiOutput("ms_preview")), wellPanel(style="background: #b7fff0",
                                                                                                                    fluidRow(
                                                                                                                      column(2,uiOutput("ms_codonnum")),
                                                                                                                      column(2,uiOutput("ms_m_o"))
                                                                                                                    ),
                                                                                                                    fluidRow(
                                                                                                                      column(2,p("Aminoacid:"),uiOutput("ms_aa")),
                                                                                                                      column(2, actionButton("ms_m_apply", "Apply"))
                                                                                                                    )
                                                                       ),                         
                                                                       fluidRow(column(2,actionButton(inputId = "ms_n_p", 'Next'), br()))
                                                              ),
                                                              tabPanel("Results", 
                                                                       fluidRow(column(12, h2("Legend")), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", "Prefix", "</span>",
                                                                                                                                "<span style=\"background-color: #fc9191; font-size: large;\">", "Restriction Enzyme", "</span>",
                                                                                                                                "<span style=\"background-color: #e2d544; font-size: large;\">", "Suffix", "</span>",
                                                                                                                                "<span style=\"background-color: #9fff8e; font-size: large;\">", "Vector", "</span>",
                                                                                                                                "<span style=\"background-color: #76fcb7; font-size: large;\">", "Overhang", "</span>",
                                                                                                                                "<span style=\"background-color: #f9ffd1; font-size: large;\">", "Extra", "</span>",
                                                                                                                                "<span style=\"background-color: #a5f1ff; font-size: large;\">", "Binding Sequence", "</span>",
                                                                                                                                sep="")))),br(),uiOutput("ms_primers"))
                                                              
                                                    
                                                  )
                                        )
     )
  )

)

# Define server logic
server <- function(input, output, session) {
  rv<-reactiveValues()
  rv$restriction_sites<-c()
  rv$sp_mutations<-list()
  #########APPLY CONFIG TEMPLATES############
  ####TEMPLATE SELECTION#####
  generic_template_selection<-function(prefix) {
    observeEvent(input[[paste(prefix, "template", sep="_")]], {
      if(input[[paste(prefix, "template", sep="_")]] == "1") {
        print("pAGM9121 selected")
        updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv0")
        updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bbsi")
        updateTextInput(session, paste(sep="_", prefix, "v1"), value = "CTCA")
        updateTextInput(session, paste(sep="_", prefix, "v2"), value = "CTCG")
        updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv")
      }
      if(input[[paste(prefix, "template", sep="_")]] == "2") {
        print("pAGM22082 selected")
        updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
        updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
        updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
        updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
        updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv")        
      }
      if(input[[paste(prefix, "template", sep="_")]] == "3") {
        print("pICH86988 selected")
        updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
        updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
        updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
        updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
        updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv")
      }
      if(input[[paste(prefix, "template", sep="_")]] == "c") {
        print("custom selected")
      }
    })}
  ####RESTRICTION ENZYME SELECTION#####
  generic_re_selection<-function(prefix) {observeEvent(input[[paste(prefix,"re_enzyme_selection",sep="_")]], {
    if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bbsi") {
      updateTextInput(session, paste(prefix,"re_enzyme",sep="_"), value = "GAAGAC")
    }
    if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bsai") {
      updateTextInput(session,  paste(prefix,"re_enzyme",sep="_"), value = "GGTCTC")
    }
  })
  }
  generic_re_selection_lvl0<-function(prefix) {observeEvent(input[[paste(prefix,"lvl0_re_enzyme_selection",sep="_")]], {
    if(input[[paste(prefix, "level", sep="_")]] == "lv2") {
    if(input[[paste(prefix,"lvl0_re_enzyme_selection",sep="_")]] == "bbsi") {
      updateTextInput(session, paste(prefix,"lvl0_re_enzyme",sep="_"), value = "GAAGAC")
      print("bbsi")
    }
    if(input[[paste(prefix,"lvl0_re_enzyme_selection",sep="_")]] == "bsai") {
      updateTextInput(session,  paste(prefix,"lvl0_re_enzyme",sep="_"), value = "GGTCTC")
    }
    }
  })
  }
  
  ####OLD: REPLACE ME ! TEMPLATE SELECTION#####
  observeEvent(input$template, {
    if(input$template == "1") {
      print("pAGM9121 selected")
      updateSelectInput(session, "level", selected = "lv0")
      updateSelectInput(session, "re_enzyme_selection", selected = "bbsi")
      updateTextInput(session, "v1", value = "CTCA")
      updateTextInput(session, "v2", value = "CTCG")
      updateSelectInput(session, "cuf", selected =  "e_coli_316407.csv")
    }
    if(input$template == "2") {
      print("pAGM22082 selected")
      updateSelectInput(session, "level", selected = "lv2")
      updateSelectInput(session, "re_enzyme_selection", selected = "bsai")
      updateTextInput(session, "v1", value = "AATG")
      updateTextInput(session, "v2", value = "AAGC")
      updateSelectInput(session, "cuf", selected =  "e_coli_316407.csv")        
    }
    if(input$template == "3") {
      print("pICH86988 selected")
      updateSelectInput(session, "level", selected = "lv2")
      updateSelectInput(session, "re_enzyme_selection", selected = "bsai")
      updateTextInput(session, "v1", value = "AATG")
      updateTextInput(session, "v2", value = "AAGC")
      updateSelectInput(session, "cuf", selected =  "e_coli_316407.csv")
    }
    if(input$template == "c") {
      print("custom selected")
    }
  })
  ####RESTRICTION ENZYME SELECTION#####
  observeEvent(input$re_enzyme_selection, {
    if(input$re_enzyme_selection == "bbsi") {
      updateTextInput(session, "re_enzyme", value = "GAAGAC")
    }
    if(input$re_enzyme_selection == "bsai") {
      updateTextInput(session, "re_enzyme", value = "GGTCTC")
    }
  })
  observeEvent(input$lvl0_re_enzyme_selection, {
    if(input$lvl0_re_enzyme_selection == "bbsi") {
      updateTextInput(session, "lvl0_re_enzyme", value = "GAAGAC")
    }
    if(input$lvl0_re_enzyme_selection == "bsai") {
      updateTextInput(session, "lvl0_re_enzyme", value = "GGTCTC")
    }
  })
  observeEvent(input$d_bsai,{
    if(input$d_bsai == T){
      rv$restriction_sites<-union(rv$restriction_sites, "GGTCTC")
    }
    else {
      rv$restriction_sites<-setdiff(rv$restriction_sites, "GGTCTC")
    }
    output$d_cu_list<-renderUI({HTML(paste(rv$restriction_sites, collapse="<br>"))})
  })
  observeEvent(input$d_bbsi,{
    if(input$d_bbsi == T){
      rv$restriction_sites<-union(rv$restriction_sites, "GAAGAC")
    }
    else {
      rv$restriction_sites<-setdiff(rv$restriction_sites, "GAAGAC")
    }
    output$d_cu_list<-renderUI({HTML(paste(rv$restriction_sites, collapse="<br>"))})
  })
  
  observeEvent(input$d_add_cu, {
    rv$restriction_sites<-union(rv$restriction_sites, input$d_cu)
    output$d_cu_list<-renderUI({HTML(paste(rv$restriction_sites, collapse="<br>"))})
  })
  observeEvent(input$d_rm_cu, {
    rv$restriction_sites<-setdiff(rv$restriction_sites, input$d_cu)
    output$d_cu_list<-renderUI({HTML(paste(rv$restriction_sites, collapse="<br>"))})
  })
  
###########################################
  #PRINT PAGES#
  
#############CONFIG###############
  mut_conf<-function(prefix){renderUI({tagList(
  br(),
  fluidRow(
    column(6,
           wellPanel(style = "background: #cce6ff",
                     h2("Primer Configuration"),
                     p("You can select a preconfigured template to set those settings in accordiance to your Golden Gate Mutagenesis."),
                     selectInput(paste(prefix,"template", sep="_"), "Pre-existing configuration:", c("pAGM9121"="1", 
                                                                              "pAGM22082 Red" = "2",
                                                                              "pICH86988" = "3",
                                                                              "custom" = "c")),
                     selectInput(paste(prefix,"level",sep="_"), "Golden Gate Level:", c("Level0" = "lv0", "Level2" = "lv2")),
                     selectInput(paste(prefix,"re_enzyme_selection", sep="_"), "Restriction Enzyme:", c("BbsI"="bbsi",
                                                                                 "BsaI"="bsai", 
                                                                                 "custom" = "c")),
                     textInput(paste(prefix,"re_enzyme", sep="_"), "Restriction Enzyme Recognition Site", value = "GAAGAC"),
                     textInput(paste(prefix, "prefix", sep="_"),  "Prefix", value = "TT"),
                     textInput(paste(prefix, "suffix",sep="_"), "Suffix", value = "AA"),
                     fluidRow(column(6, textInput(paste(prefix,"v1",sep="_"), "Vector Forward 5'", value = "CTCA")),column(6, textInput(paste(prefix,"v2",sep="_"), "Vector Reverse 3'", value = "CTCG")))
           )
    ),
    column(6,
           wellPanel(style = "background: #ffcccc",
                     h2("Algorithm Settings"),
                     p("Those settings are documentated in the GoldenMutagenesis R-Package"),
                     numericInput(paste(prefix,"binding_min_length",sep="_"), 
                                  "Minimal binding length (AA)", 
                                  value = 4),
                     numericInput(paste(prefix, "primer_length",sep="_"), 
                                  "Maximal binding sequence length (AA)", 
                                  value = 9),
                     numericInput(paste(prefix,"temperature",sep="_"), 
                                  "Target Temperature in Celsius", 
                                  value = 60),
                     selectInput(paste(prefix,"cuf",sep="_"), "Codon Usage Table", list_cu_table())
           )
    )),
  fluidRow(column(12,wellPanel(style = "background: #f7ff89", h2("Golden Gate Level Settings"), uiOutput(paste(prefix,"levelsettings", sep="_"))))
  ))
  })}
  
  levelsettings<-function(prefix){
  level<-reactive({input[[paste(prefix, "level", sep="_")]]})
  output[[paste(prefix,"levelsettings", sep="_")]]<-renderUI({
    if(level() == "lv0") {
      tagList (
        fluidRow(checkboxInput(paste(prefix,"prepare_lvl2", sep="_"), "Prepare for use in Level2", value = TRUE, width = NULL)),
        fluidRow(column(6, textInput(paste(prefix,"av1", sep="_"), "Accepting Vector Forward 5'", value = "AATG")),column(6, textInput(paste(prefix,"av2",sep="_"), "Accepting Vector Reverse 3'", value = "AAGC")))
      )
    }
    else{
      tagList(
        fluidRow(column(5,checkboxInput(paste(prefix,"add_lvl0",sep="_"), "Use Level0 to go in Level2?", value = TRUE, width = NULL))),
        fluidRow(column(6,textInput(paste(prefix,"lvl0_v1",sep="_"), "Level0 Vector Forward 5'", value = "CTCA")),column(6, textInput(paste(prefix,"lvl0_v2",sep="_"), "Level0 Vector Reverse 3'", value = "CTCG"))),
        fluidRow(column(6,selectInput(paste(prefix,"lvl0_re_enzyme_selection",sep="_"), "Restriction Enzyme:", c("BbsI"="bbsi", "BsaI"="bsai", "custom" = "c"))),
                 column(6,textInput(paste(prefix,"lvl0_re_enzyme",sep="_"), "Restriction Enzyme Sequence", value = "GAAGAC"))),
        fluidRow(column(6,textInput(paste(prefix,"lvl0_prefix",sep="_"), "Prefix", value = "TT")),
                 column(6,textInput(paste(prefix,"lvl0_suffix",sep="_"), "Suffix", value = "AA")))       
      )      
    }
  })}
  
  ####REPLACE ME AS SOON AS POSSIBLE !####
  output$mut_conf<-renderUI({tagList(
    br(),
    fluidRow(
      column(6,
             wellPanel(style = "background: #cce6ff",
                       h2("Primer Configuration"),
                       p("You can select a preconfigured template to set those settings in accordiance to your Golden Gate Mutagenesis."),
                       selectInput("template", "Pre-existing configuration:", c("pAGM9121"="1", 
                                                                                "pAGM22082 Red" = "2",
                                                                                "pICH86988" = "3",
                                                                                "custom" = "c")),
                       selectInput("level", "Golden Gate Level:", c("Level0" = "lv0", "Level2" = "lv2")),
                       selectInput("re_enzyme_selection", "Restriction Enzyme:", c("BbsI"="bbsi",
                                                                                   "BsaI"="bsai", 
                                                                                   "custom" = "c")),
                       textInput("re_enzyme", "Restriction Enzyme Recognition Site", value = "GAAGAC"),
                       textInput("prefix", "Prefix", value = "TT"),
                       textInput("suffix", "Suffix", value = "AA"),
                       fluidRow(column(6, textInput("v1", "Vector Forward 5'", value = "CTCA")),column(6, textInput("v2", "Vector Reverse 3'", value = "CTCG")))
             )
      ),
      column(6,
             wellPanel(style = "background: #ffcccc",
                       h2("Algorithm Settings"),
                       p("Those settings are documentated in the GoldenMutagenesis R-Package"),
                       numericInput("binding_min_length", 
                                    "Minimal binding length (AA)", 
                                    value = 4),
                       numericInput("primer_length", 
                                    "Maximal binding sequence length (AA)", 
                                    value = 9),
                       numericInput("temperature", 
                                    "Target Temperature in Celsius", 
                                    value = 60),
                       selectInput("cuf", "Codon Usage Table", list_cu_table())
             )
      )),
    fluidRow(column(12,wellPanel(style = "background: #f7ff89", h2("Golden Gate Level Settings"), uiOutput("levelsettings")))
             
    ))
  })
###########################################
  
  ####################################
  ###########DOMESTICATION#############
  #####################################
  #########INPUT##########
  observeEvent(input$d_n_i, {
    if(input$d_seq == "") {
      shinyalert("No Sequence!", "You have not entered a sequence. The default value will be used!", type = "warning")
      updateTextAreaInput(session, "d_seq", value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
    }
    updateTabsetPanel(session, "d_s_p", "Configuration")
  })
  #########Processing##########
  observeEvent(input$d_n_c, {
      if(length(rv$restriction_sites)==0) {
        shinyalert("No Recognition Site!", "You have not selected any recognition site for domestication. Nothing to do!", type = "error")
    } else {
      updateTabsetPanel(session, "d_s_p", "Preview and Selection")
      rv$dom_mutations<-list()
      for (rs in rv$restriction_sites) {
        mutations<-domesticate(input$d_seq, rs, cuf = input$cuf)
        if(length(mutations)>0){
          rv$dom_mutations<-c(rv$dom_mutations, mutations)
        }
      }
      output$preview<-renderUI(print_sequence_dom(sequence = input$d_seq, mutations = rv$dom_mutations))
      }
  })
  #########CONFIGURATION##########
  level<-reactive({input$level})
    output$levelsettings<- renderUI({
    if(level() == "lv0") {
    tagList (
      fluidRow(checkboxInput("prepare_lvl2", "Prepare for use in Level2", value = TRUE, width = NULL)),
      fluidRow(column(6, textInput("av1", "Accepting Vector Forward 5'", value = "AATG")),column(6, textInput("av2", "Accepting Vector Reverse 3'", value = "AAGC")))
    )
    }
    else{
      tagList(
        fluidRow(column(5,checkboxInput("add_lvl0", "Use Level0 to go in Level2?", value = TRUE, width = NULL))),
        fluidRow(column(6, textInput("lvl0_v1", "Level0 Vector Forward 5'", value = "CTCA")),column(6, textInput("lvl0_v2", "Level0 Vector Reverse 3'", value = "CTCG"))),
        fluidRow(column(6,selectInput("lvl0_re_enzyme_selection", "Restriction Enzyme:", c("BbsI"="bbsi", "BsaI"="bsai", "custom" = "c"))),
        column(6,textInput("lvl0_re_enzyme", "Restriction Enzyme Sequence", value = "GAAGAC"))),
        fluidRow(column(6,textInput("lvl0_prefix", "Prefix", value = "TT")),
        column(6,textInput("lvl0_suffix", "Suffix", value = "AA")))       
      )      
    }
 })
 ############PREVIEW AND SELECTION#############
    output$codonnum<-renderUI(
      selectInput("d_codonpos", "Aminoacid Position", choices = 1:length(translate(s2c(input$d_seq))))
    )
    observeEvent(input$d_codonpos, {
      output$d_aa<-renderUI({
        HTML(aaa(translate(s2c(sequence_check(input$d_seq)[as.numeric(input$d_codonpos)]))))}
      )
      if(length(rv$dom_mutations)>0){
        positions_aa<-c()
        for(i in 1:length(rv$dom_mutations)) {
          position_aa<-as.numeric(rv$dom_mutations[[i]][1])
          positions_aa<-c(positions_aa, position_aa)
        }
        #print(positions_aa)
        if(as.numeric(input$d_codonpos) %in% positions_aa) {
          #print("check")
          updateCheckboxInput(session, "d_sm", value=T)
        }
        else {
          updateCheckboxInput(session, "d_sm", value=F)
        }
      }
    })
    observeEvent(input$d_sm_apply, {
      codon<-as.numeric(input$d_codonpos)
      positions_aa<-c()
      selected<-input$d_sm
      if(length(rv$dom_mutations)>0){
        for(i in 1:length(rv$dom_mutations)) {
          position_aa<-as.numeric(rv$dom_mutations[[i]][1])
          positions_aa<-c(positions_aa, position_aa)
        }
        if(selected==F){
          ind<-which(positions_aa==codon)
          if(length(ind) > 0) { 
            rv$dom_mutations<-list.remove(rv$dom_mutations, ind)
          }
        }
        if(selected==T){
          ind<-which(positions_aa==codon)
          if(length(ind) == 0) {
            rv$dom_mutations<-list.append(rv$dom_mutations, c(codon, translate(s2c(input$d_seq))[codon]))
          }
        }
      } else {
        if(selected==T){
            rv$dom_mutations<-list.append(rv$dom_mutations, c(codon, translate(s2c(input$d_seq))[codon]))
        }
      }
      print(rv$dom_mutations)
      output$preview<-renderUI(print_sequence_dom(sequence = input$d_seq, mutations = rv$dom_mutations))
    })
    
  ###########RESULTS################
    observeEvent(input$d_n_p, {
      if(length(rv$dom_mutations)==0) {
        shinyalert("No Mutations!", "You do not need to domesticate!", type = "info")
      }
      else{
        updateTabsetPanel(session, "d_s_p", "Results")
        rv$dom_primers<-mutate_spm(input$d_seq, prefix = input$prefix, restriction_enzyme = input$re_enzyme, suffix = input$suffix, vector = c(input$v1, input$v2), replacements = rv$dom_mutations,  binding_min_length = input$binding_min_length, target_temp = input$temperature, cuf = input$cuf, primer_length = input$primer_length)
        if(input$level=="lv0") {
          if(input$prepare_lvl2==TRUE) {
            print("prepare lv0")
            rv$dom_primers<-primer_prepare_level(rv$dom_primers, vector = c(input$av1, input$av2))
          }
        }
        if(input$level=="lv2") {
          if(input$add_lvl0 == TRUE){
            print("prepare lv2")
            rv$dom_primers<-primer_add_level(primerset = rv$dom_primers, prefix = input$lvl0_prefix,suffix =  input$lvl0_suffix, restriction_enzyme = lvl0_re_enzyme, vector = c(input$lvl0_v1, input$lvl0_v2))
          }
        }
        print(rv$dom_primers)
        output$dom_primers<-renderText(paste(capture.output(print_primer(rv$dom_primers)), collapse = "<br>"))
        output$dom_primers<-renderUI(print_primer_fancy(rv$dom_primers))
      }
    })
    observeEvent(input$d_r_sp, {
      updateTextInput(session, "sp_seq", value = rv$dom_primers@oldsequence)
      rv$sp_mutations<-rv$dom_mutations
      updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
      updateTabsetPanel(session, "s_p", selected = "Configuration")
    })
    observeEvent(input$d_r_sps, {
      updateTextInput(session, "sp_seq", value = rv$dom_primers@newsequence)
      updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
      updateTabsetPanel(session, "s_p", selected = "Configuration")
    })
    observeEvent(input$d_r_ms, {
      updateTextInput(session, "ms_seq", value = rv$dom_primers@newsequence)
      updateTabsetPanel(session, "MainNav", selected = "Saturation-Mutagenesis")
      updateTabsetPanel(session, "ms_p", selected = "Configuration")
    })
    ####################################
    ######Single Point Mutagenesis#######
    #####################################
    #########INPUT##########
    observeEvent(input$sp_n_i, {
      if(input$sp_seq == "") {
        shinyalert("No Sequence!", "You have not entered a sequence. The default value will be used!", type = "warning")
        updateTextAreaInput(session, "sp_seq", value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
      }
      updateTabsetPanel(session, "s_p", "Configuration")
    })
    ########Configuration#######
    output$sp_mut_conf<-mut_conf("sp")
    generic_template_selection("sp")
    generic_re_selection("sp")
    generic_re_selection_lvl0("sp")
    levelsettings("sp")
    observeEvent(input$sp_n_c, {
        updateTabsetPanel(session, "s_p", "Preview and Selection")
        output$sp_preview<-renderUI(print_sequence_dom(sequence = input$sp_seq, mutations = rv$sp_mutations))
    })
    ############PREVIEW AND SELECTION#############
    output$sp_codonnum<-renderUI(
      selectInput("sp_codonpos", "Aminoacid Position", choices = 1:length(translate(s2c(input$sp_seq))))
    )
    output$sp_m_o<-renderUI(
      selectInput("sp_m", "Aminoacid Replacement", choices = c("none", aaa()))
    )
    observeEvent(input$sp_codonpos, {
      output$sp_aa<-renderUI({
        HTML(aaa(translate(s2c(sequence_check(input$sp_seq)[as.numeric(input$sp_codonpos)]))))}
      )
      if(length(rv$sp_mutations)>0){
        positions_aa<-c()
        names_aa<-c()
        for(i in 1:length(rv$sp_mutations)) {
          position_aa<-as.numeric(rv$sp_mutations[[i]][1])
          name_aa<-rv$sp_mutations[[i]][2]
          positions_aa<-c(positions_aa, position_aa)
          names_aa<-c(names_aa, name_aa)
        }
        #print(positions_aa)
        if(as.numeric(input$sp_codonpos) %in% positions_aa) {
          #print("check")
          ind<-which(positions_aa == as.numeric(input$sp_codonpos))
          updateSelectInput(session, inputId = "sp_m", selected = aaa(names_aa[ind]))
        }
        else{
          updateSelectInput(session, inputId = "sp_m", selected = "none")
        }
      }
    })
    
     observeEvent(input$sp_m_apply, {
       codon<-as.numeric(input$sp_codonpos)
       positions_aa<-c()
       selected<-input$sp_m
       if(length(rv$sp_mutations)>0){
           positions_aa<-c()
           names_aa<-c()
           for(i in 1:length(rv$sp_mutations)) {
             position_aa<-as.numeric(rv$sp_mutations[[i]][1])
             name_aa<-rv$sp_mutations[[i]][2]
             positions_aa<-c(positions_aa, position_aa)
             names_aa<-c(names_aa, name_aa)
           }
         }
         if(selected=="none"){
           ind<-which(positions_aa==codon)
           if(length(ind) > 0) { 
             rv$sp_mutations<-list.remove(rv$sp_mutations, ind)
           }
         }
         if(selected!="none"){
           ind<-which(positions_aa==codon)
           if(length(ind) == 0) {
             rv$sp_mutations<-list.append(rv$sp_mutations, c(codon, a(selected)))
           }
           else{
             rv$sp_mutations<-list.remove(rv$sp_mutations, ind)
             rv$sp_mutations<-list.append(rv$sp_mutations, c(codon, a(selected)))
           }
         }
       print(rv$sp_mutations)
         output$sp_preview<-renderUI(print_sequence_dom(sequence = input$sp_seq, mutations = rv$sp_mutations))
     })
     ###########RESULTS################
     observeEvent(input$sp_n_p, {
       if(length(rv$sp_mutations)==0) {
         shinyalert("No Mutations!", "You should enter some mutations!", type = "info")
       }
       else{
         updateTabsetPanel(session, "s_p", "Results")
         rv$sp_primers<-mutate_spm(input$sp_seq, prefix = input$sp_prefix, restriction_enzyme = input$sp_re_enzyme, suffix = input$sp_suffix, vector = c(input$sp_v1, input$sp_v2), replacements = rv$sp_mutations,  binding_min_length = input$sp_binding_min_length, target_temp = input$sp_temperature, cuf = input$sp_cuf, primer_length = input$sp_primer_length)
         if(input$sp_level=="lv0") {
           if(input$sp_prepare_lvl2==TRUE) {
             print("prepare lv0")
             rv$sp_primers<-primer_prepare_level(rv$sp_primers, vector = c(input$sp_av1, input$sp_av2))
           }
         }
         if(input$sp_level=="lv2") {
           if(input$sp_add_lvl0 == TRUE){
             print("prepare lv2")
             rv$sp_primers<-primer_add_level(primerset = rv$sp_primers, prefix = input$sp_lvl0_prefix, suffix =  input$sp_lvl0_suffix, restriction_enzyme = input$sp_lvl0_re_enzyme, vector = c(input$sp_lvl0_v1, input$sp_lvl0_v2))
           }
         }
         print(rv$sp_primers)
         #output$sp_primers<-renderText(paste(capture.output(print_primer(rv$dom_primers)), collapse = "<br>"))
         output$sp_primers<-renderUI(print_primer_fancy(rv$sp_primers))
       }
     })
     ####################################
     ######Saturation Mutagenesis#######
     #####################################
     #########INPUT##########
     observeEvent(input$ms_n_i, {
       if(input$ms_seq == "") {
         shinyalert("No Sequence!", "You have not entered a sequence. The default value will be used!", type = "warning")
         updateTextAreaInput(session, "ms_seq", value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
       }
       updateTabsetPanel(session, "ms_p", "Configuration")
     })
     ########Configuration#######
     output$ms_mut_conf<-mut_conf("ms")
     generic_template_selection("ms")
     generic_re_selection("ms")
     generic_re_selection_lvl0("ms")
     levelsettings("ms")
     observeEvent(input$ms_n_c, {
       updateTabsetPanel(session, "ms_p", "Preview and Selection")
       rv$ms_mutations<-list()
       output$ms_preview<-renderUI(print_sequence_dom(sequence = input$ms_seq, mutations = rv$ms_mutations))
     })
     ############PREVIEW AND SELECTION#############
     output$ms_codonnum<-renderUI(
       selectInput("ms_codonpos", "Aminoacid Position", choices = 1:length(translate(s2c(input$ms_seq))))
     )
     output$ms_m_o<-renderUI(
       selectInput("ms_m", "Saturation ", choices = c("none", "NNN", "NNK", "NNS", "NDT", "DBK", "NRT"))
     )
     observeEvent(input$ms_codonpos, {
       output$ms_aa<-renderUI({
         HTML(aaa(translate(s2c(sequence_check(input$ms_seq)[as.numeric(input$ms_codonpos)]))))}
       )
       if(length(rv$ms_mutations)>0){
         positions_aa<-c()
         names_aa<-c()
         for(i in 1:length(rv$ms_mutations)) {
           position_aa<-as.numeric(rv$ms_mutations[[i]][1])
           name_aa<-rv$ms_mutations[[i]][2]
           positions_aa<-c(positions_aa, position_aa)
           names_aa<-c(names_aa, name_aa)
         }
         #print(positions_aa)
         if(as.numeric(input$ms_codonpos) %in% positions_aa) {
           #print("check")
           ind<-which(positions_aa == as.numeric(input$ms_codonpos))
           updateSelectInput(session, inputId = "ms_m", selected = names_aa[ind])
         }
         else{
           updateSelectInput(session, inputId = "ms_m", selected = "none")
         }
       }
     })
     
     observeEvent(input$ms_m_apply, {
       codon<-as.numeric(input$ms_codonpos)
       positions_aa<-c()
       selected<-input$ms_m
       if(length(rv$ms_mutations)>0){
         positions_aa<-c()
         names_aa<-c()
         for(i in 1:length(rv$ms_mutations)) {
           position_aa<-as.numeric(rv$ms_mutations[[i]][1])
           name_aa<-rv$ms_mutations[[i]][2]
           positions_aa<-c(positions_aa, position_aa)
           names_aa<-c(names_aa, name_aa)
         }
       }
       if(selected=="none"){
         ind<-which(positions_aa==codon)
         if(length(ind) > 0) { 
           rv$ms_mutations<-list.remove(rv$ms_mutations, ind)
         }
       }
       if(selected!="none"){
         ind<-which(positions_aa==codon)
         if(length(ind) == 0) {
           rv$ms_mutations<-list.append(rv$ms_mutations, c(codon, selected))
         }
         else{
           rv$ms_mutations<-list.remove(rv$ms_mutations, ind)
           rv$ms_mutations<-list.append(rv$ms_mutations, c(codon, selected))
         }
       }
       print(rv$ms_mutations)
       output$ms_preview<-renderUI(print_sequence_dom(sequence = input$ms_seq, mutations = rv$ms_mutations))
     })
     ###########RESULTS################
     observeEvent(input$ms_n_p, {
       if(length(rv$ms_mutations)==0) {
         shinyalert("No Mutations!", "You should enter some mutations!", type = "info")
       }
       else{
         updateTabsetPanel(session, "ms_p", "Results")
         rv$ms_primers<-msd_mutate(input$ms_seq, prefix = input$ms_prefix, restriction_enzyme = input$ms_re_enzyme, suffix = input$ms_suffix, vector = c(input$ms_v1, input$ms_v2), replacements = rv$ms_mutations,  binding_min_length = input$ms_binding_min_length, target_temp = input$ms_temperature, primer_length = input$ms_primer_length)
         if(input$ms_level=="lv0") {
           if(input$ms_prepare_lvl2==TRUE) {
             print("prepare lv0")
             rv$ms_primers<-primer_prepare_level(rv$ms_primers, vector = c(input$ms_av1, input$ms_av2))
           }
         }
         if(input$ms_level=="lv2") {
           if(input$ms_add_lvl0 == TRUE){
             print("prepare lv2")
             rv$ms_primers<-primer_add_level(primerset = rv$ms_primers, prefix = input$ms_lvl0_prefix, suffix =  input$ms_lvl0_suffix, restriction_enzyme = input$ms_lvl0_re_enzyme, vector = c(input$ms_lvl0_v1, input$ms_lvl0_v2))
           }
         }
         print(rv$ms_primers)
         #output$sp_primers<-renderText(paste(capture.output(print_primer(rv$dom_primers)), collapse = "<br>"))
         output$ms_primers<-renderUI(print_primer_fancy(rv$ms_primers))
       }
     })
     ####################################
     ######Quick Quality Control #######
     #####################################
     observeEvent(input$qqc_add, {
       rv$qqc_mutations<-union(rv$qqc_mutations, input$qqc_pos)
       print(paste(rv$qqc_mutations, collapse = " "))
       output$qqc_mut_display<-renderUI(paste(rv$qqc_mutations, collapse = " "))
     })
     
     observeEvent(input$qqc_remove, {
       rv$qqc_mutations<-setdiff(rv$qqc_mutations, input$qqc_pos)
       print(paste(rv$qqc_mutations, collapse = " "))
       output$qqc_mut_display<-renderUI(paste(rv$qqc_mutations, collapse = " "))
     })
     observeEvent(input$qqc_next,{
       updateTabsetPanel(session, "qqc", "Results")
       rv$plots<-base_distribution_shiny(input$qqc_seq, input$qqc_file$datapath, replacements = rv$qqc_mutations)
       output$plotout <- renderUI({
         image_output_list <- 
           lapply(1:length(rv$plots),
                  function(i)
                  {
                    imagename = paste0("image", i)
                    imageOutput(imagename)
                  })
         
         do.call(tagList, image_output_list)
       })
       observe({
         for (i in 1:length(rv$plots))
         {
           local({
             my_i <- i
             imagename = paste0("image", my_i)
             output[[imagename]] <- 
               renderImage({
                 list(src = rv$plots[my_i],
                      alt = "Image failed to render")
               }, deleteFile = FALSE)
           })
         }
       })
     })
     

}
# Run the application 
shinyApp(ui=ui, server=server)

