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
library(plotly)
source("helpers.R")
library("shiny")
library("DT")
library(shinyjs)
library("plyr")

## Read IPB coorporate identity
ipbheader <- HTML(readLines("ipbheader.html"))
ipbfooter <- HTML(readLines("ipbfooter.html"))

ui <- fluidPage(
  useShinyalert(),
  useShinyjs(),
  # IPB header 
  tags$head(tags$link(rel = "stylesheet", type = "text/css", 
                      href = "css/ipb-styles.css")),
  # Application title
  #titlePanel("GoldenMutagenesis Webtool Beta"),
  
  # Application title
  titlePanel(title=ipbheader, windowTitle = "GoldenMutagenesis"),
    navlistPanel(id="MainNav", widths = c(2,9), fluid = F,
     tabPanel("Welcome", h3("Welcome to the GoldenMutagenesis Webtool."),br(), h4("Please cite our publication"),
              h5("Golden Mutagenesis: An efficient multi-site-saturation mutagenesis approach by Golden Gate cloning with automated primer design"),
              h6("Pascal Püllmann, Chris Ulpinnis, Sylvestre Marillonnet, Ramona Gruetzner, Steffen Neumann & Martin J. Weissenborn "),
              br(),shiny::a("Publication as enhanced PDF on Springer Nature",href="https://rdcu.be/bMfta"),br(),br() ,p("AppVersion: 2019-07-29")), "Pre- and Postprocessing",
        tabPanel("Domestication", 
          wellPanel(width = 15,
            tabsetPanel(id="domestication",type="tabs",
                        tabPanel("Sequence Input", 
                                 wellPanel(style="background: #ffffff",uiOutput("domestication_input_panel"))),
                        tabPanel("Configuration", 
                                 wellPanel(style="background: #ffffff",
                                          fluidRow(h1("Domestication Configuration"), 
                                          wellPanel(style = "background: #ddffdd", p("Please select the restriction enzymes you want to domesticate."),
                                                    fluidRow(column(6,
                                                                    fluidRow(column(3,checkboxInput("d_bsai", "BsaI",value = T)), column(3,checkboxInput("d_bbsi", "BbsI", value = T))),br(),
                                                                    fluidRow(column(10,textInput("d_cu", "Custom Recognition Site"))),br(),
                                                                    fluidRow(column(2,actionButton("d_add_cu", "Add")), column(2, actionButton("d_rm_cu", "Remove")))),
                                                             column(6,wellPanel(h4("Selected Restriction Sites"), uiOutput("d_cu_list")))
                                                    )),
                                          h1("Mutagenesis Configuration"),uiOutput("domestication_mut_conf")
                                        )
                                  )
                ),
                tabPanel("Preview and Selection",
                         wellPanel(style="background: #ffffff",uiOutput("domestication_preview_complete"))),
                tabPanel("Results", 
                         wellPanel(style="background: #ffffff",uiOutput("domestication_primer_complete"),br()),
                         fluidRow(column(5, actionButton("domestication_result_spm", "PointMutagenesis: Keep mutations and old sequence"))),fluidRow(column(3, actionButton("domestication_result_spm_sequence", "PointMutagenesis: Keep new sequence"))), fluidRow(column(3, actionButton("domestication_result_msd_sequence", "SaturationMutagenesis: Keep new sequence")))
                         )
                )
          )),
     tabPanel("Quick-Quality-Control", 
              wellPanel(width = 15,
              tabsetPanel(id="qqc", type="tabs", tabPanel("Configuration",
                          wellPanel(style="background: #ffffff",  uiOutput("qqc_input_panel")),
                          wellPanel(style="background: #ffffff", fileInput("qqc_file", h4("Select .ab1/.abf"), multiple = FALSE, accept = c(".ab1", ".abf"), width = NULL)),
                          wellPanel(style="background: #ffffff", 
                                    h4("Mutation Positions"), 
                                    fluidRow(column(6, uiOutput("qqc_mutation_table"))),
                                    #fluidRow(column(3,numericInput("qqc_pos", "Sequence Position", value=1)), column(1,actionButton("qqc_add","Add")),column(1,actionButton("qqc_remove","Remove"))),
                                    fluidRow(uiOutput("qqc_mut_display"))),
                          fluidRow(column(3,actionButton("qqc_next", "Next")))
              ), tabPanel("Results", uiOutput("plots"))
              )
              )
            ),
     "Mutagenesis",
     tabPanel("Point-Mutagenesis", 
              wellPanel(width = 15,
                        tabsetPanel(id="spm",type="tabs",
                                    tabPanel("Sequence Input",
                                             wellPanel(style="background: #ffffff",uiOutput("spm_input_panel"))
                                    ),                 
                                    tabPanel("Configuration", wellPanel(style="background: #ffffff",h1("Mutagenesis Configuration") , uiOutput("spm_mut_conf"))
                                      ),
                                    tabPanel("Preview and Selection",
                                             wellPanel(style="background: #ffffff",uiOutput("spm_preview_complete"))
                                    ),
                                    tabPanel("Results",
                                             wellPanel(style="background: #ffffff", uiOutput("spm_primer_complete"))
                                    )
                        )
              )
      ),
     tabPanel("Saturation-Mutagenesis", wellPanel(width = 15,
                                                  tabsetPanel(id="msd",type="tabs",
                                                              tabPanel("Sequence Input",
                                                                       wellPanel(style="background: #ffffff",uiOutput("msd_input_panel"))
                                                              ),
                                                              tabPanel("Configuration", wellPanel(style="background: #ffffff",h1("Mutagenesis Configuration") , uiOutput("msd_mut_conf"))
                                                                       ),
                                                              tabPanel("Preview and Selection",
                                                                       wellPanel(style="background: #ffffff",uiOutput("msd_preview_complete"))
                                                                       ),
                                                              tabPanel("Results",
                                                                       wellPanel(style="background: #ffffff",uiOutput("msd_primer_complete"))
                                                                       )
                                                  )
                                        )
     )
  ),
  # Application footer
  fluidRow(ipbfooter)
)

# Define server logic
server <- function(input, output, session) {
  rv<-reactiveValues()
  rv$restriction_sites<-c()
  rv$sp_mutations<-list()
  rv$qqc_mutations_df<-data.frame()
  rv$qqc_mutations<-c()
  #########APPLY CONFIG TEMPLATES############
  source("server_templates.R", local=T)$value
  source("ui_templates.R", local = T)$value
 
 #####################################
 ###########DOMESTICATION#############
 #####################################
 #########INPUT##########
  generic_sequence_input("domestication", default_value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
  generic_process_input("domestication", default_value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
  generic_process_fasta_input("domestication")
  #########CONFIGURATION##########
  generic_mut_conf("domestication")
  generic_template_selection("domestication")
  generic_re_selection("domestication")
  generic_re_selection_lvl0("domestication")
  levelsettings("domestication")
  ###Off-target restriction site search configuration
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
  ###
  
  #########Processing##########
  observeEvent(input$domestication_configuration_next, {
    if(length(rv$restriction_sites)==0) {
      shinyalert("No Recognition Site!", "You have not selected any recognition site for domestication. Nothing to do!", type = "error")
    } else {
      updateTabsetPanel(session, "domestication", "Preview and Selection")
      rv$domestication_mutations<-list()
      for (rs in rv$restriction_sites) {
        mutations<-domesticate(rv$domestication_input_sequence, rs, cuf = input$domestication_cuf)
        if(length(mutations)>0){
          rv$domestication_mutations<-c(rv$domestication_mutations, mutations)
        }
      }
      generic_simple_preview_logic("domestication")
      generic_preview("domestication")  
    }
  })
  ############PREVIEW AND SELECTION#############
  generic_simple_selection("domestication")
  #generic_simple_preview_logic("domestication")
  
  ###########RESULTS################
  generic_results("domestication")
  generic_primer_output("domestication")
  #Does not work correclty, because UIOutput is triggered by clicking on the tabset. Therefore the TextAreaInput is not exisiting at this step. 
  observeEvent(input$domestication_result_spm, {
    updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
    #generic_sequence_input("spm")
    updateTextAreaInput(session, "spm_input_sequence", value = rv$domestication_primers@oldsequence)
    rv$spm_mutations<-rv$domestication_mutations
    updateTabsetPanel(session, "spm", selected = "Configuration")
    rv[["spm_input_sequence"]]<-rv$domestication_primers@oldsequence
  })
  observeEvent(input$domestication_result_spm_sequence, {
    generic_sequence_input("spm")
    updateTextAreaInput(session, "spm_input_sequence", value = as.character(rv$domestication_primers@newsequence))
    updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
    updateTabsetPanel(session, "spm", selected = "Configuration")
    rv[["spm_input_sequence"]]<-rv$domestication_primers@newsequence
  })
  observeEvent(input$domestication_result_msd_sequence, {
    generic_sequence_input("msd")
    updateTextAreaInput(session, "msd_input_sequence", value = rv$domestication_primers@newsequence)
    updateTabsetPanel(session, "MainNav", selected = "Saturation-Mutagenesis")
    updateTabsetPanel(session, "msd", selected = "Configuration")
    rv[["msd_input_sequence"]]<-rv$domestication_primers@newsequence
  })
  ####################################
  ######Single Point Mutagenesis#######
  #####################################
  #########INPUT##########
  generic_sequence_input("spm")
  generic_process_input("spm")
  generic_process_fasta_input("spm")
  ########Configuration#######
  generic_mut_conf("spm")
  generic_template_selection("spm")
  generic_re_selection("spm")
  generic_re_selection_lvl0("spm")
  levelsettings("spm")
    observeEvent(input$spm_configuration_next, {
        updateTabsetPanel(session, "spm", "Preview and Selection")
        #output$spm_preview<-renderUI(print_sequence(sequence = input$sp_seq, mutations = rv$sp_mutations))+
      generic_complex_preview_logic("spm")
      generic_preview("spm")
    })
    ############PREVIEW AND SELECTION#############
    generic_complex_selection("spm")
    
    ###########RESULTS################
    generic_results("spm")
    generic_primer_output("spm")
    
   ####################################
   ######Saturation Mutagenesis#######
   #####################################
   #########INPUT##########
    generic_sequence_input("msd", default_value = "ATGTCTCAGGTTCAGAGTGGCATTTTGCCAGAACATTGCCGCGCGGCGATTTGGATCGAAGCCAACGTGAAAGGGGAAGTTGACGCCCTGCGTGCGGCCAGTAAAACATTTGCCGACAAACTGGCAACTTTTGAAGCGAAATTCCCGGACGCGCATCTTGGTGCGGTGGTTGCCTTTGGTAACAACACCTGGCGCGCTCTGAGCGGCGGCGTTGGGGCAGAAGAGCTGAAAGATTTTCCGGGCTACGGTAAAGGCCTTGCGCCGACGACCCAGTTCGATGTGTTGATCCACATTCTTTCTCTGCGTCACGACGTAAACTTCTCTGTCGCCCAGGCGGCGATGGAAGCCTTTGGTGACTGCATTGAAGTGAAAGAAGAGATCCACGGCTTCCGTTGGGTTGAAGAGCGTGACCTGAGCGGCTTTGTTGACGGTACGGAAAACCCGGCGGGTGAAGAGACGCGTCGCGAAGTGGCGGTTATCAAAGACGGCGTGGATGCGGGCGGCAGCTATGTGTTTGTCCAGCGTTGGGAACACAACCTGAAGCAGCTCAACCGGATGAGCGTTCACGATCAGGAGATGGTGATCGGGCGCACCAAAGAGGCCAACGAAGAGATCGACGGCGACGAACGTCCGGAAACCTCTCACCTCACCCGCGTTGATCTGAAAGAAGATGGCAAAGGGCTGAAGATTGTTCGCCAGAGCCTGCCGTACGGCACTGCCAGTGGCACTCACGGTCTGTACTTCTGCGCCTACTGCGCGCGTCTGCATAACATTGAGCAGCAACTGCTGAGCATGTTTGGCGATACCGATGGTAAGCGTGATGCGATGTTGCGTTTCACCAAACCGGTAACCGGCGGCTATTATTTCGCACCGTCGCTGGACAAGTTGATGGCGCTGTAA")
    generic_process_input("msd", default_value = "ATGTCTCAGGTTCAGAGTGGCATTTTGCCAGAACATTGCCGCGCGGCGATTTGGATCGAAGCCAACGTGAAAGGGGAAGTTGACGCCCTGCGTGCGGCCAGTAAAACATTTGCCGACAAACTGGCAACTTTTGAAGCGAAATTCCCGGACGCGCATCTTGGTGCGGTGGTTGCCTTTGGTAACAACACCTGGCGCGCTCTGAGCGGCGGCGTTGGGGCAGAAGAGCTGAAAGATTTTCCGGGCTACGGTAAAGGCCTTGCGCCGACGACCCAGTTCGATGTGTTGATCCACATTCTTTCTCTGCGTCACGACGTAAACTTCTCTGTCGCCCAGGCGGCGATGGAAGCCTTTGGTGACTGCATTGAAGTGAAAGAAGAGATCCACGGCTTCCGTTGGGTTGAAGAGCGTGACCTGAGCGGCTTTGTTGACGGTACGGAAAACCCGGCGGGTGAAGAGACGCGTCGCGAAGTGGCGGTTATCAAAGACGGCGTGGATGCGGGCGGCAGCTATGTGTTTGTCCAGCGTTGGGAACACAACCTGAAGCAGCTCAACCGGATGAGCGTTCACGATCAGGAGATGGTGATCGGGCGCACCAAAGAGGCCAACGAAGAGATCGACGGCGACGAACGTCCGGAAACCTCTCACCTCACCCGCGTTGATCTGAAAGAAGATGGCAAAGGGCTGAAGATTGTTCGCCAGAGCCTGCCGTACGGCACTGCCAGTGGCACTCACGGTCTGTACTTCTGCGCCTACTGCGCGCGTCTGCATAACATTGAGCAGCAACTGCTGAGCATGTTTGGCGATACCGATGGTAAGCGTGATGCGATGTTGCGTTTCACCAAACCGGTAACCGGCGGCTATTATTTCGCACCGTCGCTGGACAAGTTGATGGCGCTGTAA")
    generic_process_fasta_input("msd")
    ########Configuration#######
    generic_mut_conf("msd")
    generic_template_selection("msd")
    generic_re_selection("msd")
    generic_re_selection_lvl0("msd")
    levelsettings("msd")
    observeEvent(input$msd_configuration_next, {
      updateTabsetPanel(session, "msd", "Preview and Selection")
      #output$spm_preview<-renderUI(print_sequence(sequence = input$sp_seq, mutations = rv$sp_mutations))
      generic_complex_preview_logic("msd", spm=F)
    
      generic_preview("msd")
    })
    ############PREVIEW AND SELECTION#############
    generic_complex_selection("msd", spm=F)
    
    ###########RESULTS################
    generic_results("msd", spm=F)
    generic_primer_output("msd")
    
     ####################################
     ######Quick Quality Control #######
     #####################################
    generic_sequence_input("qqc", button=F,default_value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA") 
    observeEvent(input[[paste("qqc", "sequence_file", sep="_")]],{
      fasta<-read.fasta(file=input[[paste("qqc", "sequence_file", sep="_")]]$datapath, seqtype = "DNA", as.string = T, seqonly = T)
      updateTextAreaInput(session, paste("qqc", "input_sequence", sep="_"), value = fasta[[1]])
      rv[[paste(sep="_", "qqc", "input_sequence")]]<-fasta[[1]]
      sequence_check(rv[[paste("qqc", "input_sequence", sep="_")]])
    })
    observeEvent(input$qqc_input_sequence,{
      if(nchar(input$qqc_input_sequence)>2) {
        codonseq<-splitseq(s2c(input$qqc_input_sequence))
        #codonseq<-sequence_check(input$qqc_sequence_input)
        rv$qqc_mutations_df<-data.frame(Mutations=as.numeric(c()), Codon<-as.character(c()), AminoAcid<-as.character(c()), stringsAsFactors = FALSE)
        colnames(rv$qqc_mutations_df)<-c("Mutation Position", "Codon", "Amino Acid Residue")
        

          my.insert.callback.qqc <- function(data, row) {
            if((data[row,1]>length(codonseq))){
              stop(paste("Pleaser enter a value between ", "1 and ", length(codonseq), ".", sep=""))
            }
            if(data[row,1]<1){
              stop(paste("Pleaser enter a value between ", "1 and ", length(codonseq), ".", sep=""))
            }
            codon<-codonseq[data[row,1]]
            data[row,2]<-codon
            data[row,3]<-translate(s2c(codon))
            rv$qqc_mutations_df <- data
            return(data)
          }
          
          my.update.callback.qqc <- function(data, olddata, row) {
            if((data[row,1]>length(codonseq))){
              stop(paste("Pleaser enter a value between ", "1 and ", length(codonseq), ".", sep=""))
            }
            if(data[row,1]<1){
              stop(paste("Pleaser enter a value between ", "1 and ", length(codonseq), ".", sep=""))
            }
            data[row,2] <- codonseq[data[row,1]]
            data[row,3] <- translate(s2c(data[row,2]))
            rv$qqc_mutations_df <- data
            return(data)
          }
          my.delete.callback.qqc <- function(data, olddata, row) {
            data <- data[-c(row),]
            rv$qqc_mutations_df <- data
            return(data)
          }
          
          DTedit::dtedit(input, output,
                     name = 'qqc_mutation_table',
                     thedata = rv$qqc_mutations_df,
                     edit.cols = c("Mutation Position"),
                     input.types = c("Mutation Position" ="numericInput"),
                     title.add = "Add Mutation",
                     label.add = "Add Mutation",
                     callback.update = my.update.callback.qqc,
                     callback.insert = my.insert.callback.qqc,
                     callback.delete = my.delete.callback.qqc,
                     datatable.options = list(searching=FALSE, buttons=c("csv")),
                     defaultPageLength = 5)
      }
      else{
        output$qqc_mutation_table<-renderUI(p("Please enter a sequence first."))
      }
    })
    

     observeEvent(input$qqc_next,{
       rv$plots<-NULL
       rv$plotlist<-NULL
       rv$qqc_mutations<-as.vector(rv$qqc_mutations_df[,1], mode = "numeric")
       updateTabsetPanel(session, "qqc", "Results")
       base_distribution_shiny(input$qqc_input_sequence, input$qqc_file$datapath, replacements = rv$qqc_mutations)
        rv$plots<-renderUI({
          output_list<-lapply(rv[["plotlist"]], function(i){renderPlotly(rv[[i]])})
          output_list<-lapply(output_list, function(i){attr(i, 'outputArgs')<-list(height="700px"); return(i)})
          return(do.call(tagList, output_list))
          })
        output$plots<-rv$plots
     })
     

}
# Run the application 
shinyApp(ui=ui, server=server)

