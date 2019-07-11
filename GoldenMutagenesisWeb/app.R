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


## Read IPB coorporate identity
ipbheader <- HTML(readLines("ipbheader.html"))
ipbfooter <- HTML(readLines("ipbfooter.html"))

ui <- fluidPage(
  useShinyalert(),
  # IPB header 
  tags$head(tags$link(rel = "stylesheet", type = "text/css", 
                      href = "css/ipb-styles.css")),
  # Application title
  #titlePanel("GoldenMutagenesis Webtool Beta"),
  
  # Application title
  titlePanel(ipbheader, windowTitle = "GoldenMutagenesis"),
    navlistPanel(id="MainNav", widths = c(4, 7),
     tabPanel("Welcome", h4("Welcome to the GoldenMutagenesis Webtool."), h4("Please select the desired task.")), "Pre- and Postprocessing",
        tabPanel("Domestication", 
          mainPanel(width = 15,
            tabsetPanel(id="domestication",type="tabs",
                        tabPanel("Sequence Input", 
                        uiOutput("domestication_input_panel")),
                        tabPanel("Configuration", 
                                 fluidRow(h1("Domestication Configuration"), 
                                          wellPanel(style = "background: #ddffdd", p("Please select the restriction enzymes you want to domesticate."),
                                                    fluidRow(column(6,
                                                                    fluidRow(column(3,checkboxInput("d_bsai", "BsaI",value = T)), column(3,checkboxInput("d_bbsi", "BbsI", value = T))),br(),
                                                                    fluidRow(column(10,textInput("d_cu", "Custom Recognition Site"))),br(),
                                                                    fluidRow(column(2,actionButton("d_add_cu", "Add")), column(2, actionButton("d_rm_cu", "Remove")))),
                                                             column(6,h4("Restriction Sites"), uiOutput("d_cu_list"))
                                                    )),
                                 h1("Mutagenesis Configuration"),uiOutput("domestication_mut_conf")
                   ) 
                ),
                tabPanel("Preview and Selection",
                         uiOutput("domestication_preview_complete")),
                tabPanel("Results", 
                         uiOutput("domestication_primer_complete"),br(),fluidRow(column(5, actionButton("domestication_result_spm", "PointMutagenesis: Keep Mutations and Sequence"))),fluidRow(column(3, actionButton("domestication_result_spm_sequence", "PointMutagenesis: Keep Sequence"))), fluidRow(column(3, actionButton("domestication_result_msd_sequence", "SaturationMutagenesis: Keep Sequence"))))
                )
          )),
     tabPanel("Quick-Quality-Control", 
              tabsetPanel(id="qqc", type="tabs", tabPanel("Configuration",
                          wellPanel(style="background: #c9ffe6",  uiOutput("qqc_input_panel")),
                          wellPanel(style="background: #75ffbf", fileInput("qqc_file", h4("Select .ab1/.abf"), multiple = FALSE, accept = c(".ab1", ".abf"), width = NULL)),
                          wellPanel(style="background: #1df28f", h4("Mutation Positions"), fluidRow(column(3,numericInput("qqc_pos", "Sequence Position", value=1)), column(1,actionButton("qqc_add","Add")),column(1,actionButton("qqc_remove","Remove"))),
                                    fluidRow(uiOutput("qqc_mut_display"))),
                          fluidRow(actionButton("qqc_next", "Next"))
              ), tabPanel("Results", uiOutput("plots"))
              )
            ),
     "Mutagenesis",
     tabPanel("Point-Mutagenesis", 
              mainPanel(width = 15,
                        tabsetPanel(id="spm",type="tabs",
                                    tabPanel("Sequence Input",
                                      uiOutput("spm_input_panel")
                                    ),                 
                                    tabPanel("Configuration", h1("Mutagenesis Configuration") , uiOutput("spm_mut_conf")
                                      ),
                                    tabPanel("Preview and Selection",
                                             wellPanel(uiOutput("spm_preview_complete"))
                                    ),
                                    tabPanel("Results",
                                             uiOutput("spm_primer_complete")
                                    )
                        )
              )
      ),
     tabPanel("Saturation-Mutagenesis", mainPanel(width = 15,
                                                  tabsetPanel(id="msd",type="tabs",
                                                              tabPanel("Sequence Input",
                                                                       uiOutput("msd_input_panel")
                                                              ),
                                                              tabPanel("Configuration", h1("Mutagenesis Configuration") , uiOutput("msd_mut_conf")
                                                                       ),
                                                              tabPanel("Preview and Selection",
                                                                       wellPanel(uiOutput("msd_preview_complete"))
                                                                       ),
                                                              tabPanel("Results",
                                                                       uiOutput("msd_primer_complete")
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
  #########APPLY CONFIG TEMPLATES############
  source("server_templates.R", local=T)$value
  source("ui_templates.R", local = T)$value
 
 #####################################
 ###########DOMESTICATION#############
 #####################################
 #########INPUT##########
  generic_sequence_input("domestication", default_value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA")
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
        mutations<-domesticate(input$domestication_input_sequence, rs, cuf = input$domestication_cuf)
        if(length(mutations)>0){
          rv$domestication_mutations<-c(rv$domestication_mutations, mutations)
        }
      }
      generic_preview("domestication")  
    }
  })
  ############PREVIEW AND SELECTION#############
  generic_simple_selection("domestication")
  generic_simple_preview_logic("domestication")
  
  ###########RESULTS################
  generic_results("domestication")
  generic_primer_output("domestication")
  observeEvent(input$domestication_result_spm, {
    updateTextInput(session, "spm_input_sequence", value = rv$domestication_primers@oldsequence)
    rv$spm_mutations<-rv$domestication_mutations
    updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
    updateTabsetPanel(session, "spm", selected = "Configuration")
  })
  observeEvent(input$domestication_result_spm_sequence, {
    updateTextInput(session, "spm_input_sequence", value = rv$domestication_primers@newsequence)
    updateTabsetPanel(session, "MainNav", selected = "Point-Mutagenesis")
    updateTabsetPanel(session, "spm", selected = "Configuration")
  })
  observeEvent(input$domestication_result_msd_sequence, {
    updateTextInput(session, "msd_input_sequence", value = rv$domestication_primers@newsequence)
    updateTabsetPanel(session, "MainNav", selected = "Saturation-Mutagenesis")
    updateTabsetPanel(session, "msd", selected = "Configuration")
  })
  ####################################
  ######Single Point Mutagenesis#######
  #####################################
  #########INPUT##########
  generic_sequence_input("spm")
  ########Configuration#######
  generic_mut_conf("spm")
  generic_template_selection("spm")
  generic_re_selection("spm")
  generic_re_selection_lvl0("spm")
  levelsettings("spm")
    observeEvent(input$spm_configuration_next, {
        updateTabsetPanel(session, "spm", "Preview and Selection")
        #output$spm_preview<-renderUI(print_sequence(sequence = input$sp_seq, mutations = rv$sp_mutations))
      generic_preview("spm")
    })
    ############PREVIEW AND SELECTION#############
    generic_complex_selection("spm")
    generic_complex_preview_logic("spm")
    
    ###########RESULTS################
    generic_results("spm")
    generic_primer_output("spm")
    
   ####################################
   ######Saturation Mutagenesis#######
   #####################################
   #########INPUT##########
    generic_sequence_input("msd", default_value = "ATGTCTCAGGTTCAGAGTGGCATTTTGCCAGAACATTGCCGCGCGGCGATTTGGATCGAAGCCAACGTGAAAGGGGAAGTTGACGCCCTGCGTGCGGCCAGTAAAACATTTGCCGACAAACTGGCAACTTTTGAAGCGAAATTCCCGGACGCGCATCTTGGTGCGGTGGTTGCCTTTGGTAACAACACCTGGCGCGCTCTGAGCGGCGGCGTTGGGGCAGAAGAGCTGAAAGATTTTCCGGGCTACGGTAAAGGCCTTGCGCCGACGACCCAGTTCGATGTGTTGATCCACATTCTTTCTCTGCGTCACGACGTAAACTTCTCTGTCGCCCAGGCGGCGATGGAAGCCTTTGGTGACTGCATTGAAGTGAAAGAAGAGATCCACGGCTTCCGTTGGGTTGAAGAGCGTGACCTGAGCGGCTTTGTTGACGGTACGGAAAACCCGGCGGGTGAAGAGACGCGTCGCGAAGTGGCGGTTATCAAAGACGGCGTGGATGCGGGCGGCAGCTATGTGTTTGTCCAGCGTTGGGAACACAACCTGAAGCAGCTCAACCGGATGAGCGTTCACGATCAGGAGATGGTGATCGGGCGCACCAAAGAGGCCAACGAAGAGATCGACGGCGACGAACGTCCGGAAACCTCTCACCTCACCCGCGTTGATCTGAAAGAAGATGGCAAAGGGCTGAAGATTGTTCGCCAGAGCCTGCCGTACGGCACTGCCAGTGGCACTCACGGTCTGTACTTCTGCGCCTACTGCGCGCGTCTGCATAACATTGAGCAGCAACTGCTGAGCATGTTTGGCGATACCGATGGTAAGCGTGATGCGATGTTGCGTTTCACCAAACCGGTAACCGGCGGCTATTATTTCGCACCGTCGCTGGACAAGTTGATGGCGCTGTAA")
   ########Configuration#######
    generic_mut_conf("msd")
    generic_template_selection("msd")
    generic_re_selection("msd")
    generic_re_selection_lvl0("msd")
    levelsettings("msd")
    observeEvent(input$msd_configuration_next, {
      updateTabsetPanel(session, "msd", "Preview and Selection")
      #output$spm_preview<-renderUI(print_sequence(sequence = input$sp_seq, mutations = rv$sp_mutations))
      generic_preview("msd")
    })
    ############PREVIEW AND SELECTION#############
    generic_complex_selection("msd", spm=F)
    generic_complex_preview_logic("msd", spm=F)
    
    ###########RESULTS################
    generic_results("msd", spm=F)
    generic_primer_output("msd")
    
     ####################################
     ######Quick Quality Control #######
     #####################################
    generic_sequence_input("qqc", button=F,default_value = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA") 
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
       base_distribution_shiny(input$qqc_input_sequence, input$qqc_file$datapath, replacements = rv$qqc_mutations)
        #browser()
        rv$plots<-renderUI({
          output_list<-lapply(rv[["plotlist"]], function(i){renderPlotly(rv[[i]])})
          output_list<-lapply(output_list, function(i){attr(i, 'outputArgs')<-list(height="700px"); return(i)})
          return(do.call(tagList, output_list))
          })
        #browser()
        output$plots<-rv$plots
        
       # output$plotout <- renderUI({
       #   image_output_list <- 
       #     lapply(1:length(rv$plots),
       #            function(i)
       #            {
       #              imagename = paste0("image", i)
       #              imageOutput(imagename)
       #            })
       #   
       #   do.call(tagList, image_output_list)
       # })
       # observe({
       #   for (i in 1:length(rv$plots))
       #   {
       #     local({
       #       my_i <- i
       #       imagename = paste0("image", my_i)
       #       output[[imagename]] <- 
       #         renderImage({
       #           list(src = rv$plots[my_i],
       #                alt = "Image failed to render")
       #         }, deleteFile = FALSE)
       #     })
       #   }
       # })
     })
     

}
# Run the application 
shinyApp(ui=ui, server=server)

