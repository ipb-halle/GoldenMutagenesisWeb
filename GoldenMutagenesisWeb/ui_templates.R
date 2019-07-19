####GENRIC SEQUENCE INPUT#####
generic_sequence_input<-function(prefix, button=T , default_value="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA") {
  output[[paste(prefix, "input_panel", sep="_")]]<-renderUI({tagList(
                            h4("Sequence Input"),
                            br(), 
                            p("You can paste a sequence into the textbox or upload your own fasta file."),
                            tabsetPanel(type="pills",
                                        tabPanel("Manual Input",
                                                 textAreaInput(paste(prefix, "input_sequence", sep="_"),label = "Paste in your sequence", cols=60, rows = 10, resize = "both",
                                                               placeholder = default_value)
                                        ),
                                        tabPanel("FASTA Upload",
                                                 p("If your fasta file has more than one sequence, only the first one will be used."),
                                                 fileInput(paste(prefix, "sequence_file", sep=""), h4("Select .FA/.FASTA"), multiple = FALSE, accept = c(".fa", ".fasta"), width = NULL)
                                        ),
                                        if(button==T){
                                          fluidRow(column(2,actionButton(inputId = paste(prefix, "sequence_next", sep="_"), 'Next')))
                                        }
                                        else{
                                          htmlOutput("")
                                        }
                            )
                  )})
}

####MUTAGENESIS TEMPLATE GENERATION####
#############CONFIG###############
generic_mut_conf<-function(prefix){output[[paste(prefix, "mut_conf", sep="_")]]<-renderUI({tagList(
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
                     numericInput(paste(prefix, "binding_max_length",sep="_"), 
                                  "Maximal binding sequence length (AA)", 
                                  value = 9),
                     numericInput(paste(prefix,"temperature",sep="_"), 
                                  "Target temperature in Celsius", 
                                  value = 60),
                     numericInput(paste(prefix, "replacement_range", sep="_"),
                                  "Maximal distance between [AA] mutation sites to be integrated into one reverse primer",
                                  value = 3),
                     numericInput(paste(prefix, "fragment_min_size", sep="_"),
                                  "Minimum fragment size in BP",
                                  value = 100),
                     selectInput(paste(prefix,"cuf",sep="_"), "Codon Usage Table", list_cu_table())
           )
    )),
  fluidRow(column(12,wellPanel(style = "background: #f7ff89", h2("Golden Gate Level Settings"), uiOutput(paste(prefix,"levelsettings", sep="_"))))
  ), fluidRow(column(2,actionButton(inputId = paste(prefix, "configuration_next", sep="_"), 'Next'), br()))
  )
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

####SELECTION################
#############################
###########SIMPLE############
generic_preview<-function(prefix){
  if(!(paste(prefix, "mutations", sep="_") %in% names(rv))){
    rv[[paste(prefix, "mutations", sep="_")]]<-list()
  }
  output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = input[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
}
generic_simple_selection<-function(prefix){
  output[[paste(prefix, "codonnum", sep="_")]]<-renderUI(
    selectInput(paste(prefix,"codonpos", sep="_"), "Aminoacid Position", choices = 1:length(translate(s2c(input[[paste(prefix, "input_sequence", sep="_")]]))))
  )
  output[[paste(prefix, "preview_complete", sep="_")]]<-renderUI({tagList(
         wellPanel(uiOutput(paste(prefix, "preview", sep="_"))), wellPanel(style="background: #b7fff0",
                                                   fluidRow(
                                                     column(2,uiOutput(paste(prefix,"codonnum",sep="_")))
                                                   ),
                                                   fluidRow(
                                                     column(2,p("Aminoacid:"),uiOutput(paste(prefix,"aa", sep="_"))),
                                                     column(2,checkboxInput(paste(prefix, "sm", sep="_"), label = "Silent Mutation?",value = F)),
                                                     column(2, actionButton(paste(prefix,"sm_apply",sep="_"), "Apply"))
                                                   )
         ),                         
         fluidRow(column(2,actionButton(inputId = paste(prefix, "selection_next", sep="_"), 'Next'), br()))
  )})}
###########COMPLEX############
generic_complex_selection<-function(prefix, spm=T){
  output[[paste(prefix, "codonnum", sep="_")]]<-renderUI(
    selectInput(paste(prefix,"codonpos", sep="_"), "Aminoacid Position", choices = 1:length(translate(s2c(input[[paste(prefix, "input_sequence", sep="_")]]))))
  )
  if(spm==T){
    output[[paste(prefix, "m_o", sep="_")]]<-renderUI(selectInput(paste(prefix, "m", sep="_"), "Aminoacid Replacement", choices = c("none", aaa())))
  } else {
    output[[paste(prefix, "m_o", sep="_")]]<-renderUI(selectInput(paste(prefix, "m", sep="_"), "Saturation ", choices = c("none", "NNN", "NNK", "NNS", "NDT", "DBK", "NRT")))
  }
  output[[paste(prefix, "preview_complete", sep="_")]]<-renderUI({tagList(
    wellPanel(uiOutput(paste(prefix, "preview", sep="_"))),wellPanel(style="background: #b7fff0",
              fluidRow(
                column(2,uiOutput(paste(prefix,"codonnum", sep="_"))),
                column(2,uiOutput(paste(prefix, "m_o", sep="_")))
              ),
              fluidRow(
                column(2,p("Aminoacid:"),uiOutput(paste(prefix, "aa", sep="_"))),
                column(2, actionButton(paste(prefix,"m_apply", sep="_"), "Apply"))
              )
    ),
    fluidRow(column(2,actionButton(inputId = paste(prefix,"selection_next",sep="_"), 'Next'), br()))
  )})
}


#####RESULTS###############
generic_primer_output<-function(prefix){
  output[[paste(prefix, "primer_complete", sep="_")]]<-renderUI({tagList(
  fluidRow(column(12, h2("Legend")), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", "Prefix", "</span>",
                                                           "<span style=\"background-color: #fc9191; font-size: large;\">", "Restriction Enzyme", "</span>",
                                                           "<span style=\"background-color: #e2d544; font-size: large;\">", "Suffix", "</span>",
                                                           "<span style=\"background-color: #9fff8e; font-size: large;\">", "Vector", "</span>",
                                                           "<span style=\"background-color: #76fcb7; font-size: large;\">", "Overhang", "</span>",
                                                           "<span style=\"background-color: #f9ffd1; font-size: large;\">", "Extra", "</span>",
                                                           "<span style=\"background-color: #a5f1ff; font-size: large;\">", "Binding Sequence", "</span>",
                                                           sep="")))),br(),uiOutput(paste(prefix,"primers", sep="_"))
  )})}

