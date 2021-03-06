####GENRIC SEQUENCE INPUT#####
generic_sequence_input<-function(prefix, button=T , default_value="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA") {
  output[[paste(prefix, "input_panel", sep="_")]]<-renderUI({tagList(
                            h4("Sequence Input"),
                            br(), 
                            p("You can paste a sequence into the textbox or upload your own fasta file."),
                            tabsetPanel(type="pills",
                                        tabPanel("Manual Input", wellPanel(
                                                 textAreaInput(paste(prefix, "input_sequence", sep="_"),label = "Paste in your sequence", cols=60, rows = 10, resize = "both",
                                                               placeholder = default_value),
                                          if(button==T){
                                            fluidRow(column(2,actionButton(inputId = paste(prefix, "sequence_next", sep="_"), 'Next')))
                                          }
                                          else{
                                             htmlOutput("")
                                          }
                                        )),
                                        tabPanel("FASTA Upload", wellPanel(
                                                 p("If your fasta file has more than one sequence, only the first one will be used."),
                                                 fluidRow(column(5,fileInput(paste(prefix, "sequence_file", sep="_"), h4("Select .FA/.FASTA"), multiple = FALSE, accept = c(".fa", ".fasta"), width = NULL))),
                                                 if(button==T){
                                                   fluidRow(column(2,actionButton(inputId = paste(prefix, "sequence_fasta_next", sep="_"), 'Next')))
                                                 }
                                                 else{
                                                   htmlOutput("")
                                                 }
                                        ))
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
                     p("You can select a preconfigured template to set those settings in accordiance to your Golden Gate Mutagenesis."),br(),
                     p("You can change the defaut values by selecting custom as pre-configuration."),
                     fluidRow(column(5,selectInput(paste(prefix,"template", sep="_"), "Pre-existing configuration:", c("pAGM9121"="1", 
                                                                                                     "pAGM22082 Red" = "2",
                                                                                                     "pICH86988" = "3",
                                                                                                     "pPAP001" = "4",
                                                                                                     "pPAP002" = "5",
                                                                                                     "pAGT572_Nemo" = "6",
                                                                                                     "custom" = "c"))),
                     column(4, style = "margin-top: 25px;",
                            actionButton(paste(prefix, "link", sep="_"), label = "View on Addgene", style=" background-image: url(./img/addgene.jpg); background-position: left; background-size: contain; background-repeat: no-repeat; padding-left: 40px;")),
                     column(2, style = "margin-top: 25px;",
                            actionButton(paste(prefix, "publication_link", sep="_"), label = "Publication", icon = icon("file-alt", lib = "font-awesome")))),
                     
                     selectInput(paste(prefix,"level",sep="_"), "Golden Gate Level:", c("Level0" = "lv0", "Level2" = "lv2")),
                     selectInput(paste(prefix,"re_enzyme_selection", sep="_"), "Restriction Enzyme:", c("BbsI"="bbsi",
                                                                                                        "BsaI"="bsai", 
                                                                                                        "custom" = "c")),
                     bsTooltip(id = paste(prefix,"re_enzyme",sep="_"), title = "Recognition site sequence of the respective restriction enzyme."),
                     textInput(paste(prefix,"re_enzyme", sep="_"), "Restriction Enzyme Recognition Site", value = "GAAGAC"),
                     bsTooltip(id = paste(prefix,"prefix",sep="_"), title = "Additional nucleobases in 5' position of the recognition site."),
                     textInput(paste(prefix, "prefix", sep="_"),  "Prefix", value = "TT"),
                     bsTooltip(id = paste(prefix,"suffix",sep="_"), title = "Spacer nucleotides matching the cleavage pattern of the enzyme."),
                     textInput(paste(prefix, "suffix",sep="_"), "Suffix", value = "AA"),
                     bsTooltip(id = paste(prefix,"v1",sep="_"), title = "Four basepair overhangs complementary to the created overhangs in the acceptor vector."),
                     bsTooltip(id = paste(prefix,"v2",sep="_"), title = "Four basepair overhangs complementary to the created overhangs in the acceptor vector."),
                     fluidRow(column(6, textInput(paste(prefix,"v1",sep="_"), "Vector Forward 5'", value = "CTCA")),column(6, textInput(paste(prefix,"v2",sep="_"), "Vector Reverse 3'", value = "CTCG")))
           )
    ),
    column(6,
           wellPanel(style = "background: #ffcccc",
                     h2("Algorithm Settings"),
                     p("Those settings are documentated in the GoldenMutagenesis R-Package. You do not need to change the default values in most applications."),
                     bsTooltip(id = paste(prefix,"binding_min_length",sep="_"), title = "The minimal threshold value of the length of the template binding sequence in amino acid residues"),
                     numericInput(paste(prefix,"binding_min_length",sep="_"), 
                                  "Minimal binding length (AA)", 
                                  value = 4),
                     bsTooltip(id = paste(prefix,"binding_max_length",sep="_"), title = "Maximal length of the binding sequence in amino acid residues"),
                     numericInput(paste(prefix, "binding_max_length",sep="_"), 
                                  "Maximal binding sequence length (AA)", 
                                  value = 9),
                     bsTooltip(id = paste(prefix,"temperature",sep="_"), title = "Melting temperature of the binding sequence in °C"),
                     numericInput(paste(prefix,"temperature",sep="_"), 
                                  "Target temperature in Celsius", 
                                  value = 60),
                     bsTooltip(id = paste(prefix,"replacement_range",sep="_"), title = "Maximum distance in amino acid residues between two randomization sites to be incoporated into a single primer (reverse, end of the fragment) - has a cascading effect for following mutations."),
                     numericInput(paste(prefix, "replacement_range", sep="_"),
                                  "Maximal distance between [AA] mutation sites to be integrated into one reverse primer",
                                  value = 3),
                     bsTooltip(id = paste(prefix,"fragment_min_size",sep="_"), title = "Minimal size of a generated gene fragment in base pairs."),
                     numericInput(paste(prefix, "fragment_min_size", sep="_"),
                                  "Minimum fragment size in BP",
                                  value = 100),
                     bsTooltip(id = paste(prefix,"cuf",sep="_"), title = "The Codon Usage Table which is being used to select the codon for an exchanged amino acid."),
                     selectInput(paste(prefix,"cuf",sep="_"), "Codon Usage Table", list_cu_table())
           )
    )),
  fluidRow(column(12,wellPanel(style = "background: #f7ff89", h2("Golden Gate Level Settings"), uiOutput(paste(prefix,"levelsettings", sep="_"))))
  ), fluidRow(column(2,actionButton(inputId = paste(prefix, "configuration_next", sep="_"), 'Next'), br()))
  )
})}

levelsettings<-function(prefix){
  level<-reactive({input[[paste(prefix, "level", sep="_")]]})
  template<-reactive({input[[paste(prefix, "template", sep="_")]]})
  output[[paste(prefix,"levelsettings", sep="_")]]<-renderUI({
    if(level() == "lv0") {
      if(template() == "c") {
      tagList (
        fluidRow(checkboxInput(paste(prefix,"prepare_lvl2", sep="_"), "Prepare for use in Level2", value = TRUE, width = NULL)),
        fluidRow(column(6, textInput(paste(prefix,"av1", sep="_"), "Accepting Vector Forward 5'", value = "AATG")),column(6, textInput(paste(prefix,"av2",sep="_"), "Accepting Vector Reverse 3'", value = "AAGC")))
      ) 
      } else {
        tagList (
          fluidRow(checkboxInput(paste(prefix,"prepare_lvl2", sep="_"), "Prepare for use in Level2", value = TRUE, width = NULL)),
          disabled(fluidRow(column(6, textInput(paste(prefix,"av1", sep="_"), "Accepting Vector Forward 5'", value = "AATG")),column(6, textInput(paste(prefix,"av2",sep="_"), "Accepting Vector Reverse 3'", value = "AAGC"))))
        ) 
       }
    }
    else{
        if(template() == "c") {
          tagList(
            fluidRow(column(5,checkboxInput(paste(prefix,"add_lvl0",sep="_"), "Use Level0 to go in Level2?", value = TRUE, width = NULL))),
            fluidRow(column(6,textInput(paste(prefix,"lvl0_v1",sep="_"), "Level0 Vector Forward 5'", value = "CTCA")),column(6, textInput(paste(prefix,"lvl0_v2",sep="_"), "Level0 Vector Reverse 3'", value = "CTCG"))),
            fluidRow(column(6,selectInput(paste(prefix,"lvl0_re_enzyme_selection",sep="_"), "Restriction Enzyme:", c("BbsI"="bbsi", "BsaI"="bsai", "custom" = "c"))),
                    disabled(column(6,textInput(paste(prefix,"lvl0_re_enzyme",sep="_"), "Restriction Enzyme Sequence", value = "GAAGAC")))),
            disabled(fluidRow(column(6,textInput(paste(prefix,"lvl0_prefix",sep="_"), "Prefix", value = "TT")),
                     column(6,textInput(paste(prefix,"lvl0_suffix",sep="_"), "Suffix", value = "AA"))))       
          ) 
        }
        else{
          tagList(
            fluidRow(column(5,checkboxInput(paste(prefix,"add_lvl0",sep="_"), "Use Level0 to go in Level2?", value = TRUE, width = NULL))),
            disabled(fluidRow(column(6,textInput(paste(prefix,"lvl0_v1",sep="_"), "Level0 Vector Forward 5'", value = "CTCA")),column(6, textInput(paste(prefix,"lvl0_v2",sep="_"), "Level0 Vector Reverse 3'", value = "CTCG")))),
            disabled(fluidRow(column(6,selectInput(paste(prefix,"lvl0_re_enzyme_selection",sep="_"), "Restriction Enzyme:", c("BbsI"="bbsi", "BsaI"="bsai", "custom" = "c"))),
                     column(6,textInput(paste(prefix,"lvl0_re_enzyme",sep="_"), "Restriction Enzyme Sequence", value = "GAAGAC")))),
            disabled(fluidRow(column(6,textInput(paste(prefix,"lvl0_prefix",sep="_"), "Prefix", value = "TT")),
                     column(6,textInput(paste(prefix,"lvl0_suffix",sep="_"), "Suffix", value = "AA"))))       
          ) 
        }
    }
  })  
}

####SELECTION################
#############################
###########SIMPLE############
generic_preview<-function(prefix){
  if(!(paste(prefix, "mutations", sep="_") %in% names(rv))){
    rv[[paste(prefix, "mutations", sep="_")]]<-list()
  }
  output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
}
generic_simple_selection<-function(prefix){
  output[[paste(prefix, "codonnum", sep="_")]]<-renderUI(
    selectInput(paste(prefix,"codonpos", sep="_"), "Aminoacid Position", choices = 1:length(translate(s2c(rv[[paste(prefix, "input_sequence", sep="_")]]))))
  )
  output[[paste(prefix, "preview_complete", sep="_")]]<-renderUI({tagList(
         wellPanel(uiOutput(paste(prefix, "preview", sep="_"))), wellPanel(style="background: #ffffff",
                                                   fluidRow(
                                                     column(5, uiOutput(paste(prefix, "mutation_table", sep="_")))
                                                   )

         ),                         
         fluidRow(column(2,actionButton(inputId = paste(prefix, "selection_next", sep="_"), 'Next'), br()))
  )})}
###########COMPLEX############
generic_complex_selection<-function(prefix, spm=T){
  output[[paste(prefix, "preview_complete", sep="_")]]<-renderUI({tagList(
    wellPanel(uiOutput(paste(prefix, "preview", sep="_"))),wellPanel(style="background: #ffffff; width=fit-content",
              fluidRow(
                column(8, p("Please note that a valid ORF is required. Thus, you can not edit the start and stop codons."))
              ),
              fluidRow(
                column(6,uiOutput(paste(prefix, "mutation_table", sep="_"))),
                bsTooltip(id = paste(prefix,"mutation_table",sep="_"), title = "Select the desired codons for saturation mutagenesis. TGG is available for using the 22c trick."),
              )),
    fluidRow(column(2,actionButton(inputId = paste(prefix,"selection_next",sep="_"), 'Next'), br()))
  )})
}


#####RESULTS###############
generic_primer_output<-function(prefix){
  output[[paste(prefix, "primer_complete", sep="_")]]<-renderUI({tagList(
  fluidRow(column(4, mydlB(paste(sep="_", prefix, "dl_report_pdf"), icon = "dna", lib="font-awesome", label = "Download Primer Report (PDF)")),
  column(4, mydlB(paste(sep="_", prefix, "dl_protocol_pdf"),icon="flask", lib="font-awesome", label = "Download Mutagenesis Protocol (PDF)"))),
  fluidRow(column(4, mydlB(paste(sep="_", prefix, "dl_report_txt"), icon="file-alt", lib="font-awesome", label = "Download Primer Report (TXT)"))),
  fluidRow(column(12, h2("Legend")), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large;\">", "Prefix", "</span>",
                                                           "<span style=\"background-color: #fc9191; font-size: large;\">", "Restriction Enzyme", "</span>",
                                                           "<span style=\"background-color: #e2d544; font-size: large;\">", "Suffix", "</span>",
                                                           "<span style=\"background-color: #9fff8e; font-size: large;\">", "Vector", "</span>",
                                                           "<span style=\"background-color: #76fcb7; font-size: large;\">", "Overhang", "</span>",
                                                           "<span style=\"background-color: #f9ffd1; font-size: large;\">", "Extra", "</span>",
                                                           "<span style=\"background-color: #a5f1ff; font-size: large;\">", "Binding Sequence", "</span>",
                                                           sep="")))),br(),uiOutput(paste(prefix,"primers", sep="_"))
  )})}

