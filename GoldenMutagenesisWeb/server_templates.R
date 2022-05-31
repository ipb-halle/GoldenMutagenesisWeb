#INPUT SEQUENCE PROCESSING
generic_process_input<-function(prefix, next_panel="Configuration", default_value="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA"){
  observeEvent(input[[paste(prefix, "sequence_next", sep="_")]], {
    updateTextAreaInput(session, paste(prefix, "input_sequence", sep="_"), value = str_remove_all(input[[paste(prefix, "input_sequence", sep="_")]], "\\s+"))
    rv[[paste(sep="_", prefix, "input_sequence")]]<-str_remove_all(input[[paste(prefix, "input_sequence", sep="_")]], "\\s+")
    if(input[[paste(prefix, "input_sequence", sep="_")]] == "") {
      shinyalert("No Sequence!", "You have not entered a sequence. The default value will be used!", type = "warning")
      updateTextAreaInput(session, paste(prefix, "input_sequence", sep="_"), value = default_value)
      rv[[paste(sep="_", prefix, "input_sequence")]]<-default_value
    } else{
    sequence_check(rv[[paste(prefix, "input_sequence", sep="_")]])
    #rv[[paste(sep="_", prefix, "input_sequence")]]<-input[[paste(sep="_", prefix, "input_sequence")]]
    }
    showTab(prefix, next_panel, select = T)
  })
}
generic_process_fasta_input<-function(prefix,next_panel="Configuration"){
  observeEvent(input[[paste(prefix, "sequence_fasta_next", sep="_")]],{
    fasta<-read.fasta(file=input[[paste(prefix, "sequence_file", sep="_")]]$datapath, seqtype = "DNA", as.string = T, seqonly = T)
    updateTextAreaInput(session, paste(prefix, "input_sequence", sep="_"), value = fasta[[1]])
    rv[[paste(sep="_", prefix, "input_sequence")]]<-fasta[[1]]
    sequence_check(rv[[paste(prefix, "input_sequence", sep="_")]])
    showTab(prefix, next_panel, select = T)
  })
}

#MUTAGENESIS TEMPLATE SELECTION
generic_template_selection<-function(prefix) {
  observeEvent(input[[paste(prefix, "template", sep="_")]], {
    if(input[[paste(prefix, "template", sep="_")]] == "1") {
      print("pAGM9121 selected")
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv0")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bbsi")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "CTCA")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "CTCG")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv")
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/51833/','_blank')"
    }
    if(input[[paste(prefix, "template", sep="_")]] == "2") {
      print("pAGM22082 selected")
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv") 
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/117225/','_blank')"
    }
    if(input[[paste(prefix, "template", sep="_")]] == "3") {
      print("pICH86988 selected")
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "e_coli_316407.csv")
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/48076/','_blank')"
    }
    if(input[[paste(prefix, "template", sep="_")]] == "4") {
      print("pPAP001 selected")
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "Pichia_pastoris_4922.csv") 
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/153488/','_blank')"
      rv[[paste(sep="_", prefix, "publication_button")]]<-"window.open('https://www.biorxiv.org/content/10.1101/2020.07.22.216432v1','_blank')"
      show(id = paste(sep="_", prefix, "publication_link"), anim = T)
    }
    if(input[[paste(prefix, "template", sep="_")]] == "5") {
      print("pPAP002 selected")
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "Pichia_pastoris_4922.csv") 
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/153489/','_blank')"
      rv[[paste(sep="_", prefix, "publication_button")]]<-"window.open('https://www.biorxiv.org/content/10.1101/2020.07.22.216432v1','_blank')"
      show(id = paste(sep="_", prefix, "publication_link"), anim = T)
    }
    if(input[[paste(prefix, "template", sep="_")]] == "6") {
      print("pAGT572_Nemo selected")
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv2")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bsai")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "s_cerevisiae_4932.csv") 
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/153487/','_blank')"
      rv[[paste(sep="_", prefix, "publication_button")]]<-"window.open('https://www.biorxiv.org/content/10.1101/2020.07.22.216432v1','_blank')"
      show(id = paste(sep="_", prefix, "publication_link"), anim = T)
    }
    if(input[[paste(prefix, "template", sep="_")]] == "7") {
      print("pICH41308_CDS1 selected")
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv0")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bbsi")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "AAGC")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "arabidopsis.csv")
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/47998/','_blank')"
    }
    if(input[[paste(prefix, "template", sep="_")]] == "8") {
      print("pAGM1287_CDS1ns selected")
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      disable(paste(sep="_", prefix, "level"))
      disable(paste(sep="_", prefix, "re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      disable(paste(sep="_", prefix, "lvl0_v1"))
      disable(paste(sep="_", prefix, "lvl0_v2"))
      disable(paste(sep="_", prefix, "v1"))
      disable(paste(sep="_", prefix, "v2"))
      disable(paste(sep="_", prefix, "av1"))
      disable(paste(sep="_", prefix, "av2"))
      updateSelectInput(session, paste(sep="_", prefix, "level"), selected = "lv0")
      updateSelectInput(session, paste(sep="_", prefix, "re_enzyme_selection"), selected = "bbsi")
      updateTextInput(session, paste(sep="_", prefix, "v1"), value = "AATG")
      updateTextInput(session, paste(sep="_", prefix, "v2"), value = "CGAA")
      updateSelectInput(session, paste(sep="_", prefix, "cuf"), selected =  "arabidopsis.csv")
      rv[[paste(sep="_", prefix, "addgene_button")]]<-"window.open('https://www.addgene.org/47996/','_blank')"
    }
    if(input[[paste(prefix, "template", sep="_")]] == "c") {
      hide(id = paste(sep="_", prefix, "link"), anim = T)
      hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
      print("custom selected")
      enable(paste(sep="_", prefix, "level"))
      enable(paste(sep="_", prefix, "re_enzyme_selection"))
      enable(paste(sep="_", prefix, "lvl0_re_enzyme_selection"))
      enable(paste(sep="_", prefix, "v1"))
      enable(paste(sep="_", prefix, "v2"))
      enable(paste(sep="_", prefix, "av1"))
      enable(paste(sep="_", prefix, "av2"))
      enable(paste(sep="_", prefix, "lvl0_v1"))
      enable(paste(sep="_", prefix, "lvl0_v2"))
    }
    else{
      show(id = paste(sep="_", prefix, "link"), anim = T)
      #hide(id = paste(sep="_", prefix, "publication_link"), anim = T)
    }
  })
      onclick(paste(sep="_", prefix, "link"), runjs(rv[[paste(sep="_", prefix, "addgene_button")]]))
      onclick(paste(sep="_", prefix, "publication_link"), runjs(rv[[paste(sep="_", prefix, "publication_button")]]))
}
#RESTRICTION ENZYME SELECTION
generic_re_selection<-function(prefix) {observeEvent(input[[paste(prefix,"re_enzyme_selection",sep="_")]], {
  if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bbsi") {
    updateTextInput(session, paste(prefix,"re_enzyme",sep="_"), value = "GAAGAC")
    updateTextInput(session, paste(prefix, "suffix", sep="_"), value = "AA")
    updateTextInput(session, paste(prefix, "prefix", sep="_"), value = "TT")
    disable(paste(sep="_", prefix, "re_enzyme"))
    disable(paste(sep="_", prefix, "suffix"))
    disable(paste(sep="_", prefix, "prefix"))

  }
  if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bsai") {
    updateTextInput(session, paste(prefix,"re_enzyme",sep="_"), value = "GGTCTC")
    updateTextInput(session, paste(prefix, "suffix", sep="_"), value = "A")
    updateTextInput(session, paste(prefix, "prefix", sep="_"), value = "TT")
    disable(paste(sep="_", prefix, "re_enzyme"))
    disable(paste(sep="_", prefix, "suffix"))
    disable(paste(sep="_", prefix, "prefix"))
  }
  if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "c"){
    enable(paste(sep="_", prefix, "re_enzyme"))
    enable(paste(sep="_", prefix, "suffix"))
    enable(paste(sep="_", prefix, "prefix"))
  }
})
}
generic_re_selection_lvl0<-function(prefix) {observeEvent(input[[paste(prefix,"lvl0_re_enzyme_selection",sep="_")]], {
  if(input[[paste(prefix, "level", sep="_")]] == "lv2") {
    generic_re_selection(paste(prefix,"lvl0",sep="_"))
  }
})
}
#Preview
generic_simple_preview_logic<-function(prefix){
  codonseq<-splitseq(s2c(rv[[paste(prefix, "input_sequence", sep="_")]]))
  if(length(rv[[paste(sep="_", prefix, "mutations")]])>0) {
    pos<-sapply(rv[[paste(prefix, "mutations", sep="_")]], function(x) as.numeric(x[1]))
    codon<-sapply(rv[[paste(prefix, "mutations", sep="_")]], function(x) as.character(x[2]))
    rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=pos, Codon=codonseq[pos], AminoAcid=codon, stringsAsFactors = FALSE)       
  } else {
    rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=as.numeric(c()), Codon=as.character(c()), AminoAcid=as.character(c()), stringsAsFactors = FALSE)
  }
  colnames(rv[[paste(sep="_", prefix, "mutations_df")]])<-c("Mutation Position", "Codon", "Amino Acid Residue")
  my.insert.callback.local <- function(data, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    codon<-codonseq[data[row,1]]
    data[row,2]<-as.character(codon)
    data[row,3]<-translate(s2c(codon))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], x[3])))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  
  my.update.callback.local <- function(data, olddata, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    data[row,2] <- codonseq[data[row,1]]
    data[row,3] <- translate(s2c(data[row,2]))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], x[3])))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  my.delete.callback.local <- function(data, olddata, row) {
    data <- data[-c(row),]
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], x[3])))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  DTedit::dtedit(input, output,
                 name = paste(prefix, "mutation_table", sep="_"),
                 thedata = rv[[paste(sep="_", prefix, "mutations_df")]],
                 edit.cols = c("Mutation Position"),
                 input.types = c("Mutation Position" ="numericInput"),
                 title.add = "Add Mutation",
                 label.add = "Add Mutation",
                 callback.update = my.update.callback.local,
                 callback.insert = my.insert.callback.local,
                 callback.delete = my.delete.callback.local,
                 datatable.options = list(searching=FALSE, buttons=c("csv")),
                 show.copy = F,
                 defaultPageLength = 5)
}

generic_complex_preview_logic<-function(prefix, spm=T){
  codonseq<-splitseq(s2c(rv[[paste(prefix, "input_sequence", sep="_")]]))
  if((length(rv[[paste(sep="_", prefix, "mutations")]])>0) & is.null(rv[[paste(prefix, "mutations_df", sep="_")]])) {
    pos<-sapply(rv[[paste(prefix, "mutations", sep="_")]], function(x) as.numeric(x[1]))
    codon<-sapply(rv[[paste(prefix, "mutations", sep="_")]], function(x) as.character(x[2]))
    if(spm==T){
      rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=pos, Codon=codonseq[pos], AminoAcid=aaa(codon), Replacement=factor(aaa(codon), levels=aaa()), stringsAsFactors = FALSE) 
    } else{
      rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=pos, Codon=codonseq[pos], AminoAcid=aaa(codon), Replacement=factor(aaa(codon), levels=c("NNN", "NNK", "NNS", "NDT", "DBK", "NRT", "VHG", "VRK", "NYC", "KST", "TGG")), stringsAsFactors = FALSE) 
    }
  } else if(is.null(rv[[paste(prefix, "mutations_df", sep="_")]])){
    if(spm==T){
      rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=as.numeric(c()), Codon=as.character(c()), AminoAcid=as.character(c()),Replacement=factor(as.character(c()), levels=aaa()), stringsAsFactors = FALSE)
    } else {
      rv[[paste(sep="_", prefix, "mutations_df")]]<-data.frame(Mutations=as.numeric(c()), Codon=as.character(c()), AminoAcid=as.character(c()),Replacement=factor(as.character(c()), levels=c("NNN", "NNK", "NNS", "NDT", "DBK", "NRT", "VHG", "VRK", "NYC", "KST", "TGG")), stringsAsFactors = FALSE)
    }
  }
  colnames(rv[[paste(sep="_", prefix, "mutations_df")]])<-c("Mutation Position", "Codon", "Amino Acid Residue", "Mutation")
  my.insert.callback.local <- function(data, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    codon<-codonseq[data[row,1]]
    data[row,2]<-as.character(codon)
    data[row,3]<-aaa(translate(s2c(codon)))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], a(as.character(x[,4])))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  
  my.update.callback.local <- function(data, olddata, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    data[row,2] <- codonseq[data[row,1]]
    data[row,3] <- aaa(translate(s2c(data[row,2])))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], a(as.character(x[,4])))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  my.delete.callback.local <- function(data, olddata, row) {
    data <- data[-c(row),]
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], a(as.character(x[,4])))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  my.insert.callback.msd <- function(data, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    codon<-codonseq[data[row,1]]
    data[row,2]<-as.character(codon)
    data[row,3]<-aaa(translate(s2c(codon)))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], as.character(x[,4]))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  
  my.update.callback.msd <- function(data, olddata, row) {
    if((data[row,1]>length(codonseq)-1)){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    if(data[row,1]<2){
      stop(paste("Please enter a value between ", "2 and ", length(codonseq)-1, ".", sep=""))
    }
    data[row,2] <- codonseq[data[row,1]]
    data[row,3] <- aaa(translate(s2c(data[row,2])))
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], as.character(x[,4]))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  my.delete.callback.msd <- function(data, olddata, row) {
    data <- data[-c(row),]
    rv[[paste(sep="_", prefix, "mutations_df")]] <- data
    rv[[paste(prefix, "mutations", sep="_")]]<-alply(rv[[paste(sep="_", prefix, "mutations_df")]], 1, function(x) as.character(c(x[1], as.character(x[,4]))))
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = rv[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
    return(data)
  }
  if(spm==T){
    DTedit::dtedit(input, output,
                     name = paste(prefix, "mutation_table", sep="_"),
                     thedata = rv[[paste(sep="_", prefix, "mutations_df")]],
                     edit.cols = c("Mutation Position", "Mutation"),
                     input.types = c("Mutation Position" ="numericInput", "Mutation"="selectInput"),
                     title.add = "Add Mutation",
                     label.add = "Add Mutation",
                     callback.update = my.update.callback.local,
                     callback.insert = my.insert.callback.local,
                     callback.delete = my.delete.callback.local,
                     datatable.options = list(searching=FALSE, buttons=c("csv")),
                     show.copy = F,
                     defaultPageLength = 5)
  } else {
    DTedit::dtedit(input, output,
                   name = paste(prefix, "mutation_table", sep="_"),
                   thedata = rv[[paste(sep="_", prefix, "mutations_df")]],
                   edit.cols = c("Mutation Position", "Mutation"),
                   input.types = c("Mutation Position" ="numericInput", "Mutation"="selectInput"),
                   title.add = "Add Mutation",
                   label.add = "Add Mutation",
                   callback.update = my.update.callback.msd,
                   callback.insert = my.insert.callback.msd,
                   callback.delete = my.delete.callback.msd,
                   datatable.options = list(searching=FALSE, buttons=c("csv")),
                   show.copy = F,
                   defaultPageLength = 5)
  }
}

#Results
generic_results<-function(prefix, panel="Results", spm=T){
observeEvent(input[[paste(prefix, "selection_next", sep="_")]], {
#browser()
     if(length(rv[[paste(prefix, "mutations", sep="_")]])==0) {
      shinyalert("No Mutations!", "No selected mutations, nothing to do!", type = "info")
    }
    else{
      showTab(prefix, panel, select = T)
      print(input[[paste(prefix, "input_sequence", sep="_")]])
      print(rv[[paste(prefix, "mutations", sep="_")]])
      if(spm==T){
        rv[[paste(prefix, "primers", sep="_")]]<-mutate_spm(rv[[paste(prefix, "input_sequence", sep="_")]], prefix = input[[paste(prefix, "prefix", sep="_")]], restriction_enzyme = input[[paste(prefix, "re_enzyme", sep="_")]], suffix = input[[paste(prefix, "suffix", sep="_")]], vector = c(input[[paste(prefix, "v1", sep="_")]], input[[paste(prefix, "v2", sep="_")]]), replacements = rv[[paste(prefix, "mutations", sep="_")]],  binding_min_length = input[[paste(prefix, "binding_min_length", sep="_")]], target_temp = input[[paste(prefix, "temperature", sep="_")]], cuf = input[[paste(prefix, "cuf", sep="_")]], binding_max_length = input[[paste(prefix, "binding_max_length", sep="_")]], replacement_range = input[[paste(prefix, "replacement_range", sep="_")]],  fragment_min_size = input[[paste(prefix, "fragment_min_size", sep="_")]])
      } else {
        rv[[paste(prefix, "primers", sep="_")]]<-mutate_msd(rv[[paste(prefix, "input_sequence", sep="_")]], prefix = input[[paste(prefix, "prefix", sep="_")]], restriction_enzyme = input[[paste(prefix, "re_enzyme", sep="_")]], suffix = input[[paste(prefix, "suffix", sep="_")]], vector = c(input[[paste(prefix, "v1", sep="_")]], input[[paste(prefix, "v2", sep="_")]]), replacements = rv[[paste(prefix, "mutations", sep="_")]],  binding_min_length = input[[paste(prefix, "binding_min_length", sep="_")]], target_temp = input[[paste(prefix, "temperature", sep="_")]], binding_max_length = input[[paste(prefix, "binding_max_length", sep="_")]], replacement_range = input[[paste(prefix, "replacement_range", sep="_")]], fragment_min_size = input[[paste(prefix, "fragment_min_size", sep="_")]])
      }
      if(input[[paste(prefix, "level", sep="_")]]=="lv0") {
        if(input[[paste(prefix, "prepare_lvl2", sep="_")]]==TRUE) {
          rv[[paste(prefix, "primers", sep="_")]]<-primer_prepare_level(rv[[paste(prefix, "primers", sep="_")]], vector = c(input[[paste(prefix, "av1", sep="_")]], input[[paste(prefix, "av2", sep="_")]]))
        }
      }
      if(input[[paste(prefix, "level", sep="_")]]=="lv2") {
        if(input[[paste(prefix, "add_lvl0", sep="_")]] == TRUE){
          rv[[paste(prefix, "primers", sep="_")]]<-primer_add_level(primerset = rv[[paste(prefix, "primers", sep="_")]], prefix = input[[paste(prefix, "lvl0_prefix", sep="_")]], suffix =  input[[paste(prefix, "lvl0_suffix", sep="_")]], restriction_enzyme = input[[paste(prefix, "lvl0_re_enzyme", sep="_")]], vector = c(input[[paste(prefix, "lvl0_v1", sep="_")]], input[[paste(prefix, "lvl0_v2", sep="_")]]))
        }
      }
      output[[paste(prefix, "primers", sep="_")]]<-renderUI(print_primer_fancy(rv[[paste(prefix, "primers", sep="_")]]))
      output[[paste(sep="_", prefix, "dl_report_pdf")]]<-download_report(rv[[paste(prefix, "primers", sep="_")]])
      output[[paste(sep="_", prefix, "dl_protocol_pdf")]]<-downloadHandler(filename = "protocol.pdf", content = function(file){file.copy("www/files/protocol.pdf", file)})
      output[[paste(sep="_", prefix, "dl_report_txt")]]<-download_report_txt(rv[[paste(prefix, "primers", sep="_")]])
      }
  })
}

#QQC

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
    shinyalert(title="The alignment is in reverse direction.", type = "info", showCancelButton = F, showConfirmButton = T)
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
  if(length(sums_row)==0){
    shinyalert(paste("Could not find mutations for any positions.", "The tracematrix was empty because of the internal threshold. This means the basecall quality in the sequencing file ist too low.", sep="\n"), type = "warning")
  }
  tracematrix_subject<-as.data.frame(tracematrix_subject[sums_row,])
  if(length(setdiff(replacements, unique(ceiling(pattern_pos/3))))){
    shinyalert(paste("Could not find mutations for position(s):",paste(sapply(setdiff(replacements, unique(ceiling(pattern_pos/3))), function(i){paste(aaa(translate(s2c(input_sequence))[i]), i, sep=" ")}), collapse = "\n"), sep="\n"), type = "warning")
  }
  local({
    for(element in sums_row) {
      # plotting as pie chart
      sliceit <- dplyr::slice (tracematrix_subject,element)
      slices <- as.numeric(sliceit)
      lbls <- c("Adenine", "Cytosine", "Guanine", "Thymine")
      if(reverse==T) {
        lbls <- c("Thymine", "Guanine", "Cytosine", "Adenine")
      }
      pct <- round(slices/sum(slices)*100)
      rv[[paste0(pattern_pos[element],"_",subject_pos[element])]]<-layout(plot_ly(data.frame(cbind(slices, lbls)), labels=~lbls, values=~slices, type="pie", marker = list(colors=brewer.pal(4, "Spectral")), sort=F, height=700),title=list(text=paste("Peak intensity distribution for",aaa(translate(s2c(input_sequence))[ceiling(pattern_pos[element]/3)]) ,as.character(ceiling(pattern_pos[element]/3)), "\nPosition", pattern_pos[element], "(Template) -", subject_pos[element], "(Sequencing)", sep=" "), font=list(size=16)), margin=list(t=120, b=100))
      rv[["plotlist"]]<-c(rv[["plotlist"]], paste0(pattern_pos[element],"_",subject_pos[element]))
    }
  })
}
