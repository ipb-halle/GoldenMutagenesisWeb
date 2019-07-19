#INPUT SEQUENCE PROCESSING
generic_process_input<-function(prefix, next_panel="Configuration", default_value="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACGATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGTCGACAAGCTTGCGGCCGCACTCGAGTGA"){
  observeEvent(input[[paste(prefix, "sequence_next", sep="_")]], {
    updateTextAreaInput(session, paste(prefix, "input_sequence", sep="_"), value = str_remove_all(input[[paste(prefix, "input_sequence", sep="_")]], "\\s+"))
    if(input[[paste(prefix, "input_sequence", sep="_")]] == "") {
      shinyalert("No Sequence!", "You have not entered a sequence. The default value will be used!", type = "warning")
      updateTextAreaInput(session, paste(prefix, "input_sequence", sep="_"), value = default_value)
    }
    sequence_check(input[[paste(prefix, "input_sequence", sep="_")]])
    updateTabsetPanel(session, prefix, next_panel)
  })
}

#MUTAGENESIS TEMPLATE SELECTION
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
#RESTRICTION ENZYME SELECTION
generic_re_selection<-function(prefix) {observeEvent(input[[paste(prefix,"re_enzyme_selection",sep="_")]], {
  if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bbsi") {
    updateTextInput(session, paste(prefix,"re_enzyme",sep="_"), value = "GAAGAC")
    updateTextInput(session, paste(prefix, "suffix", sep="_"), value = "AA")
    updateTextInput(session, paste(prefix, "prefix", sep="_"), value = "TT")
    
  }
  if(input[[paste(prefix,"re_enzyme_selection",sep="_")]] == "bsai") {
    updateTextInput(session, paste(prefix,"re_enzyme",sep="_"), value = "GGTCTC")
    updateTextInput(session, paste(prefix, "suffix", sep="_"), value = "A")
    updateTextInput(session, paste(prefix, "prefix", sep="_"), value = "TT")
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
  observeEvent(input[[paste(prefix, "codonpos", sep="_")]], {
    output[[paste(prefix, "aa", sep="_")]]<-renderUI({
      HTML(aaa(translate(s2c(sequence_check(input[[paste(prefix, "input_sequence", sep="_")]])[as.numeric(input[[paste(prefix, "codonpos", sep="_")]])]))))}
    )
    if(length(rv[[paste(prefix, "mutations", sep="_")]])>0){
      positions_aa<-c()
      for(i in 1:length(rv[[paste(prefix, "mutations", sep="_")]])) {
        position_aa<-as.numeric(rv[[paste(prefix, "mutations", sep="_")]][[i]][1])
        positions_aa<-c(positions_aa, position_aa)
      }
      #print(positions_aa)
      if(as.numeric(input[[paste(prefix, "codonpos", sep="_")]]) %in% positions_aa) {
        #print("check")
        updateCheckboxInput(session, paste(prefix, "sm", sep="_"), value=T)
      }
      else {
        updateCheckboxInput(session, paste(prefix, "sm", sep="_"), value=F)
      }
    }
  })
  observeEvent(input[[paste(prefix, "sm_apply", sep="_")]], {
    codon<-as.numeric(input[[paste(prefix, "codonpos", sep="_")]])
    positions_aa<-c()
    selected<-input[[paste(prefix, "sm", sep="_")]]
    if(length(rv[[paste(prefix, "mutations", sep="_")]])>0){
      for(i in 1:length(rv[[paste(prefix, "mutations", sep="_")]])) {
        position_aa<-as.numeric(rv[[paste(prefix, "mutations", sep="_")]][[i]][1])
        positions_aa<-c(positions_aa, position_aa)
      }
      if(selected==F){
        ind<-which(positions_aa==codon)
        if(length(ind) > 0) { 
          rv[[paste(prefix, "mutations", sep="_")]]<-list.remove(rv[[paste(prefix, "mutations", sep="_")]], ind)
        }
      }
      if(selected==T){
        ind<-which(positions_aa==codon)
        if(length(ind) == 0) {
          rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, translate(s2c(input[[paste(prefix, "input_sequence", sep="_")]]))[codon]))
        }
      }
    } else {
      if(selected==T){
        rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, translate(s2c(input[[paste(prefix, "input_sequence", sep="_")]]))[codon]))
      }
    }
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = input[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
  })
}
generic_complex_preview_logic<-function(prefix, spm=T){
  observeEvent(input[[paste(prefix, "codonpos", sep="_")]], {
    output[[paste(prefix, "aa", sep="_")]]<-renderUI({
        HTML(aaa(translate(s2c(sequence_check(input[[paste(prefix, "input_sequence", sep="_")]])[as.numeric(input[[paste(prefix, "codonpos", sep="_")]])]))))}
    )
    if(length(rv[[paste(prefix, "mutations", sep="_")]])>0){
      positions_aa<-c()
      names_aa<-c()
      for(i in 1:length(rv[[paste(prefix, "mutations", sep="_")]])) {
        position_aa<-as.numeric(rv[[paste(prefix, "mutations", sep="_")]][[i]][1])
        name_aa<-rv[[paste(prefix, "mutations", sep="_")]][[i]][2]
        positions_aa<-c(positions_aa, position_aa)
        names_aa<-c(names_aa, name_aa)
      }
      #print(positions_aa)
      if(as.numeric(input[[paste(prefix, "codonpos", sep="_")]]) %in% positions_aa) {
        #print("check")
        ind<-which(positions_aa == as.numeric(input[[paste(prefix, "codonpos", sep="_")]]))
        if(spm==T){
          updateSelectInput(session, inputId = paste(prefix,"m", sep="_"), selected = aaa(names_aa[ind]))
        }
        else{
          updateSelectInput(session, inputId = paste(prefix,"m", sep="_"), selected = names_aa[ind])
        }
      }
      else{
        updateSelectInput(session, inputId = paste(prefix, "m", sep="_"), selected = "none")
      }
    }
  })
  
  observeEvent(input[[paste(prefix, "m_apply", sep="_")]], {
    #browser()
    codon<-as.numeric(input[[paste(prefix, "codonpos", sep = "_")]])
    positions_aa<-c()
    selected<-input[[paste(prefix, "m", sep="_")]]
    if(length(rv[[paste(prefix, "mutations", sep="_")]])>0){
      positions_aa<-c()
      names_aa<-c()
      for(i in 1:length(rv[[paste(prefix, "mutations", sep="_")]])) {
        position_aa<-as.numeric(rv[[paste(prefix, "mutations", sep="_")]][[i]][1])
        name_aa<-rv[[paste(prefix, "mutations", sep="_")]][[i]][2]
        positions_aa<-c(positions_aa, position_aa)
        names_aa<-c(names_aa, name_aa)
      }
    }
    if(selected=="none"){
      ind<-which(positions_aa==codon)
      if(length(ind) > 0) { 
        rv[[paste(prefix, "mutations", sep="_")]]<-list.remove(rv[[paste(prefix, "mutations", sep="_")]], ind)
      }
    }
    if(selected!="none"){
      ind<-which(positions_aa==codon)
      if(length(ind) == 0) {
        if(spm==T){
          rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, a(selected)))
        }
        else {
          rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, selected))
        }
      }
      else{
        rv[[paste(prefix, "mutations", sep="_")]]<-list.remove(rv[[paste(prefix, "mutations", sep="_")]], ind)
        if(spm==T){
          rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, a(selected)))
        } else{
          rv[[paste(prefix, "mutations", sep="_")]]<-list.append(rv[[paste(prefix, "mutations", sep="_")]], c(codon, selected))
        }
      }
    }
    output[[paste(prefix, "preview", sep="_")]]<-renderUI(print_sequence(sequence = input[[paste(prefix, "input_sequence", sep="_")]], mutations = rv[[paste(prefix, "mutations", sep="_")]]))
  })
  }
#Results
generic_results<-function(prefix, panel="Results", spm=T){
observeEvent(input[[paste(prefix, "selection_next", sep="_")]], {
    if(length(rv[[paste(prefix, "mutations", sep="_")]])==0) {
      shinyalert("No Mutations!", "No selected mutations, nothing to do!", type = "info")
    }
    else{
      updateTabsetPanel(session, prefix, panel)
      #print(input[[paste(prefix, "input_sequence", sep="_")]])
      if(spm==T){
        rv[[paste(prefix, "primers", sep="_")]]<-mutate_spm(input[[paste(prefix, "input_sequence", sep="_")]], prefix = input[[paste(prefix, "prefix", sep="_")]], restriction_enzyme = input[[paste(prefix, "re_enzyme", sep="_")]], suffix = input[[paste(prefix, "suffix", sep="_")]], vector = c(input[[paste(prefix, "v1", sep="_")]], input[[paste(prefix, "v2", sep="_")]]), replacements = rv[[paste(prefix, "mutations", sep="_")]],  binding_min_length = input[[paste(prefix, "binding_min_length", sep="_")]], target_temp = input[[paste(prefix, "temperature", sep="_")]], cuf = input[[paste(prefix, "cuf", sep="_")]], binding_max_length = input[[paste(prefix, "binding_max_length", sep="_")]], replacement_range = input[[paste(prefix, "replacement_range", sep="_")]],  fragment_min_size = input[[paste(prefix, "fragment_min_size", sep="_")]])
      } else {
        rv[[paste(prefix, "primers", sep="_")]]<-mutate_msd(input[[paste(prefix, "input_sequence", sep="_")]], prefix = input[[paste(prefix, "prefix", sep="_")]], restriction_enzyme = input[[paste(prefix, "re_enzyme", sep="_")]], suffix = input[[paste(prefix, "suffix", sep="_")]], vector = c(input[[paste(prefix, "v1", sep="_")]], input[[paste(prefix, "v2", sep="_")]]), replacements = rv[[paste(prefix, "mutations", sep="_")]],  binding_min_length = input[[paste(prefix, "binding_min_length", sep="_")]], target_temp = input[[paste(prefix, "temperature", sep="_")]], binding_max_length = input[[paste(prefix, "binding_max_length", sep="_")]], replacement_range = input[[paste(prefix, "replacement_range", sep="_")]], fragment_min_size = input[[paste(prefix, "fragment_min_size", sep="_")]])
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
      #lbls <- paste(lbls, pct) # add percents to labels
      #lbls <- paste(lbls,"%",sep="") # ad % to labels
      #print(element)
      #file<-tempfile(fileext = ".png")
      #plotlist<-c(plotlist,file)
      #png(file)
      #pie(slices,labels = lbls, col=brewer.pal(4,"Spectral"),main = paste("Peak intensity distribution for \nPosition", pattern_pos[element], "(Template) -", subject_pos[element], "(Sequencing)", sep=" "))
      #browser()
      rv[[paste0(pattern_pos[element],"_",subject_pos[element])]]<-layout(plot_ly(data.frame(cbind(slices, lbls)), labels=~lbls, values=~slices, type="pie", colors=brewer.pal(4, "Spectral"), height=700),title=list(text=paste("Peak intensity distribution for",aaa(translate(s2c(input_sequence))[ceiling(pattern_pos[element]/3)]) ,as.character(ceiling(pattern_pos[element]/3)), "\nPosition", pattern_pos[element], "(Template) -", subject_pos[element], "(Sequencing)", sep=" "), font=list(size=16)), margin=list(t=120, b=100))
      rv[["plotlist"]]<-c(rv[["plotlist"]], paste0(pattern_pos[element],"_",subject_pos[element]))
      #dev.off()
    }
    #print(plotlist)
  })
}
