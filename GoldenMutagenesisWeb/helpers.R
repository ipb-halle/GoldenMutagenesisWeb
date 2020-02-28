mydlB<-function (outputId, label = "Download", icon="download", lib = "font-awesome" , class = NULL, ...) 
{
  aTag <- tags$a(id = outputId, class = paste("btn btn-default shiny-download-link", 
                                              class), href = "", target = "_blank", download = NA, 
                 icon(icon, lib=lib), label, ...)
}

sequence_check<-function(input_sequence){
  input_sequence<-str_to_upper(input_sequence)
  #print(input_sequence)
  if(nchar(input_sequence)%%3!=0) {
    shinyalert("Error", paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "), "error", closeOnEsc = F, showConfirmButton = F)
    stop(paste("The length of the sequence is no factor of 3. Please check your sequence.", "The length of the sequence was:", nchar(input_sequence),  sep=" "))
  }
  codon_seq<-splitseq(s2c(input_sequence))
  met<-which(str_detect(codon_seq, "ATG"))
  if(length(met) == 0) {
    shinyalert("Error", "No Methionine in the provided sequence. Stopping here. Please check the provided sequence.", "error", closeOnEsc = F, showConfirmButton = F)
    stop("No Methionine in the provided sequence. Stopping here. Please check the provided sequence.")
  }
  
  if(min(met) != 1){
    shinyalert("Warning", paste("No Methionine at first codon found! Please check the provided sequence! Took codon #", min(met), "as start.", sep=" "), closeOnEsc = F, showConfirmButton = T)
    warning(paste("No Methionine at first codon found! Please check the provided sequence! Took codon #", min(met), "as start.", sep=" "))
    codon_seq<-codon_seq[min(met):length(codon_seq)]
  } #else(codon_seq<-codon_seq[-1])
  
  stop<-which(str_detect(codon_seq, "(TAA)|(TGA)|(TAG)"))
  if(length(stop) == 0) {
    shinyalert("Error", "No stop codon in the provided sequence. Stopping here. Please check the provided sequence!", "error", closeOnEsc = F, showConfirmButton = F)
    stop("No stop codon in the provided sequence. Stopping here. Please check the provided sequence!")
  }
  
  if(max(stop) != length(codon_seq)) {
    shinyalert("Warning", paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "), closeOnEsc = F, showConfirmButton = T)
    warning(paste("There is no stop codon at the end of the sequence. Please check the provided sequence! Took codon #", max(stop), "as end.", sep= " "))  
    codon_seq<-codon_seq[1:max(stop)]
  }# else {
  #codon_seq <- codon_seq[-length(codon_seq)]
  #}
  return(codon_seq)
}

print_sequence<-function(sequence, mutations) {
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
    if((i %% 25 == 0) & (i != length(aa_seq))){
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
                     fluidRow(h5("Primer Forward"), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@prefix, "</span>",
                                                                          "<span style=\"background-color: #fc9191; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@restriction_enzyme, "</span>",
                                                                          "<span style=\"background-color: #e2d544; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@suffix, "</span>",
                                                                          "<span style=\"background-color: #9fff8e; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@vector, "</span>",
                                                                          "<span style=\"background-color: #76fcb7; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@overhang, "</span>",
                                                                          "<span style=\"background-color: #f9ffd1; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@extra, "</span>",
                                                                          "<span style=\"background-color: #a5f1ff; font-size: large; word-break: break-all\">", primerset@primers[[i]][[1]]@binding_sequence, "</span>",
                                                                          sep="")))),
                     fluidRow(column(4, h5("Melting temperature binding site"), round(primerset@primers[[i]][[1]]@temperature, digits=3)), column(3, h5("Temperature difference to setting"), round(primerset@primers[[i]][[1]]@difference, digits=3))),
                     fluidRow(h5("Primer Reverse"), column(12, HTML(paste("<span style=\"background-color: #fcfc92; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@prefix, "</span>",
                                                                          "<span style=\"background-color: #fc9191; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@restriction_enzyme, "</span>",
                                                                          "<span style=\"background-color: #e2d544; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@suffix, "</span>",
                                                                          "<span style=\"background-color: #9fff8e; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@vector, "</span>",
                                                                          "<span style=\"background-color: #76fcb7; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@overhang, "</span>",
                                                                          "<span style=\"background-color: #f9ffd1; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@extra, "</span>",
                                                                          "<span style=\"background-color: #a5f1ff; font-size: large; word-break: break-all\">", primerset@primers[[i]][[2]]@binding_sequence, "</span>",
                                                                          sep="")))),
                     fluidRow(column(4, h5("Melting temperature binding site"), round(primerset@primers[[i]][[2]]@temperature, digits = 3)), column(3, h5("Temperature difference to forward primer"), round(primerset@primers[[i]][[2]]@difference, digits = 3)))
    )
    renderlist<-list.append(renderlist, panel)
  }
  return(renderlist)
}
download_report<-function(primerset) {
  return(
    downloadHandler(
      filename="report.pdf",
      content = function(file) {
        show_modal_spinner("spring", text="Your download is being prepared!")
        tempReport<-file.path(tempdir(), "primer.Rmd")
        tempReportsu<-file.path(tempdir(), "primer_subunit.Rmd")
        tempReportsty<-file.path(tempdir(), "analysis.sty")
        tempReportprimers<-file.path(tempdir(), "primers.Rdata")
        primers<-primerset
        save(primers, file=tempReportprimers)
        file.copy("primer.Rmd", tempReport, overwrite = TRUE)
        file.copy("primer_subunit.Rmd", tempReportsu, overwrite = TRUE)
        file.copy("analysis.sty", tempReportsty, overwrite = TRUE)
        file.copy("figures", tempdir(), recursive = T)
        params <- list(path = tempReportprimers)
        rmarkdown::render(tempReport, output_file = file, output_dir = tempdir(), params = params,knit_root_dir=tempdir(),envir = new.env(parent = globalenv()))
        #knitr::knitr(input=tempReport, output=file, envir=new.env(parent = globalenv()))
        remove_modal_spinner()
      }
    )
  )
}
download_report_txt<-function(primerset) {
  return(
    downloadHandler(
      filename="report.txt",
      content = function(file) {
        sink(file, append = F, split = F)
        print_primer(primerset)
        sink()
      }
    )
  )
}
