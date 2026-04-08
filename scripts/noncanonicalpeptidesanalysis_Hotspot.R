# ABOUTME: Part of the DIA-NN Level 1 pipeline toolchain.
# ABOUTME: Analyses non-canonical hotspot-SNV peptides; outputs normalised tables and PDFs.

library(tidyverse)
library(RColorBrewer)
library(data.table)


# Build Protein ID → Gene name lookup from UniProt FASTA headers (GN= field)
args <- commandArgs(trailingOnly = TRUE)
proteome_file <- if (length(args) >= 1) args[1] else "human_canonical_proteome.fasta"
headers <- grep("^>", readLines(proteome_file), value = TRUE)
headersplit =   str_split(headers, "\\|")
Protein_ID = rep("",length(headersplit) )
Gene_name = rep("",length(headersplit) )
for(cont in 1:length(headersplit)){
  Protein_ID[cont] = headersplit[[cont]][2]
  Gene_name[cont] = sub(".*GN=([^ ]+).*", "\\1", headersplit[[cont]][3])
}
Id_genematch = data.frame(Gene = Gene_name, Protein_ID = Protein_ID)


# Hotspot entries are distinguished from gene fusions by ":" in the header
# (e.g., Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3)
noncanonical_peptides = data.frame(fread("non_canonical_peptide_headers.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1),, drop = FALSE]


# Extract peptide sequence: penultimate "_"-delimited field (last field is charge state).
# Loop needed because protein IDs may themselves contain underscores.
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}


outputDIANN =data.frame(fread("Reports/report_peptidoforms.pr_matrix.strict.tsv", fill=TRUE), check.names=FALSE)

numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]

# CPM normalisation
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
#normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000
normalizednumericoutputDIANN = numericoutputDIANN

colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones])))
normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN, check.names = FALSE)

metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)


normalizednumericoutputDIANN_selection = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence%in% noncanonical_peptides_sequenceonly,]


# When multiple mutations share a peptide, pick the most frequent one (frequency encoded after ":" in protein ID)
mutationssharingpeptide = str_split(normalizednumericoutputDIANN_selection$Protein.Ids, ";")

for(nshare in seq_along(mutationssharingpeptide)){
  thismutationssharingpeptide = mutationssharingpeptide[[nshare]]
  if(length(thismutationssharingpeptide) > 1){
    mostcommonmut_pos = which.max(as.numeric(str_split_fixed(thismutationssharingpeptide, ":", 2)[,2]))
    mostcommonmut = thismutationssharingpeptide[mostcommonmut_pos]
    normalizednumericoutputDIANN_selection$Protein.Group[nshare] = mostcommonmut
    theidandmut = str_split_fixed(mostcommonmut, "_", 2)
    thegene = Id_genematch$Gene[Id_genematch$Protein_ID == theidandmut[1]]
    normalizednumericoutputDIANN_selection$Genes[nshare] =paste0(thegene,"_",theidandmut[2] )
  }
}


dir.create("Peptidomics_Results", showWarnings = FALSE)
write.table(normalizednumericoutputDIANN_selection, "Peptidomics_Results/hotspot_peptides.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )

# Convert to character so these columns are excluded by is.numeric() downstream
normalizednumericoutputDIANN_selection$Proteotypic = as.character(normalizednumericoutputDIANN_selection$Proteotypic)
normalizednumericoutputDIANN_selection$Precursor.Charge = as.character(normalizednumericoutputDIANN_selection$Precursor.Charge)


normalizednumericoutputDIANN_selection$Gene_and_mut = apply(cbind(normalizednumericoutputDIANN_selection$Genes, normalizednumericoutputDIANN_selection$Stripped.Sequence  ), 1, paste, collapse= "_")

numeric_cols <- sapply(normalizednumericoutputDIANN_selection, is.numeric)
selected_normalizednumericoutputDIANN <- normalizednumericoutputDIANN_selection[, c(names(normalizednumericoutputDIANN_selection)[numeric_cols], "Gene_and_mut")]


myColors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E5C100",
  "#A65628", "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3",
  "#66C2A5", "#0033A0", "#F4A6D7", "#FC8D62", "#8DD3C7", "#FFFFB3",
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
  "#D9D9D9", "#BC80BD"
)


# Find the canonical counterpart peptide for each hotspot mutant.
# Three cases depending on how the mutation affects tryptic cleavage:
#   Branch A: Alt is K/R → mutation created a cleavage site; canonical is longer
#   Branch B: Ref is K/R → mutation destroyed a cleavage site; canonical is shorter
#   Branch C: standard missense → same-length canonical, fuzzy-matched (agrep distance=1)
noncanonicalpeptides = normalizednumericoutputDIANN_selection$Stripped.Sequence

sequencesmatching_samelength_canonical_SNV = c()
Genenames_sequencesmatching_samelength_canonical_SNV = c()
for(cont in 1:length(noncanonicalpeptides)){
  THERef <- NULL
  THElocation <- NULL
  thismut = str_split_fixed( normalizednumericoutputDIANN_selection$Genes[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)

  if((Alt == "R" | Alt == "K") & str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont]) ){
    # Branch A
    mutatedAaremoved = substring(noncanonicalpeptides[cont], 1, (nchar(noncanonicalpeptides[cont])-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    if(length(potentialnonmut) > 0){
      for(npotentialnonmut in 1:length(potentialnonmut)){
        Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
        Mutpept = str_split(mutatedAaremoved, "")[[1]]
        option1 = Mutpept == Refpept[1: length(Mutpept)]
        if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
          THERef= potentialnonmut[npotentialnonmut]
          THElocation = locationofpotentialnonmut[npotentialnonmut]
        }
      }
      if(!is.null(THERef)){
        sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
        Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
      }
    }

  } else if ((Ref == "R" | Ref == "K" )& str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont])  ){
    # Branch B
    if (Alt[[1]] == "*") {
      pattern_to_use = "\\*"
    } else {
      pattern_to_use = Alt[[1]]
    }
    fragmentsofpeptide = str_split(noncanonicalpeptides[cont], pattern_to_use)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    potentialnonmut =potentialnonmut [nchar(potentialnonmut) <  nchar(noncanonicalpeptides[cont])]
    if(length(potentialnonmut) > 0 ){
    potentialnonmut = potentialnonmut[order(nchar(potentialnonmut))]
    for(npotentialnonmut in 1:length(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(noncanonicalpeptides[cont], "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
        break()
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
    }

  } else {
    # Branch C
    potentialmatches = agrep(noncanonicalpeptides[cont], normalizednumericoutputDIANN$Stripped.Sequence,max.distance = 1,fixed = T)
    sequencesmatching= normalizednumericoutputDIANN$Stripped.Sequence[potentialmatches]
    Genenamematching = normalizednumericoutputDIANN$Genes[potentialmatches]
    sequencesmatching_samelength = sequencesmatching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    Genenamematching_samelength = Genenamematching[nchar(sequencesmatching) == nchar(noncanonicalpeptides[cont])]
    sequencesmatching_samelength_canonical = sequencesmatching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    Genenamematching_samelength_canonical = Genenamematching_samelength[!sequencesmatching_samelength %in% noncanonicalpeptides]
    if(length(sequencesmatching_samelength_canonical)!= 0){
      mutthiscase  = normalizednumericoutputDIANN_selection$Genes[cont]
      mutationchangefull = str_split_fixed(mutthiscase, "_", 3)
      mutationchange = c(substring(mutationchangefull[,2],1,1),substring(mutationchangefull[,3],1,1) )
      sequencesmatching_samelength_canonicalfiltered = c()
      Genes_sequencesmatching_samelength_canonicalfiltered = c()
      for(matches in 1:length(sequencesmatching_samelength_canonical)){

        chars1 <- strsplit(sequencesmatching_samelength_canonical[matches], "")[[1]]
        chars2 <- strsplit(noncanonicalpeptides[cont], "")[[1]]

        diff_pos <- which(chars1 != chars2)
        if(chars1[diff_pos] == mutationchange[1] &  chars2[diff_pos] == mutationchange[2]){
          sequencesmatching_samelength_canonicalfiltered = c(sequencesmatching_samelength_canonicalfiltered, sequencesmatching_samelength_canonical[matches])
          Genes_sequencesmatching_samelength_canonicalfiltered = c(Genes_sequencesmatching_samelength_canonicalfiltered, Genenamematching_samelength_canonical[matches])

        }
      }
      sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, sequencesmatching_samelength_canonicalfiltered)
      Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, Genes_sequencesmatching_samelength_canonicalfiltered)
    }
  }
}


canonicalpeptidesfromSNV = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence %in% sequencesmatching_samelength_canonical_SNV,]
canonicalpeptidesfromSNV$Gene_and_mut = apply(cbind(canonicalpeptidesfromSNV$Genes, canonicalpeptidesfromSNV$Stripped.Sequence  ), 1, paste, collapse= "_")

canonicalpeptidesfromSNV$Proteotypic = as.character(canonicalpeptidesfromSNV$Proteotypic)
canonicalpeptidesfromSNV$Precursor.Charge = as.character(canonicalpeptidesfromSNV$Precursor.Charge)

numeric_cols <- sapply(canonicalpeptidesfromSNV, is.numeric)
selected_canonicalpeptidesfromSNV <- canonicalpeptidesfromSNV[, c(names(canonicalpeptidesfromSNV)[numeric_cols], "Gene_and_mut")]

selected_normalizednumericoutputDIANN$Canon = FALSE
selected_canonicalpeptidesfromSNV$Canon = TRUE

noncanonandcanon = rbind(selected_normalizednumericoutputDIANN, selected_canonicalpeptidesfromSNV)


noncanonandcanon$Label = str_split_fixed(noncanonandcanon$Gene_and_mut, ";",2)[,1]
Labelcanon = str_split_fixed(noncanonandcanon$Label, "_", 2)[,1]
noncanonandcanon$Label[noncanonandcanon$Canon] = Labelcanon[noncanonandcanon$Canon]


# Deduplicate: keep the row with the highest total intensity per label (dominant charge state)
numeric_cols <- names(noncanonandcanon)[sapply(noncanonandcanon, is.numeric)]
filtered_df <- data.frame(noncanonandcanon %>%
                            rowwise() %>%
                            mutate(Total = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
                            group_by(Gene_and_mut, Canon, Label) %>%
                            filter(Total == max(Total, na.rm = TRUE)) %>%
                            ungroup() %>%
                            select(-Total),
                          check.names = FALSE
)
noncanonandcanon = filtered_df
noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
meltselected_normalizednumericoutputDIANN = data.frame(pivot_longer(noncanonandcanon, cols = where(is.numeric), names_to = "variable", values_to = "value"), check.names = FALSE)
meltselected_normalizednumericoutputDIANN$variable = factor(meltselected_normalizednumericoutputDIANN$variable , levels = unique(meltselected_normalizednumericoutputDIANN$variable)[length(unique(meltselected_normalizednumericoutputDIANN$variable )):1])


noncanonandcanon$Sequence = ""
for(cont in 1:length(sequencesmatching_samelength_canonical_SNV)){
  whichone = agrep(Genenames_sequencesmatching_samelength_canonical_SNV[cont], noncanonandcanon$Gene_and_mut, max.distance = 1)
  noncanonandcanon$Sequence[whichone] = sequencesmatching_samelength_canonical_SNV[cont]
}
noncanonandcanon$Sequence[noncanonandcanon$Canon == "FALSE"] = ""

noncanonandcanon = noncanonandcanon[,!colnames(noncanonandcanon) %in% c("Proteotypic", "Precursor.Charge")]


# PDF 1: one plot per gene (all mutations for that gene)
pdf("Peptidomics_Results/hotspot_by_gene.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))
noncanonLabel = noncanonLabel[!duplicated(str_split_fixed(noncanonLabel,"_",2)[,1])]

noncanonandcanon$Label[noncanonandcanon$Canon == TRUE] = paste0(noncanonandcanon$Label[noncanonandcanon$Canon == TRUE], "_", noncanonandcanon$Sequence[noncanonandcanon$Canon == TRUE] )

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  thecanon = noncanonandcanon[agrep(the_sequence, noncanonandcanon$Sequence, 1),]
  thenoncanon = noncanonandcanon[agrep(the_sequence, thenoncanonsequences, 1),]
  thenoncanon2 = noncanonandcanon[grep(the_name, thenoncanonnames, 1),]
  thisselectionandcanon = unique(rbind(thenoncanon,thenoncanon2 , thecanon))
  thisselectionandcanon$Label = as.character(thisselectionandcanon$Label)
  meltselected_normalizednumericoutputDIANN_thisprot = data.frame(pivot_longer(thisselectionandcanon, cols = where(is.numeric), names_to = "variable", values_to = "value"), check.names = FALSE)
  meltselected_normalizednumericoutputDIANN_thisprot$variable = factor(meltselected_normalizednumericoutputDIANN_thisprot$variable , levels = unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)[length(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable )):1])

  # dummy row keeps "Mutated" shape in the legend even when no mutant rows appear in this subset
  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
    Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],
    Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],
    value = NA

  )
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot, dummy_row)
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),]
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3, na.rm = TRUE)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Raw Intensity")+
    scale_color_manual(values = myColors) + ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}
dev.off()


# PDF 2: one plot per individual mutation
pdf("Peptidomics_Results/hotspot_by_mutation.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))

for(cont in 1:length(noncanonLabel)){
  THERef <- NULL
  THElocation <- NULL
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  thismut = str_split_fixed( thisselection$Gene_and_mut[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)

  if(Alt == "R" | Alt == "K" ){
    # Branch A
    mutatedAaremoved = substring(the_sequence, 1, (nchar(the_sequence)-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]
    if(length(potentialnonmut) > 0){
      for(npotentialnonmut in 1:length(potentialnonmut)){
        Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
        Mutpept = str_split(mutatedAaremoved, "")[[1]]
        option1 = Mutpept == Refpept[1: length(Mutpept)]
        if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
          THERef= potentialnonmut[npotentialnonmut]
          THElocation = locationofpotentialnonmut[npotentialnonmut]
        }
      }
    }
    if(!is.null(THERef)){
      thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]
    } else {
      thecanon = noncanonandcanon[0,]
    }

  } else if (Ref == "R" | Ref == "K" ){
    # Branch B
    fragmentsofpeptide = str_split(the_sequence, Alt)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]

    if(length(potentialnonmut) > 0){
      for(npotentialnonmut in 1:length(potentialnonmut)){
        Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
        Mutpept = str_split(the_sequence, "")[[1]]
        option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
        option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
        if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
          THERef= potentialnonmut[npotentialnonmut]
          THElocation = locationofpotentialnonmut[npotentialnonmut]
        }
      }
    }
    if(!is.null(THERef)){
      thecanon = noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]
    } else {
      thecanon = noncanonandcanon[0,]
    }

  } else {
    # Branch C
    thecanon = noncanonandcanon[agrep(the_sequence, noncanonandcanon$Sequence, 1),]
  }
  thenoncanon = noncanonandcanon[grep(the_sequence, thenoncanonsequences, 1),]

  meltselected_normalizednumericoutputDIANN_thisprot_onemut = data.frame(pivot_longer(rbind(thecanon, thenoncanon), cols = where(is.numeric), names_to = "variable", values_to = "value"), check.names = FALSE)

  dummy_row <- data.frame(
    Gene_and_mut = "",
    Canon = factor("TRUE", levels = c("TRUE", "FALSE")),
    Label =meltselected_normalizednumericoutputDIANN_thisprot_onemut$Label[1],
    Sequence = "",
    variable = meltselected_normalizednumericoutputDIANN_thisprot_onemut$variable[1],
    value = NA

  )
  plot_data <- rbind(meltselected_normalizednumericoutputDIANN_thisprot_onemut, dummy_row)
  plot_data$Label = factor(
    plot_data$Label,
    levels = unique(plot_data$Label)[order(lengths(regmatches(unique(plot_data$Label), gregexpr("_", unique(plot_data$Label)))))] )
  plot_data$Status = "Mutated"
  plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
  plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
  plot_data= plot_data[order(plot_data$Status ),]
  plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))
  plot_data$variable = factor(plot_data$variable , levels = unique(plot_data$variable)[length(unique(plot_data$variable )):1])
  title = unique(gsub("_", " ", plot_data$Label))[which.max(nchar(unique(gsub("_", " ", plot_data$Label))))]
  p = ggplot(plot_data, aes(x = variable, y = value , color = Label, shape = Status))   +
    geom_point(size = 3, na.rm = TRUE)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Raw Intensity")+
    scale_color_manual(values = myColors) + ggtitle(title) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}

dev.off()


write.table(noncanonandcanon, "Peptidomics_Results/hotspot_peptides_with_canonical.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
