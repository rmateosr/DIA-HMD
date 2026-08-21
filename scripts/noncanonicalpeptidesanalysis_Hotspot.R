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


# Resolve the gene label of a wild-type counterpart peptide.
# DIA-NN's Genes field names only ONE protein a peptide is consistent with, so it mislabels in two
# ways. A peptide shared between paralogs gets whichever one DIA-NN picked: IGDFGLATVK is
# consistent with ARAF, BRAF and RAF1, and was reported as RAF1, so the BRAF V600E panel labelled
# its own wild-type counterpart RAF1. And a wild-type peptide can be reported under a *different*
# hotspot variant of the same gene, because that variant entry contains the peptide unchanged:
# unmutated LVVVGAGGVGK carries Genes = "KRAS_K117_R:1". Going through Protein.Ids and preferring
# the gene of the mutation the panel is about fixes both, and returns a bare gene symbol so the
# label no longer depends on stripping a variant suffix.
resolve_wt_gene = function(protein_ids, genes_field, preferred_gene) {
  accs = str_split(protein_ids, ";")[[1]]
  accs = sub("-.*$", "", str_split_fixed(accs, "_", 2)[, 1])
  candidates = unique(Id_genematch$Gene[match(accs, Id_genematch$Protein_ID)])
  candidates = candidates[!is.na(candidates)]
  if (length(preferred_gene) == 1 && !is.na(preferred_gene) &&
      preferred_gene %in% candidates) {
    return(preferred_gene)
  }
  if (length(candidates) >= 1) return(candidates[1])
  str_split_fixed(genes_field, "_", 2)[, 1]
}


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


# The .strict matrix has variant-peptide cells blanked where that run failed the Q.Value gate,
# the replicate requirement or fragment geometry, so a mutation the results table rejected cannot
# be plotted here. Wild-type counterpart rows are not gated, so the canonical matching below is
# unaffected.
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

# scale_color_manual aborts rather than recycling when a panel has more series than the
# palette has values, which kills the whole script. Append interpolated colours so no panel
# can exhaust it, leaving the 26 hand-picked values in front and used first.
myColors <- c(myColors, grDevices::colorRampPalette(myColors)(256))


# Find the canonical counterpart peptide for each hotspot mutant.
# Three cases depending on how the mutation affects tryptic cleavage:
#   Branch A: Alt is K/R → mutation created a cleavage site; canonical is longer
#   Branch B: Ref is K/R → mutation destroyed a cleavage site; canonical is shorter
#   Branch C: standard missense → same-length canonical, fuzzy-matched (agrep distance=1)
noncanonicalpeptides = normalizednumericoutputDIANN_selection$Stripped.Sequence

sequencesmatching_samelength_canonical_SNV = c()
Genenames_sequencesmatching_samelength_canonical_SNV = c()
# Gene of the mutation each wild-type counterpart belongs to, kept parallel to the two vectors
# above. Used only to label the counterpart; it is deliberately NOT fed to the agrep matching
# below, which needs the original Genes strings to stay long enough to match specifically.
Mutgenes_sequencesmatching_samelength_canonical_SNV = c()
# The mutant peptide each wild-type counterpart was found for, same parallel structure. Used to
# discard a counterpart whose own mutant peptide did not survive the filters.
Mutpeptides_sequencesmatching_samelength_canonical_SNV = c()
for(cont in 1:length(noncanonicalpeptides)){
  thismut = str_split_fixed( normalizednumericoutputDIANN_selection$Genes[cont], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)
  mutgene = thismut[1]

  if((Alt == "R" | Alt == "K") & str_locate_all(noncanonicalpeptides[cont], "[KR]")[[1]][,1][1] == nchar(noncanonicalpeptides[cont]) ){
    # Branch A
    mutatedAaremoved = substring(noncanonicalpeptides[cont], 1, (nchar(noncanonicalpeptides[cont])-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, normalizednumericoutputDIANN$Stripped.Sequence)
    potentialnonmut = normalizednumericoutputDIANN$Stripped.Sequence[locationofpotentialnonmut]
    # NA and not NULL, and reset per mutant peptide: all four vectors below must gain exactly one
    # element each so they stay index-aligned, and c(x, NULL) would append nothing. A NULL left
    # over from the previous iteration would instead attach the previous peptide's counterpart to
    # this one. An NA entry matches no sequence downstream, so it is inert.
    THERef = NA_character_
    THElocation = NA_integer_
    for(npotentialnonmut in seq_along(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, THERef)
    Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, normalizednumericoutputDIANN$Genes[THElocation])
    Mutgenes_sequencesmatching_samelength_canonical_SNV = c(Mutgenes_sequencesmatching_samelength_canonical_SNV, mutgene)
    Mutpeptides_sequencesmatching_samelength_canonical_SNV = c(Mutpeptides_sequencesmatching_samelength_canonical_SNV, noncanonicalpeptides[cont])

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
    # Same reasoning as Branch A: an NA sentinel reset here keeps the four vectors aligned and
    # cannot inherit the previous peptide's counterpart.
    THERef = NA_character_
    THElocation = NA_integer_
    if(length(potentialnonmut) > 0 ){
    potentialnonmut = potentialnonmut[order(nchar(potentialnonmut))]
    for(npotentialnonmut in seq_along(potentialnonmut)){
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
    Mutgenes_sequencesmatching_samelength_canonical_SNV = c(Mutgenes_sequencesmatching_samelength_canonical_SNV, mutgene)
    Mutpeptides_sequencesmatching_samelength_canonical_SNV = c(Mutpeptides_sequencesmatching_samelength_canonical_SNV, noncanonicalpeptides[cont])
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
      Mutgenes_sequencesmatching_samelength_canonicalfiltered = c()
      Mutpeptides_sequencesmatching_samelength_canonicalfiltered = c()
      for(matches in seq_along(sequencesmatching_samelength_canonical)){

        chars1 <- strsplit(sequencesmatching_samelength_canonical[matches], "")[[1]]
        chars2 <- strsplit(noncanonicalpeptides[cont], "")[[1]]

        diff_pos <- which(chars1 != chars2)
        if(chars1[diff_pos] == mutationchange[1] &  chars2[diff_pos] == mutationchange[2]){
          sequencesmatching_samelength_canonicalfiltered = c(sequencesmatching_samelength_canonicalfiltered, sequencesmatching_samelength_canonical[matches])
          Genes_sequencesmatching_samelength_canonicalfiltered = c(Genes_sequencesmatching_samelength_canonicalfiltered, Genenamematching_samelength_canonical[matches])
          Mutgenes_sequencesmatching_samelength_canonicalfiltered = c(Mutgenes_sequencesmatching_samelength_canonicalfiltered, mutgene)
          Mutpeptides_sequencesmatching_samelength_canonicalfiltered = c(Mutpeptides_sequencesmatching_samelength_canonicalfiltered, noncanonicalpeptides[cont])

        }
      }
      sequencesmatching_samelength_canonical_SNV = c(sequencesmatching_samelength_canonical_SNV, sequencesmatching_samelength_canonicalfiltered)
      Genenames_sequencesmatching_samelength_canonical_SNV = c(Genenames_sequencesmatching_samelength_canonical_SNV, Genes_sequencesmatching_samelength_canonicalfiltered)
      Mutgenes_sequencesmatching_samelength_canonical_SNV = c(Mutgenes_sequencesmatching_samelength_canonical_SNV, Mutgenes_sequencesmatching_samelength_canonicalfiltered)
      Mutpeptides_sequencesmatching_samelength_canonical_SNV = c(Mutpeptides_sequencesmatching_samelength_canonical_SNV, Mutpeptides_sequencesmatching_samelength_canonicalfiltered)
    }
  }
}


canonicalpeptidesfromSNV = normalizednumericoutputDIANN[normalizednumericoutputDIANN$Stripped.Sequence %in% sequencesmatching_samelength_canonical_SNV,]
canonicalpeptidesfromSNV$Gene_and_mut = apply(cbind(canonicalpeptidesfromSNV$Genes, canonicalpeptidesfromSNV$Stripped.Sequence  ), 1, paste, collapse= "_")

# Carried as its own column rather than folded into Genes or Gene_and_mut, because those strings
# are the pattern and the target of the agrep matching further down and must keep their original
# form. Only the panel label reads this.
canonicalpeptidesfromSNV$WT_Gene = mapply(
  resolve_wt_gene,
  canonicalpeptidesfromSNV$Protein.Ids,
  canonicalpeptidesfromSNV$Genes,
  Mutgenes_sequencesmatching_samelength_canonical_SNV[
    match(canonicalpeptidesfromSNV$Stripped.Sequence,
          sequencesmatching_samelength_canonical_SNV)],
  USE.NAMES = FALSE)

# Which mutant peptides this counterpart belongs to, so an orphan can be recognised later.
canonicalpeptidesfromSNV$WT_Partners = sapply(
  canonicalpeptidesfromSNV$Stripped.Sequence,
  function(s) paste(unique(Mutpeptides_sequencesmatching_samelength_canonical_SNV[
      which(sequencesmatching_samelength_canonical_SNV == s)]), collapse = ";"),
  USE.NAMES = FALSE)

canonicalpeptidesfromSNV$Proteotypic = as.character(canonicalpeptidesfromSNV$Proteotypic)
canonicalpeptidesfromSNV$Precursor.Charge = as.character(canonicalpeptidesfromSNV$Precursor.Charge)

numeric_cols <- sapply(canonicalpeptidesfromSNV, is.numeric)
selected_canonicalpeptidesfromSNV <- canonicalpeptidesfromSNV[, c(names(canonicalpeptidesfromSNV)[numeric_cols], "Gene_and_mut", "WT_Gene", "WT_Partners")]

selected_normalizednumericoutputDIANN$WT_Gene = ""
selected_normalizednumericoutputDIANN$WT_Partners = ""
selected_normalizednumericoutputDIANN$Mut_Sequence = normalizednumericoutputDIANN_selection$Stripped.Sequence
selected_canonicalpeptidesfromSNV$Mut_Sequence = ""
selected_normalizednumericoutputDIANN$Canon = FALSE
selected_canonicalpeptidesfromSNV$Canon = TRUE

noncanonandcanon = rbind(selected_normalizednumericoutputDIANN, selected_canonicalpeptidesfromSNV)


noncanonandcanon$Label = str_split_fixed(noncanonandcanon$Gene_and_mut, ";",2)[,1]
Labelcanon = str_split_fixed(noncanonandcanon$Label, "_", 2)[,1]
# Prefer the resolved gene, so a counterpart shared between paralogs is named after the mutation
# whose panel it appears in rather than after whichever paralog DIA-NN happened to report.
has_resolved = nzchar(noncanonandcanon$WT_Gene) & !is.na(noncanonandcanon$WT_Gene)
Labelcanon[has_resolved] = noncanonandcanon$WT_Gene[has_resolved]
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

noncanonandcanon = noncanonandcanon[,!colnames(noncanonandcanon) %in% c("Proteotypic", "Precursor.Charge", "WT_Gene")]


# Drop mutant peptides with no surviving measurement. gate_variant_cells.py empties every cell of a
# variant precursor that fails the q-value, replicate or fragment-geometry filter, so these rows
# would draw an empty panel -- and the results table reports no call for them either. Wild-type rows
# are kept whatever their values: they are reference context, not detection claims.
mutant_numeric = names(noncanonandcanon)[sapply(noncanonandcanon, is.numeric)]
has_measurement = rowSums(!is.na(noncanonandcanon[, mutant_numeric, drop = FALSE])) > 0
drop_empty = (noncanonandcanon$Canon == "FALSE") & !has_measurement
cat("Mutant peptide rows dropped for having no surviving measurement:", sum(drop_empty),
    "of", sum(noncanonandcanon$Canon == "FALSE"), "\n")
noncanonandcanon = noncanonandcanon[!drop_empty, ]


# Keep a wild-type counterpart only while at least one of the mutant peptides it was found for is
# still plotted. Otherwise the panel shows the wild type of a mutation it no longer displays.
surviving_mut = unique(noncanonandcanon$Mut_Sequence[noncanonandcanon$Canon == "FALSE"])
surviving_mut = surviving_mut[nzchar(surviving_mut)]
is_canon = noncanonandcanon$Canon == "TRUE"
keep_row = rep(TRUE, nrow(noncanonandcanon))
keep_row[is_canon] = vapply(
  strsplit(noncanonandcanon$WT_Partners[is_canon], ";"),
  function(p) any(p %in% surviving_mut), logical(1))
cat("Wild-type counterparts dropped as orphaned:", sum(!keep_row),
    "of", sum(is_canon), "\n")
noncanonandcanon = noncanonandcanon[keep_row, ]
noncanonandcanon = noncanonandcanon[, !colnames(noncanonandcanon) %in% c("WT_Partners", "Mut_Sequence")]


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
    geom_point(size = 3, na.rm = TRUE)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Peptide Signal Intensity")+
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = myColors) + ggtitle(str_split_fixed(noncanonLabel[cont], "_",2)[1]) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}
dev.off()


# PDF 2: one plot per individual mutation
pdf("Peptidomics_Results/hotspot_by_mutation.pdf", width = 10, height = 15)

noncanonLabel = as.character(unique(noncanonandcanon$Label[noncanonandcanon$Canon == FALSE]))

for(cont in 1:length(noncanonLabel)){
  noncanonandcanon$Canon = factor(noncanonandcanon$Canon, c("TRUE", "FALSE"))
  thisselection = noncanonandcanon[agrep(noncanonLabel[cont],noncanonandcanon$Label),]
  the_sequence = (str_split_fixed(thisselection$Label, "_", 4)[,4])[1]
  the_name = (str_split_fixed(thisselection$Label, "_", 4)[,1])[1]
  thenoncanonsequences =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,4])
  thenoncanonnames =  (str_split_fixed(noncanonandcanon$Label, "_", 4)[,1])
  # thisselection holds the rows for THIS label, so the panel's mutation is its first row.
  # Indexing it by the outer loop counter read a different mutation's annotation, or NA once
  # cont exceeded the subset's row count.
  thismut = str_split_fixed( thisselection$Gene_and_mut[1], "_",3)
  Ref = substring(thismut[2],1,1)
  Alt = substring(thismut[3],1,1)

  if(Alt == "R" | Alt == "K" ){
    # Branch A
    mutatedAaremoved = substring(the_sequence, 1, (nchar(the_sequence)-1))
    locationofpotentialnonmut = grep(mutatedAaremoved, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]
    # seq_along, not 1:length: an empty candidate list makes 1:length() count 1,0 and index 0
    # returns character(0), which str_split()[[1]] cannot subscript. THERef is reset each panel
    # so a panel with no match cannot inherit the previous panel's counterpart.
    THERef = NA_character_
    for(npotentialnonmut in seq_along(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(mutatedAaremoved, "")[[1]]
      option1 = Mutpept == Refpept[1: length(Mutpept)]
      if((sum(option1)== length(option1)) &  (Refpept[length(Mutpept) + 1] == Ref)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = if (is.na(THERef)) noncanonandcanon[0, ] else noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]

  } else if (Ref == "R" | Ref == "K" ){
    # Branch B
    fragmentsofpeptide = str_split(the_sequence, Alt)[[1]]
    longestfragment  =fragmentsofpeptide[which.max(nchar(fragmentsofpeptide))]
    locationofpotentialnonmut = grep(longestfragment, noncanonandcanon$Sequence)
    potentialnonmut = noncanonandcanon$Sequence[locationofpotentialnonmut]

    THERef = NA_character_
    for(npotentialnonmut in seq_along(potentialnonmut)){
      Refpept = str_split(potentialnonmut[npotentialnonmut], "")[[1]]
      Mutpept = str_split(the_sequence, "")[[1]]
      option1 = Mutpept[1:(length(Refpept)-1)] == Refpept[1:length(Refpept)-1]
      option2 = Mutpept[(length(Mutpept) - length(Refpept) + 1):length(Mutpept)] == Refpept
      if(sum(option1)  == length(option1) | sum(option2)  == length(option2)){
        THERef= potentialnonmut[npotentialnonmut]
        THElocation = locationofpotentialnonmut[npotentialnonmut]
      }
    }
    thecanon = if (is.na(THERef)) noncanonandcanon[0, ] else noncanonandcanon[grep(THERef, noncanonandcanon$Sequence),]

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
    geom_point(size = 3, na.rm = TRUE)+ theme_minimal() + coord_flip()  +xlab("Cell Line") + ylab("Peptide Signal Intensity")+
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = myColors) + ggtitle(title) +
    scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE)

  print(p)

}

dev.off()


write.table(noncanonandcanon, "Peptidomics_Results/hotspot_peptides_with_canonical.tsv", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE )
