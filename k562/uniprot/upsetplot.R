library(data.table)
library(dplyr)
library(UpSetR)
library(ggplot2)


parsimony = fread("k562/uniprot/data/comparison/parsimony_up.csv")

inclusion = fread("k562/uniprot/data/comparison/inclusion_up.csv")

exclusion = fread("k562/uniprot/data/comparison/exclusion_up.csv")

peptidecentric = fread("k562/uniprot/data/comparison/peptide_centric_up.csv")


parsimony <- parsimony[parsimony$Q_Value<=0.01,]
inclusion <- inclusion[inclusion$Q_Value<=0.01,]
exclusion <- exclusion[exclusion$Q_Value<=0.01,]
peptidecentric <- peptidecentric[peptidecentric$Q_Value<=0.01,]

peptidecentric <- peptidecentric %>%
  dplyr::mutate(Protein = strsplit(Protein, ";")) %>%
  tidyr::unnest(Protein)

parsimony_proteins <- parsimony$Protein
inclusion_proteins <- inclusion$Protein
exclusion_proteins <- exclusion$Protein
peptidecentric_proteins <- peptidecentric$Protein


peptidecentric_proteins <- unique(peptidecentric_proteins)


pars <- length(parsimony_proteins)
incl <- length(inclusion_proteins)
exc <- length(exclusion_proteins)
pepc <- length(peptidecentric_proteins)

pars_incl <- length(intersect(parsimony_proteins, inclusion_proteins))
pars_exc <- length(intersect(parsimony_proteins, exclusion_proteins))
pars_pepc <- length(intersect(parsimony_proteins, peptidecentric_proteins))
incl_exc <- length(intersect(inclusion_proteins, exclusion_proteins))
incl_pepc <- length(intersect(inclusion_proteins, peptidecentric_proteins))
excl_pepc <- length(intersect(exclusion_proteins, peptidecentric_proteins))

pars_incl_exc <- length(intersect(intersect(parsimony_proteins, inclusion_proteins), exclusion_proteins))
incl_exc_pepc <- length(intersect(intersect(inclusion_proteins, exclusion_proteins), peptidecentric_proteins))
pars_exc_pepc <- length(intersect(intersect(parsimony_proteins, exclusion_proteins), peptidecentric_proteins))
pars_incl_pepc <- length(intersect(intersect(parsimony_proteins, inclusion_proteins), peptidecentric_proteins))
incl_exc_pars_pepc <- length(intersect(intersect(intersect(inclusion_proteins, exclusion_proteins), peptidecentric_proteins), parsimony_proteins))


setdiff_pars_incl <- length(setdiff(parsimony_proteins, inclusion_proteins))
setdiff_pars_exc <- length(setdiff(exclusion_proteins, parsimony_proteins))
setdiff_pars_pepc <- length(setdiff(parsimony_proteins, peptidecentric_proteins))
setdiff_incl_exc <- length(setdiff(exclusion_proteins, inclusion_proteins))
setdiff_pepc_incl <- length(setdiff(peptidecentric_proteins, inclusion_proteins))
setdiff_excl_pepc <- length(setdiff(exclusion_proteins, peptidecentric_proteins))



# Dataset
input <- c(
  exclusion = exc,
  parsimony = pars,
  peptide.centric = pepc,
  inclusion = incl,
  "exclusion&parsimony" = pars_exc,
  "exclusion&peptide.centric" = excl_pepc,
  "exclusion&inclusion" = incl_exc,
  "parsimony&peptide.centric" = pars_pepc,
  "parsimony&inclusion" = pars_incl,
  "inclusion&peptide.centric" = incl_pepc,
  "inclusion&peptide.centric&parsimony" = pars_incl_pepc,
  "exclusion&peptide.centric&parsimony" = pars_exc_pepc,
  "exclusion&inclusion&parsimony" = pars_incl_exc,
  "exclusion&inclusion&peptide.centric" = incl_exc_pepc,
  "exclusion&inclusion&parsimony&peptide.centric" = incl_exc_pars_pepc
)



pdf("k562/uniprot/plots/upsetplot_up.pdf", width=18, height=12)


# Plot
upset(fromExpression(input),
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 3, 
      point.size = 4, 
      line.size = 2)


dev.off()

