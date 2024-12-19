
from k562.uniprot import file_converter_helper

from matplotlib_venn import venn2, venn2_circles
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt



inhouse_data_incl = file_converter_helper.convert_internal_pi(filename="k562/uniprot/data/comparison/inclusion_up.csv")

inhouse_data_pars = file_converter_helper.convert_internal_pi(filename="k562/uniprot/data/comparison/parsimony_up.csv")

inhouse_data_exc = file_converter_helper.convert_internal_pi(filename="k562/uniprot/data/comparison/exclusion_up.csv")

new_parsimony_data = []
for pars in inhouse_data_pars:
    updated_prot = file_converter_helper.convert_pyproteininference_ids_for_pia(pars[0])
    pars[0] = updated_prot
    new_parsimony_data.append(pars)

new_inclusion_data = []
for incl in inhouse_data_incl:
    updated_prot = file_converter_helper.convert_pyproteininference_ids_for_pia(incl[0])
    incl[0] = updated_prot
    new_inclusion_data.append(incl)

# PIA Inclusion
pia_inclusion = file_converter_helper.convert_pia_old(filename="k562/uniprot/data/comparison/pia_knime_inclusion_up.txt", ftype="inclusion")

# PIA Parsimony
inhouse_parsimony_proteins = [x[0] for x in new_parsimony_data]
pia_parsimony = file_converter_helper.convert_pia_new(filename="k562/uniprot/data/comparison/pia_knime_parsimony_up.txt", ftype="parsimony", pyproteininference_proteins=inhouse_parsimony_proteins)


# Percolator PI
perc_pi_data = file_converter_helper.convert_perc_pi(target_filename="k562/uniprot/data/comparison/target_protein_up.txt",
                                                     decoy_filename="k562/uniprot/data/comparison/decoy_protein_up.txt")

percpi_data_with_qvalues = file_converter_helper.calc_qvalues(perc_pi_data, decoy_symbol="DECOY_")
exc_data_with_qvalues = file_converter_helper.calc_qvalues(inhouse_data_exc, decoy_symbol="DECOY_")
inhouse_incl_data_with_qvalues = file_converter_helper.calc_qvalues(new_inclusion_data, decoy_symbol="DECOY_")
inhouse_pars_data_with_qvalues = file_converter_helper.calc_qvalues(new_parsimony_data, decoy_symbol="DECOY_")
pia_pars_data_with_qvalues = file_converter_helper.calc_qvalues(pia_parsimony, decoy_symbol="DECOY_")
pia_incl_data_with_qvalues = file_converter_helper.calc_qvalues(pia_inclusion, decoy_symbol="DECOY_")


#### PARSIMONY ####
#### PARSIMONY ####
#### PARSIMONY ####

pdf_file_parsimony = 'k562/uniprot/plots2/venn_proteins_uniprot_parsimony.pdf'

labels_parsimony = ["PIA Occams Razor", "PyProteinInference Parsimony"]

pia_pars_prots = set([x[0] for x in pia_pars_data_with_qvalues if x[-1]<0.01])
inhouse_pars_data_with_qvalues_prots = set([x[0] for x in inhouse_pars_data_with_qvalues if x[-1]<0.01])

# protein_prophet_prots = set([x[0].split(" ")[0]  for x in protein_prophet_data if "##" not in x])
# inhouse_pars_pulp_prots = set([x[0] for x in inhouse_data_pars_glpk if "##" not in x])

parsimony_intersection = pia_pars_prots.intersection(inhouse_pars_data_with_qvalues_prots)
pp_ven = PdfPages(pdf_file_parsimony)

# Subset sizes
s = (
   len(pia_pars_prots-inhouse_pars_data_with_qvalues_prots),  # Ab
   len(inhouse_pars_data_with_qvalues_prots-pia_pars_prots),  # aB
   len(parsimony_intersection),  # AB
)

v = venn2(subsets=s, set_labels=('PIA Occams Razor    ', '    PyProteinInferece'))

# Subset labels
lbl1 = v.get_label_by_id('10')
lbl2 = v.get_label_by_id('01')
lbl3 = v.get_label_by_id('11')

lbl1.set_text(str(len(pia_pars_prots-inhouse_pars_data_with_qvalues_prots)))
x1, y1 = lbl1.get_position()
lbl1.set_position((x1-.2, y1+0))

lbl2.set_text(str(len(inhouse_pars_data_with_qvalues_prots-pia_pars_prots)))
x2, y2 = lbl2.get_position()
lbl2.set_position((x1+1.25, y1+0))

lbl3.set_text(str(len(parsimony_intersection)))

# Subset colors
v.get_patch_by_id('10').set_color('#97e0ff')
v.get_patch_by_id('01').set_color('#e9d848')
v.get_patch_by_id('11').set_color('#8abe48')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.7)
v.get_patch_by_id('01').set_alpha(0.7)
v.get_patch_by_id('11').set_alpha(0.7)

#Border styles
c = venn2_circles(subsets=s, linestyle='solid')
for text in v.set_labels:
    text.set_fontsize(20)
for x in range(len(v.subset_labels)):
    if v.subset_labels[x] is not None:
        v.subset_labels[x].set_fontsize(20)

plt.title("Parsimony Protein Comparison 1% FDR")
pp_ven.savefig()
plt.show()
plt.close()
pp_ven.close()

#### EXCLUSION ####
#### EXCLUSION ####
#### EXCLUSION ####


pdf_file_exclusion = 'k562/uniprot/plots2/venn_proteins_uniprot_exclusion.pdf'

perc_pi_prots = set([x[0] for x in percpi_data_with_qvalues  if x[-1]<0.01])
inhouse_exc_prots = set([x[0]  for x in exc_data_with_qvalues  if x[-1]<0.01])

# perc_pi_prots = set([x[0] for x in perc_pi_data if "##" not in x])
# inhouse_exc_prots = set([x[0] for x in inhouse_data_exc if "##" not in x])


exclusion_intersection = perc_pi_prots.intersection(inhouse_exc_prots)
pp_ven = PdfPages(pdf_file_exclusion)

# Subset sizes
s = (
   len(perc_pi_prots-inhouse_exc_prots),  # Ab
   len(inhouse_exc_prots-perc_pi_prots),  # aB
   len(exclusion_intersection),  # AB
)

v = venn2(subsets=s, set_labels=('Percolator Protein \n Inference    ', '    PyProteinInferece'))

# Subset labels
lbl1 = v.get_label_by_id('10')
lbl2 = v.get_label_by_id('01')
lbl3 = v.get_label_by_id('11')

lbl1.set_text(str(len(perc_pi_prots-inhouse_exc_prots)))
x1, y1 = lbl1.get_position()
lbl1.set_position((x1-.2, y1+0))

lbl2.set_text(str(len(inhouse_exc_prots-perc_pi_prots)))
x2, y2 = lbl2.get_position()
lbl2.set_position((x1+1.25, y1+0))

lbl3.set_text(str(len(exclusion_intersection)))

# Subset colors
v.get_patch_by_id('10').set_color('#97e0ff')
v.get_patch_by_id('01').set_color('#e9d848')
v.get_patch_by_id('11').set_color('#8abe48')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.7)
v.get_patch_by_id('01').set_alpha(0.7)
v.get_patch_by_id('11').set_alpha(0.7)

#Border styles
c = venn2_circles(subsets=s, linestyle='solid')
for text in v.set_labels:
    text.set_fontsize(20)
for x in range(len(v.subset_labels)):
    if v.subset_labels[x] is not None:
        v.subset_labels[x].set_fontsize(20)

plt.title("Exclusion Protein Comparison 1% FDR")
pp_ven.savefig()
plt.show()
plt.close()
pp_ven.close()


#### INCLUSION ####
#### INCLUSION ####
#### INCLUSION ####


pdf_file_incl = 'k562/uniprot/plots2/venn_proteins_uniprot_inclusion.pdf'

pia_incl_prots = set([x[0] for x in pia_incl_data_with_qvalues if x[-1]<0.01])
inhouse_incl_prots = set([x[0] for x in inhouse_incl_data_with_qvalues  if x[-1]<0.01])

# fido_prots = set([x[0] for x in fido_data if "##" not in x[0] and x[2]<0.01])
# inhouse_incl_prots = set([x[0] for x in inhouse_data_incl if "##" not in x[0] and x[2]<0.01])

inclusion_intersection = pia_incl_prots.intersection(inhouse_incl_prots)
pp_ven = PdfPages(pdf_file_incl)

# Subset sizes
s = (
   len(pia_incl_prots-inhouse_incl_prots),  # Ab
   len(inhouse_incl_prots-pia_incl_prots),  # aB
   len(inclusion_intersection),  # AB
)

v = venn2(subsets=s, set_labels=('PIA Report All    ', '    PyProteinInferece'))

# Subset labels
lbl1 = v.get_label_by_id('10')
lbl2 = v.get_label_by_id('01')
lbl3 = v.get_label_by_id('11')

lbl1.set_text(str(len(pia_incl_prots-inhouse_incl_prots)))
x1, y1 = lbl1.get_position()
lbl1.set_position((x1-.2, y1+0))

lbl2.set_text(str(len(inhouse_incl_prots-pia_incl_prots)))
x2, y2 = lbl2.get_position()
lbl2.set_position((x1+1.25, y1+0))

lbl3.set_text(str(len(inclusion_intersection)))

# Subset colors
v.get_patch_by_id('10').set_color('#97e0ff')
v.get_patch_by_id('01').set_color('#e9d848')
v.get_patch_by_id('11').set_color('#8abe48')

# Subset alphas
v.get_patch_by_id('10').set_alpha(0.7)
v.get_patch_by_id('01').set_alpha(0.7)
v.get_patch_by_id('11').set_alpha(0.7)

#Border styles
c = venn2_circles(subsets=s, linestyle='solid')
for text in v.set_labels:
    text.set_fontsize(20)
for x in range(len(v.subset_labels)):
    if v.subset_labels[x] is not None:
        v.subset_labels[x].set_fontsize(20)

plt.title("Inclusion Protein Comparison 1% FDR")
pp_ven.savefig()

plt.show()
plt.close()
pp_ven.close()

