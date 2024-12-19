import math
from k562.uniprot import file_converter_helper

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr
from matplotlib.backends.backend_pdf import PdfPages


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


pia_incl_dict = {x[0]:x[-1] for x in pia_incl_data_with_qvalues}
pia_pars_dict = {x[0]:x[-1] for x in pia_pars_data_with_qvalues}


inhouse_data_incl_with_pia = []
non_matching_inclusion = []
for inc in inhouse_incl_data_with_qvalues:
    if "DECOY_" not in inc[0]:
        try:
            pia_match = pia_incl_dict[inc[0]]
            inhouse_data_incl_with_pia.append([inc[0], inc[-1], pia_match, pia_match-inc[-1]])
        except KeyError:
            non_matching_inclusion.append(inc[0])



pdf_file_inclusion = 'k562/uniprot/plots2/scatter_inclusion_final_full.pdf'
png_file_inclusion = 'k562/uniprot/plots2/scatter_inclusion_final_full.png'


pp_scatter_inclusion = PdfPages(pdf_file_inclusion)

x=np.array([x[1] for x in inhouse_data_incl_with_pia])
y=np.array([x[2] for x in inhouse_data_incl_with_pia])
plt.scatter(x,y, alpha=.25)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)
m,b = np.polyfit(x,y,1)
plt.plot(x, x, c="black")
plt.xlim([0,0.2])
plt.ylim([0,0.2])
plt.xticks(np.arange(0, 0.25, 0.05))
plt.yticks(np.arange(0, 0.25, 0.05))
plt.xlabel("PyProteinInference Inclusion FDR")
plt.ylabel("PIA Report All FDR")
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
# set the axis line width in pixels
for axis in ['left', 'bottom', "right", "top"]:
  plt.subplot().spines[axis].set_linewidth(2)
# set the parameters for both axis: label size in font points, the line tick line
# width and length in pixels
plt.subplot().tick_params(axis='both', which='major', width=2)
corrincl = spearmanr([x[1] for x in inhouse_data_incl_with_pia],[x[2] for x in inhouse_data_incl_with_pia])[0]
corrincl = math.trunc(corrincl * 10000) / 10000
plt.title("Spearman Correlation: {:.4f}".format(corrincl))
plt.savefig(png_file_inclusion, dpi=300)

pp_scatter_inclusion.savefig()
plt.show()
plt.close()
pp_scatter_inclusion.close()



inhouse_data_pars_with_pia = []
non_matching_parsimony = []
for pars in inhouse_pars_data_with_qvalues:
    if "DECOY_" not in pars[0]:

        try:
            piap_match = pia_pars_dict[pars[0]]
            inhouse_data_pars_with_pia.append([pars[0], pars[-1], piap_match])
        except KeyError:
            non_matching_parsimony.append(pars[0])


pdf_file_parsimony = 'k562/uniprot/plots2/scatter_parsimony_final_full.pdf'
png_file_parsimony = 'k562/uniprot/plots2/scatter_parsimony_final_full.png'

pp_scatter_parsimony = PdfPages(pdf_file_parsimony)

x=np.array([x[1] for x in inhouse_data_pars_with_pia])
y=np.array([x[2] for x in inhouse_data_pars_with_pia])
plt.scatter(x,y, alpha=.05)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)
m,b = np.polyfit(x,y,1)
plt.plot(x, x, c="black")
plt.xlim([0,0.2])
plt.ylim([0,0.2])
plt.xticks(np.arange(0, 0.25, 0.05))
plt.yticks(np.arange(0, 0.25, 0.05))
plt.xlabel("PyProteinInference Parsimony FDR")
plt.ylabel("PIA Occams Razor FDR")
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
# set the axis line width in pixels
for axis in ['left', 'bottom', "right", "top"]:
  plt.subplot().spines[axis].set_linewidth(2)
# set the parameters for both axis: label size in font points, the line tick line
# width and length in pixels
plt.subplot().tick_params(axis='both', which='major', width=2)
corrpars = spearmanr([x[1] for x in inhouse_data_pars_with_pia],[x[2] for x in inhouse_data_pars_with_pia])[0]
corrpars = math.trunc(corrpars * 10000) / 10000
plt.title("Spearman Correlation: {:.4f}".format(corrpars))
plt.savefig(png_file_parsimony, dpi=300)
pp_scatter_parsimony.savefig()
plt.show()
plt.close()
pp_scatter_parsimony.close()


percpi_dict = {x[0]:x[-1] for x in perc_pi_data}

inhouse_data_exc_with_percpi = []
non_matching_exclusion = []
for exc in exc_data_with_qvalues:
    if "DECOY_" not in exc[0]:

        try:
            if exc[-1]<=0.20:
                percpi_match = percpi_dict[exc[0]]
                inhouse_data_exc_with_percpi.append([exc[0], exc[-1], percpi_match])
        except KeyError:
            non_matching_exclusion.append(exc[0])


pdf_file_exclusion = 'k562/uniprot/plots2/scatter_exclusion_final_full.pdf'
png_file_exclusion = 'k562/uniprot/plots2/scatter_exclusion_final_full.png'
pp_scatter_exclusion = PdfPages(pdf_file_exclusion)

x=np.array([x[1] for x in inhouse_data_exc_with_percpi])
y=np.array([x[2] for x in inhouse_data_exc_with_percpi])
plt.scatter(x,y, alpha=.05)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)
m,b = np.polyfit(x,y,1)
plt.plot(x, x, c="black")
plt.xlim([0,0.2])
plt.ylim([0,0.2])
plt.xticks(np.arange(0, 0.25, 0.05))
plt.yticks(np.arange(0, 0.25, 0.05))
plt.xlabel("PyProteinInference Exclusion FDR")
plt.ylabel("Percolator Protein FDR")
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(left=0.2)
# set the axis line width in pixels
for axis in ['left', 'bottom', "right", "top"]:
  plt.subplot().spines[axis].set_linewidth(2)
# set the parameters for both axis: label size in font points, the line tick line
# width and length in pixels
plt.subplot().tick_params(axis='both', which='major', width=2)
correxcl = spearmanr([x[1] for x in inhouse_data_exc_with_percpi],[x[2] for x in inhouse_data_exc_with_percpi])[0]
correxcl = math.trunc(correxcl * 10000) / 10000
plt.title("Spearman Correlation: {:.4f}".format(correxcl))
plt.savefig(png_file_exclusion, dpi=300)
pp_scatter_exclusion.savefig()
plt.show()
plt.close()
pp_scatter_exclusion.close()