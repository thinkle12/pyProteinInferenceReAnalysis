## TODO Write a converter from pyproteininference to proteome standard
import collections
## TODO just need two columns, ProteinId	q-value

## TODO also, protein ids can be , separated (Peptide Centric) - Replace all ; with , and no space...


## Just loop over each file, read in parse, write out. Also needs to write out as tsv

import os
import csv

in_dir = "prest/data/pyproteininference"
filenames = os.listdir(in_dir)
filenames = [x for x in filenames if ".csv" in x]

write_directory = "prest/proteoform-standard/input"

for filename in filenames:
    basename = os.path.basename(filename)
    full_filename = os.path.join(in_dir, filename)
    complete_file = []
    with open(full_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    prot_map = {}
    if complete_file[0][6] == "Other_Potential_Identifiers":
        for row in complete_file[1:]:
            prot_map[row[0]] = row[6:]

        parsimony_group_mapper = {x:",".join(prot_map[x]) for x in prot_map.keys()}


    restricted = [[x[0],x[2]] for x in complete_file]

    fixed = []
    for prot in restricted:
        if ";" in prot[0]:
            prot[0] = prot[0].replace(";", ",")
        if "DECOY_" in prot[0]:
            pass
        else:
            fixed.append(prot)

    # Define the output TSV file path
    filename_only = basename.split(".")[0]
    tsv_output_file = os.path.join(write_directory, f"{filename_only}.tsv")

    filtered_data_fixed = [fixed[0]]
    observed_proteins = set()
    if "parsimony" in basename:
        for row in fixed[1:]:
            if parsimony_group_mapper[row[0]]:
                row[0] = row[0]+","+parsimony_group_mapper[row[0]]
            filtered_data_fixed.append(row)

        fixed = filtered_data_fixed
    if "peptide_centric" in basename:
        print("Fixing Peptide Centric...")
        # Process each row
        i = 0
        for row in fixed[1:]:
            i = i +1
            proteins = row[0].split(",")  # Split protein identifiers by commas
            if all(protein not in observed_proteins for protein in proteins):  # Check if any protein is new
                filtered_data_fixed.append(row)  # Add the row to the filtered data
                observed_proteins.update(set(proteins))  # Add all proteins to the tracking set
            else:
                print(str(proteins) + " Already added...")

        fixed = filtered_data_fixed
    # Write the fixed list to the TSV file
    with open(tsv_output_file, 'w', newline='') as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter='\t')
        tsv_writer.writerows(fixed)
