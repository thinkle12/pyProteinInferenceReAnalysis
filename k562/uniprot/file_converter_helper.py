import os
import csv
import matplotlib.pyplot as plt
import re


# Helpers to generate benchmark data

def convert_pyproteininference_ids_for_pia(protein_id):
    # Check if the protein is a decoy
    if protein_id.startswith("DECOY_"):
        return f"{protein_id}"
    else:
        # Extract the type (sp| or tr|), accession, and isoform if present
        match = re.match(r"(sp|tr)\|(\w+)(-\d+)?\|", protein_id)
        if match:
            accession = match.group(2)  # Q8NGE2 or K7EM90
            isoform = match.group(3) if match.group(3) else ""  # -2 or ""
            return f"{accession}{isoform}"
        else:
            return protein_id  # Return unmodified if it doesn't match the pattern


def convert_perc_pi(target_filename, decoy_filename):
    complete_target = []
    with open(target_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_target.append(row)
    target_table = [[x[0], float(x[3]), float(x[2]), x[-1]] for x in complete_target[1:]]


    complete_decoy = []
    with open(decoy_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_decoy.append(row)
    decoy_table = [[x[0], float(x[3]), float(x[2]), x[-1]] for x in complete_decoy[1:]]

    table = target_table+decoy_table

    table = sorted(table, key=lambda x: (x[1],x[2]))

    return table

def convert_proteinsieve(filename):
    # Convert to [protein, score, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    table = [[x[1], float(x[-1]), float(x[-3])] for x in complete_file[1:]]

    return table



def convert_fido(target_filename,decoy_filename):
    complete_target = []
    with open(target_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_target.append(row)
    target_table = [[x[0], float(x[3]), float(x[2]), x[-1]] for x in complete_target[1:]]


    complete_decoy = []
    with open(decoy_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_decoy.append(row)
    decoy_table = [[x[0], float(x[3]), float(x[2]), x[-1]] for x in complete_decoy[1:]]

    table = target_table+decoy_table

    table = sorted(table, key=lambda x: (x[1],x[-1]))

    return table


def convert_internal_pi(filename):
    # Convert to [protein, score, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    table = [[x[0], float(x[1]), float(x[2]), " ".join(x[6:])] for x in complete_file[1:]]

    return table

def convert_pia_old(filename, ftype="parsimony"):
    # Convert to [protein, score, clusterID, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar=']')
        for row in spamreader:
            complete_file.append(row)

    if ftype=="inclusion":
        table = []
        for x in complete_file[1:]:

            string_proteins = x[0]
            cleaned_string = string_proteins.strip('"[]')

            # Split by comma and strip any extra spaces
            proteins = [item.strip() for item in cleaned_string.split(',')]

            for prot in proteins:
                table.append([prot, float(x[1]), float(x[6]), float(x[-1])])

    elif ftype=="parsimony":
        table = []
        for x in complete_file[1:]:
            proteins = x[0]
            cleaned_string = proteins.strip('"[]')

            # Split by comma and strip any extra spaces
            parsed_list = [item.strip() for item in cleaned_string.split(',')]

            lead_protein = parsed_list[0]
            table.append([lead_protein, float(x[1]), float(x[6]), float(x[-1])])

    else:
        pass

    return table

def convert_pia_new(filename, ftype="parsimony", pyproteininference_proteins=None):
    # Convert to [protein, score, clusterID, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar=']')
        for row in spamreader:
            complete_file.append(row)

    if ftype=="inclusion":
        table = []
        for x in complete_file[1:]:

            string_proteins = x[0]
            cleaned_string = string_proteins.strip('"[]')

            # Split by comma and strip any extra spaces
            proteins = [item.strip() for item in cleaned_string.split(',')]

            for prot in proteins:
                table.append([prot, float(x[1]), float(x[6]), float(x[-1])])

    elif ftype=="parsimony":
        table = []
        for x in complete_file[1:]:
            proteins = x[0]
            cleaned_string = proteins.strip('"[]')

            # Split by comma and strip any extra spaces
            parsed_list = [item.strip() for item in cleaned_string.split(',')]

            matching_proteins = [x for x in parsed_list if x in pyproteininference_proteins]
            if matching_proteins:
                lead_protein = matching_proteins[0]
            else:
                lead_protein = proteins[0]

            table.append([lead_protein, float(x[1]), float(x[6]), float(x[-1])])

    else:
        pass

    return table


def convert_pia(filename, ftype="parsimony", pyproteininference_proteins=None):
    # Convert to [protein, score, clusterID, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar=']')
        for row in spamreader:
            complete_file.append(row)

    if ftype=="inclusion":
        table = []
        for x in complete_file[1:]:
            row = x[0].split(",")
            proteins = row[-2].split(";")
            for prot in proteins:
                table.append([prot, float(row[-1])])

    elif ftype=="parsimony":
        # TODO use pyproteininference_proteins here...
        table = []
        for x in complete_file[1:]:
            row = x[0].split(",")
            proteins = row[-2].split(";")
            matching_proteins = [x for x in proteins if x in pyproteininference_proteins]
            if matching_proteins:
                lead_protein = matching_proteins[0]
            else:
                lead_protein = proteins[0]
            table.append([lead_protein, float(row[-1])])

    else:
        pass

    return table


def convert_protein_prophet(filename, only_leads=True):
    # Convert to [protein, score, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    if only_leads:
        leads = []
        for rows in complete_file:
            try:
                float(rows[0])
                leads.append(rows)
            except ValueError:
                if "a" in rows[0]:
                    leads.append(rows)
                else:
                    pass
    else:
        leads = complete_file[1:]

    leads = [x for x in leads if float(x[2])>0]

    table = [[x[1], float(x[6]), float(x[2]), x[-1]] for x in leads]

    table = sorted(table, key=lambda x: (x[2],x[1]), reverse=True)

    return table

def calc_qvalues(proper_list, decoy_symbol="DECOY_"):
    # Now pick out only the lead protein identifiers
    lead_proteins = [x[0] for x in proper_list]

    # Reverse the list (best to worst) -> (worst to best)
    lead_proteins.reverse()

    fdr_list = []
    for i in range(len(lead_proteins)):
        binary_decoy_target_list = [1 if decoy_symbol in elem else 0 for elem in lead_proteins]
        total = len(lead_proteins)
        decoys = sum(binary_decoy_target_list)
        fdr = (2 * decoys) / (float(total))
        fdr_list.append(fdr)
        del lead_proteins[0]

    qvalue_list = []
    new_fdr_list = []
    for fdrs in fdr_list:
        new_fdr_list.append(fdrs)
        qvalue = min(new_fdr_list)
        qvalue_list.append(qvalue)

    qvalue_list.reverse()

    for k in range(len(proper_list)):
        proper_list[k].append(qvalue_list[k])

    return proper_list


def plot_benchmark(dfdr,efdr,pdf_file = None, title = None):
    f = plt.figure()
    plt.plot([0, 1], color = 'k')
    plt.plot([0, 0.6666666666666666], '--', color = 'k')
    plt.plot([0, 1.5], '--', color = 'k')
    plt.plot(dfdr, efdr, "-o")
    plt.xlim([0.00001, .1])
    plt.ylim([0.00001, .1])
    plt.xlabel('Decoy FDR')
    plt.ylabel('Entrapment FDR')
    plt.title(title)
    if pdf_file:
        f.savefig(pdf_file)
    plt.show()
    plt.close()


def roc_plots(proper_list, fdr0_list, pdf_file = None, title = None, true_db = None):
    import matplotlib.pyplot as plt
    main_list = []
    count_list = []
    fdrs = fdr0_list[0]
    fdrs.reverse()
    for x in range(len(proper_list)):
        if '##' not in proper_list[x][0]:
            count_list.append(proper_list[x])
        main_list.append([fdrs[x],len(count_list)])
    f = plt.figure()

    main_list = [x for x in main_list if x[0]<.2]

    plt.plot([x[0] for x in main_list], [x[1] for x in main_list], '-')

    plt.xlabel('FDR')
    plt.ylabel('Hits')
    plt.xlim([-.01, .2])
    plt.title(title)
    if pdf_file:
        f.savefig(pdf_file)
    plt.show()
    plt.close()


def convert_fido_groups(target_filename,decoy_filename):
    complete_target = []
    with open(target_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_target.append(row)
    target_table = [[x[0], int(x[1]), float(x[3]), float(x[2])] for x in complete_target[1:]]


    complete_decoy = []
    with open(decoy_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_decoy.append(row)
    decoy_table = [[x[0], int(x[1]), float(x[3]), float(x[2])] for x in complete_decoy[1:]]

    table = target_table+decoy_table

    table = sorted(table, key=lambda x: (x[2],x[-1]))

    group_tracker = []
    new_table = []
    for items in table:
        if items[1] not in group_tracker:
            group_tracker.append(items[1])
            new_table.append(items)
        else:
            pass


    return new_table

def convert_fido_and_apply_groups(target_filename,decoy_filename):
    complete_target = []
    with open(target_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_target.append(row)
    target_table = [[x[0], int(x[1]), float(x[3]), float(x[2]), x[4].split(" ")] for x in complete_target[1:]]


    complete_decoy = []
    with open(decoy_filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in spamreader:
            complete_decoy.append(row)
    decoy_table = [[x[0], int(x[1]), float(x[3]), float(x[2]), x[4].split(" ")] for x in complete_decoy[1:]]

    table = target_table+decoy_table

    # Sort by len...
    table = sorted(table, key=lambda x: len(x[-1]))

    set1 = [set(x[-1]) for x in table]

    [x.remove('') for x in set1]

    sets1 = list([frozenset(e) for e in set1])

    sets = [frozenset(e) for e in set1]
    print(str(len(sets)) + ' number of peptide sets')
    us = set()
    i = 0
    # Get all peptide sets that are not a subset...
    while sets:
        i = i + 1
        e = sets.pop()
        if any(e.issubset(s) for s in sets) or any(e.issubset(s) for s in us):
            continue
        else:
            us.add(e)
        if i % 10000 == 0:
            print("Parsed {} Peptide Sets".format(i))

    dictionary = {"us" : us, "sets" : sets1, "table":table}


    return dictionary


def convert_inhouse_and_apply_groups(filename):
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    table = [[x[0], float(x[1]), float(x[2]), set(x[6:])] for x in complete_file[1:]]

    # Sort by len...
    table = sorted(table, key=lambda x: len(x[-1]))

    set1 = [set(x[-1]) for x in table]

    sets1 = list([frozenset(e) for e in set1])

    sets = [frozenset(e) for e in set1]
    print(str(len(sets)) + ' number of peptide sets')
    us = set()
    i = 0
    # Get all peptide sets that are not a subset...
    while sets:
        i = i + 1
        e = sets.pop()
        if any(e.issubset(s) for s in sets) or any(e.issubset(s) for s in us):
            continue
        else:
            us.add(e)
        if i % 10000 == 0:
            print("Parsed {} Peptide Sets".format(i))

    dictionary = {"us" : us, "sets" : sets1, "table":table}

    return dictionary


def convert_internal_pi_inclusion(filename):
    # Convert to [protein, score, fdr/qval] Sorted by score
    complete_file = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            complete_file.append(row)

    table = [[x[0], float(x[1]), float(x[2])] for x in complete_file[1:]]

    return table