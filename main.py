from Bio import SeqIO

import tkinter as tk
from tkinter import filedialog

from dna.Library import Library
from collections import Counter

#######################################################################################################################

'''def sequence_quality_control(seq):
    if seq == "":
        return 0
    if not re.match(r"^[atcgATCG]*$", seq):
        print("Invalid bases!")
        return 1
    if not re.match(r"ATG|atg", seq):
        print("Does not start with a START codon!")
        return 1
    if not re.match(r".*(TAA$|taa$|TGA$|tga$|TAG$|tag$)", seq):
        print("Does not end with a STOP codon!")
        return 1
    if len(seq) % 3 != 0:
        print("The sequence length is not a multiple of 3!")
        return 1
    return 0'''

def browse_seq():
    tk.Tk().withdraw()
    filename = tk.filedialog.askopenfilename()
    return filename

def select_record(filename):
    records = list(SeqIO.parse(filename, "genbank"))
    if len(records) == 1:
        print(records[0].name, "selected. No other records in .gb file.\n")
        sel = records[0]
    else:
        for count, seq_record in enumerate(records):
            print(count + 1, ":", seq_record.name)
            sel = input("Select record: ")
            sel = records[sel-1]
    return sel

def select_feature(selected_record):
    count = 1
    lib_feats = dict()
    for num, feat in enumerate(selected_record.features):
        if feat.type == "C.L.":
            print(count, feat.qualifiers["label"][0])
            lib_feats[count] = feat
            count = count + 1

    while True:
        sel = input(r"Select feature to synthesize. Type 'EXIT' to close: ")
        if sel == "EXIT" or sel == "exit":
            quit()
        elif int(sel) not in lib_feats.keys():
            print("Wrong selection")
            continue
        else:
            print_selection(selected_record, lib_feats[int(sel)])
            return lib_feats[int(sel)]

def print_selection(record, feature):
    print("\nFeature",
          "\033[1;31m" + feature.qualifiers["label"][0] + "\033[0;0m",
          "selected from plasmid",
          "\033[1;31m" + record.name + "\033[0;0m\n")


if __name__ == '__main__':

    # Set filename for testing. If not set, a browser window opens
    filename = r"C:\Users\alegr\Downloads\2801-pcaggs-jedi2p-cyofp1.gb"
    #filename = "C:/Users/alegr/Downloads/2964-pcaggs_x2-cpmapplegrab2-ebfp2-caax.gb"

    if "filename" not in locals():
        filename = browse_seq()

    # Lets the user select the record (plasmid) if more than one is in the .gb file.
    selected_record = select_record(filename)

    # Lets the user select the feature from the record. The feature selected (needs to be of type "C.L.")
    # dictates the frame of reference for the numbering of the residues.
    selected_feature = select_feature(selected_record)
    nt_sequence = selected_feature.location.extract(selected_record.seq)

    # Creates the library object, and prints the parental sequence.
    library = Library(nt_sequence)
    library.print_parental()

    # Allows the insertion of mutations. When done, prints the table with the mutations.
    library.mod_mutations(1)
    library.print_lib()



###### TESTING STUFF #####################################

    '''c1 = input("\nDoes this look correct? y/n")
    while True:
        if c1 == "y" or c1 == "Y":
            print("Oro, go ahead")
            break
        elif c1 == "n" or c1 == "N":
            print("What to do?\n"
                  "0: It was actually correct\n"
                  "1: Add a mutation\n"
                  "2: Remove a mutation")
            c2 = input("Choice: ")
            if c2 == "0":
                c1 = "y"
                continue
            elif c2 == "1":
                library.mod_mutations(1)
            elif c2 == "2":
                library.mod_mutations(-1)'''



    '''# CODE FOR GRAPH
    s = Digraph('structs', filename='graph.pdf',
                node_attr={'shape': 'record'})
    s.node('target', '<f0> %s' % library.parental.nt_sequence)
    #s.view()'''


