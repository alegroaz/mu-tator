import re
from texttable import Texttable
from graphviz import Digraph

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

import tkinter as tk
from tkinter import filedialog

import pandas as pd
import numpy as np
#######################################################################################################################


class DNASeq:
    def __init__(self, nt_seq):
        self.nt_sequence = nt_seq
        self.aa_sequence = self.nt_sequence.translate()

    def print_seq(self):
        table = Texttable(100000)
        table.set_cols_align(["c" for x in range(0, len(self.aa_sequence))])
        table.add_rows([
            [x for x in range(1, len(self.aa_sequence) + 1)],
            [(self.nt_sequence[i:i + 3]) for i in range(0, len(self.nt_sequence), 3)],
            [x for x in self.aa_sequence]
        ])
        print(table.draw())

    def get_aa_sequence(self):
        return self.aa_sequence

    def get_nt_sequence(self):
        return self.nt_sequence


#######################################################################################################################


class Library:
    def __init__(self, sequence):
        self.parental = DNASeq(sequence)
        self.mutations = []
        #mutations_table = None
        #self.mutations_df = None

    def get_parental(self):
        return self.parental

    def add_mutation(self, mut):
        extracted = re.search("(\w)(\d+)(\w)", mut)
        ori = extracted.group(1)
        pos = extracted.group(2)
        new = extracted.group(3)

        self.mutations.append(pos+ori+new)
        self.mutations.sort()

    def add_mutations(self):
        print("Input desired mutations in standard format (e.g. K47A). Type STOP to stop:")

        count = 1
        while True:
            mut = input(f"Mutation {count}: ").upper()
            if mut == "STOP" or mut == "stop":
                break

            flag = mutation_quality_control(mut, library)
            if flag == 1:
                continue
            library.add_mutation(mut)
            count += 1

    def print_lib(self):
        nt_seq = self.parental.get_nt_sequence()
        aa_seq = self.parental.get_aa_sequence()
        #print(list(range(0, len(nt_seq), 3)))
        #print((list(nt_seq[i])) for i in list(range(0, len(nt_seq), 3)))
        df = pd.DataFrame(columns=range(1, len(aa_seq)+1),
                          index=["nt", "aa"])


        df.loc['nt'] = [(str(nt_seq[i:i + 3])) for i in range(0, len(nt_seq), 3)]
        df.loc['aa'] = [(str(aa_seq[i])) for i in range(0, len(aa_seq))]
        df.loc[''] = ["-" for i in range(0, len(aa_seq))]

        mutations_table = np.full((1, len(self.parental.get_aa_sequence())), "")

        for mut in self.mutations:
            position = int(re.findall("(\d+)", mut)[0]) - 1
            target = re.findall("\d\w(\w)", mut)[0]
            height, width = mutations_table.shape

            row = 0
            while True:
                if mutations_table[row][position] == "":
                    mutations_table[row][position] = target

                    break
                else:
                    row = row + 1
                    if row > height-1:
                        newline = np.full((1, len(self.parental.get_aa_sequence())), "")
                        mutations_table = np.vstack((mutations_table, newline))
        mutations_table = pd.DataFrame(mutations_table, columns=range(1, len(aa_seq) + 1))

        renamer = dict()
        for i in mutations_table.index:
            renamer[i] = "mut%s" % str(i+1)

        mutations_table = mutations_table.rename(renamer)
        df = df.append(mutations_table)
        #df = df.append(pd.DataFrame(mutations_table, columns=range(1, len(aa_seq)+1)))
        print(df.to_string())
        df.style.set_properties(**{'background-color': 'black',
                                   'color': 'green'})


#######################################################################################################################


def sequence_quality_control(seq):
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
    return 0

def mutation_quality_control(mut, library):
    if not re.match("\w\d+\w", mut):
        print("Wrong format!")
        return 1

    extracted = re.search("(\w)(\d+)(\w)", mut)
    original = extracted.group(1)
    position = int(extracted.group(2))
    mutated = extracted.group(3)

    if str(position) + original + mutated in library.mutations:
        print("Mutation already submitted!")
        return 1

    if original == mutated:
        print("You are trying to mutate an aa into the same aa!")
        return 1

    if not re.match(r"[GPAVLIMCFYWHKRQNEDST]", original):
        print("The aa you want to mutate is not an aa!")
        return 1

    if not re.match(r"[GPAVLIMCFYWHKRQNEDST]", mutated):
        print("The mutated aa you selected is not an aa!")
        return 1

    parental = library.get_parental().get_aa_sequence()

    if position > len(parental):
        print("Position exceeds sequence length!")
        return 1

    if parental[position-1] != original:
        print("The aa in position", position, "is not a", original)
        return 1

    #print(original, "in position", position, "mutated to a", mutated)

def generate_codon_table():
    bases = "tcag"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return dict(zip(codons, amino_acids))

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
    map = dict()
    for num, feat in enumerate(selected_record.features):
        if feat.type == "C.L.":
            print(count, feat.qualifiers["label"][0])
            map[count] = feat
            count = count + 1
    sel = int(input("Select feature to synthesize: "))
    return map[sel]

def print_selection(record, feature):
    print("Feature",
          "\033[1;31m" + feature.qualifiers["label"][0] + "\033[0;0m",
          "selected from plasmid",
          "\033[1;31m" + record.name + "\033[0;0m\n")


if __name__ == '__main__':

    filename = "C:/Users/alegr/Downloads/2964-pcaggs_x2-cpmapplegrab2-ebfp2-caax.gb"
    if "filename" not in locals():
        filename = browse_seq()

    #selected_record = select_record(filename)
    #selected_feature = select_feature(select_record)


    # For testing purposes these
    selected_record = list(SeqIO.parse(filename, "genbank"))[0]
    selected_feature = selected_record.features[6]
    print_selection(selected_record, selected_feature)
    #nt_sequence = "atggagacgactgtgaggtatgaacaggggtcagagctcactaaaacttcgagctctccaacagcagatgagcccacgataaagattgatgatggtcgtgatgagggtaatgaacaagacagctgttccaataccattaggagaaaaatttccccgtttgtgatgtcatttggattcagagtatttggagttgtgcttatcattgtagacatcatagtggtgattgtggatctggccatcagtgagaagaaaagaggcattagagagattcttgaaggtgtttccctggctatagcactcttcttccttgttgatgttctcatgagagtgtttgttgaaggcttcaagaactatttccggtccaaactggttactttggatgcagtcatagtagtgggcactctgctaattaatatgacctactccttctctgaccttcaaacccccgtgtgttccgaggttatgtaccccgaggacggcgccctgaagagcgagatcaagaaggggctgaggctgaaggacggcggccactacgccgccgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacatcgtcgacatcaagttggacatcgtgtcccacaacgaggactacaccatcgtggaacagtgcgaacgcgccgaggcccgccactccaccggcggcatggacgagctgtacaagggaggtacaggcgggagtctggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgaggcctttcagaccgctaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacattaagcacccagccgacatccccgactacttcaagctgtccttccccgagggcttcaggtgggagcgcgtgatgaacttcgaggacggcggcattattcacgttaaccaggactcctccctgcaggacggctacttcatctacaaggtgaagctgcgcggcaccaacttcccccccgacggccccgtaatgcagaagaagaccatgggctgggagtccgtggtcacagatcagatgccgcacatggttactcttcttcgagttctgaaaattgttatcttaataagaatatttcgcctggcttcacagaagaaacaacttgaagtggtaaccggatccggagccaccaacttcagcctgctgaagcaggcaggcgacgtggaagagaaccctggccctgtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgaggggcgagggcgagggcgatgccaccaacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgagccacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcacctacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcgtcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacttcaacagccacaacatctatatcatggccgtcaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacgtggaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacagccactacctgagcacccagtccgtgctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttccgcaccgccgccgggatcactctcggcatggacgagctgtacaagcttaagaagatgagcaaggacggcaagaagaaaaagaagaagtccaagactaagtgcgtgatcatgtaatga"

    #print(selected_record.seq)
    nt_sequence = selected_feature.location.extract(selected_record.seq)

    library = Library(nt_sequence)
    library.add_mutations()


    #print(library.mutations)

    library.print_lib()

    '''# CODE FOR GRAPH
    s = Digraph('structs', filename='graph.pdf',
                node_attr={'shape': 'record'})
    s.node('target', '<f0> %s' % library.parental.nt_sequence)
    #s.view()'''


