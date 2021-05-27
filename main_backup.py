import re
from texttable import Texttable
from graphviz import Digraph
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
#######################################################################################################################


class DNASeq:
    def __init__(self, nt_seq, codons):
        self.nt_sequence = Seq(nt_seq)
        #self.aa_sequence = self.translate(self.nt_sequence, codons)
        self.aa_sequence = self.nt_sequence.translate(cds=True)

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
    def __init__(self, parental):
        self.parental = parental
        self.mutations = []
        self.mutations_table = None

    def get_parental(self):
        return self.parental

    def add_mutation(self, mut):
        extracted = re.search("(\w)(\d+)(\w)", mut)
        ori = extracted.group(1)
        pos = extracted.group(2)
        new = extracted.group(3)

        self.mutations.append(pos+ori+new)
        self.mutations.sort()

    def print_lib(self):
        if self.mutations_table is None:
            nt_sequence = self.parental.get_nt_sequence()
            aa_sequence = self.parental.get_aa_sequence()

            self.mutations_table = Texttable(100000)
            self.mutations_table.set_cols_align(["c" for x in range(0, len(aa_sequence))])
            self.mutations_table.add_rows([
                [x for x in range(1, len(aa_sequence) + 1)],
                [(nt_sequence[i:i + 3]) for i in range(0, len(nt_sequence), 3)],
                [x for x in aa_sequence]
            ])

            depth = 0
            prev = None
            new_rows = []
            #print("LENGTH:", len(self.mutations))
            for mut in self.mutations:
                position = int(re.findall("(\d+)", mut)[0]) - 1
                print(type(position))
                target = re.findall("\d\w(\w)", mut)[0]
                #print(target)

                current = position
                if prev is None:
                    new_row_aa = ["" for i in range(0, len(aa_sequence))]
                    new_row_aa[position] = target
                    # new_row_nt = []
                    new_rows.append(new_row_aa)
                    prev = position
                    self.mutations_table.add_row(new_rows[depth])

                elif position == prev:
                    depth = depth + 1
                    new_row_aa = ["" for i in range(0, len(aa_sequence))]
                    new_row_aa[position] = target
                    new_rows.append(new_row_aa)
                    prev = position
                    self.mutations_table.add_row(new_rows[depth])


                '''if prev == current:
                    self.mutations_table.add_rows([
                        [(nt_sequence[i:i + 3]) for i in range(0, len(nt_sequence), 3)],
                        [x for x in aa_sequence]
                    ])'''


            row = []
            #for mut in self.mutations:


            print(self.mutations_table.draw())





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


if __name__ == '__main__':
    nt_sequence = "atggagacgactgtgaggtatgaacaggggtcagagctcactaaaacttcgagctctccaacagcagatgagcccacgataaagattgatgatggtcgtgatgagggtaatgaacaagacagctgttccaataccattaggagaaaaatttccccgtttgtgatgtcatttggattcagagtatttggagttgtgcttatcattgtagacatcatagtggtgattgtggatctggccatcagtgagaagaaaagaggcattagagagattcttgaaggtgtttccctggctatagcactcttcttccttgttgatgttctcatgagagtgtttgttgaaggcttcaagaactatttccggtccaaactggttactttggatgcagtcatagtagtgggcactctgctaattaatatgacctactccttctctgaccttcaaacccccgtgtgttccgaggttatgtaccccgaggacggcgccctgaagagcgagatcaagaaggggctgaggctgaaggacggcggccactacgccgccgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacatcgtcgacatcaagttggacatcgtgtcccacaacgaggactacaccatcgtggaacagtgcgaacgcgccgaggcccgccactccaccggcggcatggacgagctgtacaagggaggtacaggcgggagtctggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgaggcctttcagaccgctaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacattaagcacccagccgacatccccgactacttcaagctgtccttccccgagggcttcaggtgggagcgcgtgatgaacttcgaggacggcggcattattcacgttaaccaggactcctccctgcaggacggctacttcatctacaaggtgaagctgcgcggcaccaacttcccccccgacggccccgtaatgcagaagaagaccatgggctgggagtccgtggtcacagatcagatgccgcacatggttactcttcttcgagttctgaaaattgttatcttaataagaatatttcgcctggcttcacagaagaaacaacttgaagtggtaaccggatccggagccaccaacttcagcctgctgaagcaggcaggcgacgtggaagagaaccctggccctgtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgaggggcgagggcgagggcgatgccaccaacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgagccacggcgtgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcacctacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcgtcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacttcaacagccacaacatctatatcatggccgtcaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacgtggaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacagccactacctgagcacccagtccgtgctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttccgcaccgccgccgggatcactctcggcatggacgagctgtacaagcttaagaagatgagcaaggacggcaagaagaaaaagaagaagtccaagactaagtgcgtgatcatgtaatga"
    #nt_sequence = "atggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgcatggcttgctag"

    for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
        print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    flag = 0
    while flag == 1:
        print("Input sequence to mutate:")
        nt_sequence = input()
        flag = sequence_quality_control(nt_sequence)
    library = Library(DNASeq(nt_sequence, codon_table))
    library.parental.print_seq()
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

    print(library.mutations)

    library.print_lib()

    # CODE FOR GRAPH
    s = Digraph('structs', filename='graph.pdf',
                node_attr={'shape': 'record'})
    s.node('target', '<f0> %s' % library.parental.nt_sequence)
    #s.view()


