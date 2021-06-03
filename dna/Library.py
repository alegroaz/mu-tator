import re

import numpy as np
import pandas as pd

from dna.DNASeq import DNASeq

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

    def del_mutation(self, mut):
        extracted = re.search("(\w)(\d+)(\w)", mut)
        ori = extracted.group(1)
        pos = extracted.group(2)
        new = extracted.group(3)

        if mut in self.mutations:
            self.mutations.remove(mut)
            self.mutations.sort()
        else:
            print("The mutation you inserted is not there!")

    def mutation_quality_control(self, mut, selector):
        if selector == 1:
            if not re.match("\w\d+\w", mut):
                print("Wrong format!")
                return 1

            extracted = re.search("(\w)(\d+)(\w)", mut)
            original = extracted.group(1)
            position = int(extracted.group(2))
            mutated = extracted.group(3)

            if str(position) + original + mutated in self.mutations:
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

            parental = self.get_parental().get_aa_sequence()

            if position > len(parental):
                print("Position exceeds sequence length!")
                return 1

            if parental[position - 1] != original:
                print("The aa in position", position, "is not a", original)
                return 1

        else:
            if not re.match("\d+\w", mut):
                print("Wrong format!")
                return 1

            extracted = re.search("(\d+)(\w)", mut)
            position = int(extracted.group(1))
            mutated = extracted.group(2)

            if not re.match(r"[GPAVLIMCFYWHKRQNEDST]", mutated):
                print("The mutated aa you selected is not an aa!")
                return 1

            if not re.match(r"\w" + str(position) + mutated, mut) in self.mutations:
                print("Mutation non existent!")
                return 1

        # print(original, "in position", position, "mutated to a", mutated)

    def mod_mutations(self, selector):
        if selector == 1:
            print("\nInput desired mutations in standard format (e.g. K47A). Type STOP to stop:")

            count = 1
            while True:
                mut = input(f"Mutation {count}: ").upper()
                if mut == "STOP" or mut == "stop":
                    break

                flag = self.mutation_quality_control(mut, selector)
                if flag == 1:
                    continue
                self.add_mutation(mut)
                count += 1

        else:
            print("\nInput the mutations you want to delete in format [position][aminoacid] (e.g. 47A). Type STOP to stop:")

            count = 1
            while True:
                mut = input(f"Mutation to delete {count}: ").upper()
                if mut == "STOP" or mut == "stop":
                    break

                flag = self.mutation_quality_control(mut, selector)
                if flag == 1:
                    continue
                self.del_mutation(mut)
                count += 1

    def print_parental(self):
        print("\033[1;31mParental DNA:\033[0m")
        nt_seq = self.parental.get_nt_sequence()
        aa_seq = self.parental.get_aa_sequence()
        # print(list(range(0, len(nt_seq), 3)))
        # print((list(nt_seq[i])) for i in list(range(0, len(nt_seq), 3)))
        df = pd.DataFrame(columns=range(1, len(aa_seq) + 1),
                          index=["nt", "aa"])

        df.loc['nt'] = [(str(nt_seq[i:i + 3])) for i in range(0, len(nt_seq), 3)]
        df.loc['aa'] = [(str(aa_seq[i])) for i in range(0, len(aa_seq))]

        print(df.to_string())

    def print_lib(self):
        print("\033[1;31mMutations Table:\033[0m")
        pd.set_option('expand_frame_repr', False)
        nt_seq = self.parental.get_nt_sequence()
        aa_seq = self.parental.get_aa_sequence()
        #print(list(range(0, len(nt_seq), 3)))
        #print((list(nt_seq[i])) for i in list(range(0, len(nt_seq), 3)))
        df = pd.DataFrame(columns=range(1, len(aa_seq)+1),
                          index=["nt", "aa"])


        df.loc['nt'] = [(str(nt_seq[i:i + 3])) for i in range(0, len(nt_seq), 3)]
        df.loc['aa'] = [(str(aa_seq[i])) for i in range(0, len(aa_seq))]


        mutations_table = np.full((1, len(self.parental.get_aa_sequence())), "")

        for mut in self.mutations:
            position = int(re.findall("(\d+)", mut)[0]) - 1
            target = re.findall("\d+\w(\w)", mut)[0]
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
        if len(self.mutations) != 0:
            df.loc[''] = ["-" for i in range(0, len(aa_seq))]
            for i in mutations_table.index:
                renamer[i] = "mut%s" % str(i+1)

            mutations_table = mutations_table.rename(renamer)
            df = df.append(mutations_table)
        #df = df.append(pd.DataFrame(mutations_table, columns=range(1, len(aa_seq)+1)))
        print(df.to_string())
