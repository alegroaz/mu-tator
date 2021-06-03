import pandas as pd

class DNASeq:
    def __init__(self, nt_seq):
        self.nt_sequence = nt_seq
        self.aa_sequence = self.nt_sequence.translate()

    def print_seq(self):
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

    def get_aa_sequence(self):
        return self.aa_sequence

    def get_nt_sequence(self):
        return self.nt_sequence
