import os
from .MSA import clustalo, mafft, t_coffee
from Bio import SeqIO
from .kmed import KmedGrouping


class Clustering:
    def __init__(self, inf):
        self.inf = inf
        self.seq_dic = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))

    def __msa(self, name):
        # Name argument should be a MSA method in string format.
        # This function will build a new directory and store the alignment output in new directory.
        # Finally it cd to the new directory and return the name of the alignment file name.
        os.mkdir(name)
        exe = "self.outf = {}(self.inf, outdir='{}')".format(name, name)  # Here I can only us self otherwise in function local veriable environment it cannot define outf
        exec(exe)
        print(self.outf)
        os.chdir("./"+name)
        return


    def __sub_aln(self, k, sub_sets, name):
        # This function take sub_sets and turn it into sub-fasta and sub-alignment file
        for i in range(0, k):
            subfile = self.inf.replace(".fasta", "_set" + str(i) + ".fasta")
            with open(subfile, "w") as handle:
                for seq_id in sub_sets[i]:
                    rec = self.seq_dic[seq_id]
                    SeqIO.write(rec, handle, "fasta")
            exec("{}(subfile)".format(name))

    def clo(self, k=3):
        self.__msa("clustalo")
        seq_sets = KmedGrouping(self.outf).get_subfa(k)
        self.__sub_aln(k, seq_sets, "clustalo")

    def mft(self, k=3):
        outf = self.__msa("mafft")
        seq_sets = KmedGrouping(self.outf).get_subfa(k)
        self.__sub_aln(k, seq_sets, "mafft")

    def tcf(self, k=3):
        outf = self.__msa("t_coffee")
        seq_sets = KmedGrouping(outf).get_subfa(k)
        self.__sub_aln(k, seq_sets, "t_coffee")
