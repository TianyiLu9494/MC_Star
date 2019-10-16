from .plotting import Plotting
from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import os


class AlnPlotting:
    def __init__(self, directory):
        self.alignment = {}
        self.fatsa = {}
        for alnfile in os.listdir(directory + "aln/"):
            self.alignment[alnfile] = Plotting(directory + "aln/" + alnfile)
        for seqfile in os.listdir(directory + "seq/"):
            print(directory + "./aln/" + seqfile)
            self.fatsa[seqfile] = SeqIO.parse(directory + "aln/" + seqfile,"fasta")
        # for file in os.listdir(dir):
        #     if file[-3::] == "aln":
        #         self.alignment[file] = Plotting(dir+file)
        #     elif file[-5::] == "fasta":
        #         self.fasta = SeqIO.parse(dir+file)


    def permit_mismatch(self, window_size):
        fig, ax = plt.subplots(1, 2)
        for filename, instance in self.alignment.items():
            lab = filename[0:-13]  # this name is only fit for clostalo since "_clostal0.aln" have 13 digit
            data0 = instance.permit_edit(window_size)
            data1 = instance.permit_mismatch(window_size)
            ax[0].plot(data0["percentage"], label=lab)
            ax[0].legend()
            ax[1].plot(data1["percentage"], label=lab)
            ax[1].legend()
        plt.show()


    def heat(self, window_size):
        fig, ax = plt.subplots(1, 2)
        matrix0 = np.array([instance.perfect_matches_in_window(window_size) for instance in self.alignment.values()])
        matrix1 = np.array([instance.window_average(window_size) for instance in self.alignment.values()])
        sns.heatmap(matrix0,ax=ax[0])
        sns.heatmap(matrix1,ax=ax[1])
        plt.show()


    def test_lenth(self):
        for rec in self.fatsa["A_exon2_set3.fasta"]:
            print(rec)
            break

# def scarter(self,window_size):
