from Bio import AlignIO
from Bio.Align import AlignInfo
from statistics import mean
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Plotting:
    def __init__(self, file):
        self.alignment = AlignIO.read(file, "clustal")

    def __per_iden(self):
        # This function will return a position base per-match to consensus in a list.
        summary_align = AlignInfo.SummaryInfo(self.alignment)
        consensus = summary_align.gap_consensus()  # consensus which allowed gap
        my_pssm = summary_align.pos_specific_score_matrix(consensus)  # position specific score matrix
        ts = len(self.alignment)  # ts for short for total sequences of alignment
        idty = []  # means average identity, for storing the positional identity matching to consensus.
        consensus = list(consensus)
        for i in range(0, len(consensus)):
            # This part if to convert X in consensus to most common nucleotide. In case if want to return that.
            # if consensus[i] == "X":
            #     nude = max([(value, key) for key, value in my_pssm[i].items()])
            #     consensus[i] = nude[1]  # replace the X
            #     idty.append(nude[0] / ts)
            # else:
            score = max(my_pssm[i].values()) / ts
            idty.append(score)
        return idty

    def window_average(self, k):
        # This function calculate average identity of (i,i+k-1) window.
        # Since in python last position are ignored, window in python code is [i:i+k]
        pos_idty = self.__per_iden()
        win_avg = []
        for i in range(0, len(pos_idty) - k + 1):
            score = mean(pos_idty[i:i + k])
            win_avg.append(score)
        return win_avg

    def perfect_matches_in_window(self, k):
        # This function calculate how many perfect match in each sliding window.
        # It return a list of numbers of perfect matches in sliding window, allowing to draw histogram.
        sw = []  # shot for sliding window
        idty = self.__per_iden()
        for i in range(0, len(idty) - k + 1):  # Position from 0 to end - k
            count = 0
            for j in idty[i:i + k]:
                if j == 1:
                    count += 1
            sw.append(count)
        return sw

    def permit_mismatch(self, k):
        # This function calculate numbers of windows allowing different numbers of mismatch
        sw = np.array(self.perfect_matches_in_window(k))  # Turn list into array so it allows us to easy filter windows by score.
        total = len(sw)  # Number of windows in total
        lower, upper = (min(sw), max(sw))
        count = {"windows": [],
                 "percentage": []}  # Store information in dictionary for turning it into pandas DataFrame
        index = []
        for i in range(upper, lower - 1, -1):  # In decrement order
            index.append(k - i)  # index means how many mismatch are allowed
            # Count numbers of windows score higher then i, and append to list "count"
            sumup = len(sw[sw >= i])
            count["windows"].append(sumup)
            count["percentage"].append(sumup / total)
        return pd.DataFrame(count, index=index)

    def permit_edit(self, k):
        win_avg = np.array(self.window_average(k))
        total = len(win_avg)
        lower, upper = (min(win_avg),max(win_avg))
        step = (upper-lower)/k
        count = {"windows": [], "percentage": []}
        index = []
        for i in np.arange(lower,upper,step):
            index.append(1 - i)
            sum_up = len(win_avg[win_avg >= i])
            count["windows"].append(sum_up)
            count["percentage"].append(sum_up / total)
        return pd.DataFrame(count, index=index)

    def avgw_plt(self, k=20):
        win_avg = self.window_average(k)
        plt.scatter(range(len(win_avg)), win_avg)
        plt.xlabel('position')
        plt.ylabel('average percentage identity of {}bp window'.format(k))
        plt.show()

    def win_cout(self, k=20):
        sw = self.perfect_matches_in_window(k)
        sns.countplot(x=sw)
        plt.show()

    def heat_win(self, k=20):
        sw = self.perfect_matches_in_window(k)
        bar = np.array([sw])  # turning the list into a 2D array with only one row
        sns.heatmap(bar)
        plt.show()

    def heat_avg(self, k=20):
        avg_i = []
        # plt.plot(figsize=(10,2))
        iden = self.__per_iden()
        for i in range(0, len(iden) - k + 1):
            score = mean(iden[i:i + k])
            avg_i.append(score)
        sns.heatmap(np.array([avg_i]))
        plt.show()
