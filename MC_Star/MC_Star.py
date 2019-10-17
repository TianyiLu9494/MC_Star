'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) TianyiLu, 16 Oct 2019 
License     : MIT 
Maintainer  : lutianyi9494@icloud.com 
Portability : POSIX

The program reads one or more input FASTA files. For each file it computes a
variety of statistics, and then prints a summary of the statistics as output.
'''

from argparse import ArgumentParser
import os
import sys
import logging
import pkg_resources
from Bio import SeqIO

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
DEFAULT_WIN_SIZE = 20
DEFAULT_VERBOSE = False
DEFAULT_K = 5
HEADER = 'FILENAME\tNUMSEQ\tTOTAL\tMIN\tAVG\tMAX'
PROGRAM_NAME = "MC_Star"

try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more FASTA files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument('-w',
                        '--winsize',
                        metavar='N',
                        type=int,
                        default=DEFAULT_WIN_SIZE,
                        help='The window size to screening primer ability in plotting (default {})'.format(
                            DEFAULT_WIN_SIZE))
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('allele_file',
                        nargs='1',
                        metavar='FASTA_FILE',
                        type=str,
                        help='Input alleles file in FASTA format')
    parser.add_argument('-k',
                        nargs='1',
                        metavar='N',
                        default=DEFAULT_K,
                        type=int,
                        help='The number of clusters to generate (default {})'.format(
                            DEFAULT_K))
    return parser.parse_args()


def mafft(inf, outdir=None):
    '''
    Run the mafft by biopython wrapper.
    :param inf: The input file.
    :param outdir: The output directory.
    :return: Output filename.
    '''
    from Bio.Align.Applications import MafftCommandline
    from io import StringIO
    from Bio import AlignIO
    mafft_cline = MafftCommandline("mafft", input=inf)
    print(mafft_cline)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    outfile = inf.replace('.fasta', '_mafft.aln')
    if not outdir:
        AlignIO.write(align, outfile, "clustal")
    else:
        AlignIO.write(align, outdir + '/' + outfile, "clustal")
    return outfile


class Clustering(object):
    def __init__(self, inf, k):
        self.inf = inf
        self.k = k
        self.seq_dic = SeqIO.to_dict(SeqIO.parse(inf, "fasta"))

    def msa(self):
        """
        Make a new directory outdir after the MSA software chosen and k(# of Clusters).
        Running MSA and store the output in the directory.
        :return: Whole set alignment name.
        """
        outdir = 'Mafft_' + str(self.k)
        os.mkdir(outdir)
        os.system("cp {file} {dir}/{file}".format(file=self.inf, dir=outdir))
        outf = mafft(self.inf, outdir)
        print(outf)
        os.chdir("./" + outdir)
        return outf

    def kmed(self):
        seq_sets = KmedGrouping(self.outf).get_subfa(self.k)
        self.__sub_aln(k, seq_sets, "mafft")

    def sub_aln(self, sub_sets):
        # This function take sub_sets and turn it into sub-fasta and sub-alignment file
        for i in range(0, self.k):
            subfile = self.inf.replace(".fasta", "_set" + str(i) + ".fasta")
            with open(subfile, "w") as handle:
                for seq_id in sub_sets[i]:
                    rec = self.seq_dic[seq_id]
                    SeqIO.write(rec, handle, "fasta")
            mafft(subfile)

class KmedGrouping:
    def __init__(self, filename):
        from Bio import AlignIO
        self.filename = filename
        self.aln = AlignIO.read(filename, "clustal")
        self.ns = len(self.aln)

    def __get_dm(self):
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        import numpy as np
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(self.aln)
        dm_array = np.zeros(shape=(self.ns, self.ns))
        for row in range(0, self.ns):
            for cln in range(0, self.ns):
                if cln > row:
                    dm_array[row, cln] = dm[cln][row]
                else:
                    dm_array[row, cln] = dm[row][cln]
        return dm_array

    def __kmed(self, k):
        from pyclustering.cluster.kmedoids import kmedoids
        import random
        initial_medoids = random.sample(range(0, self.ns), k)
        distant_matrix = self.__get_dm()
        kmedoids_instance = kmedoids(distant_matrix, initial_medoids, data_type='distance_matrix', itermax=200)
        kmedoids_instance.process()
        clusters = kmedoids_instance.get_clusters()
        return clusters

    def get_subaln(self, k):
        # extract sub-alignment
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        import os
        clusters = self.__kmed(k)
        new_dir = "./sub_aln_{}/".format(k)
        os.mkdir(new_dir)
        os.chdir(new_dir)
        for i in range(0, len(clusters)):
            file_name = self.filename.replace(".aln", "_sub_aln" + str(i))
            print(file_name)
            with open(file_name, "w") as handle:
                sub_aln = [self.aln[rec] for rec in clusters[i]]
                sub_aln = MultipleSeqAlignment(sub_aln)
                AlignIO.write(sub_aln, handle, "clustal")

    def get_subfa(self, k):
        # extract sub-sequences-set in fasta format
        id_set = []
        clusters = self.__kmed(k)
        for i in range(0, len(clusters)):
            id_cluster = [self.aln[rec].id for rec in clusters[i]]
            id_set.append(id_cluster)
        return id_set


def process_files(options):
    '''Compute and print FastaStats for each input FASTA file specified on the
    command line. If no FASTA files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    if options.allele_file:
        fasta_filename = options.allele_file
        logging.info("Processing FASTA file from %s", fasta_filename)
        try:
            fasta_file = open(fasta_filename)
        except IOError as exception:
            exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
        else:
            with fasta_file:
                instancee = Clustering(fasta_file,options.k)
                wsaln = instancee.msa()
                seq_sets = KmedGrouping(wsaln).get_subfa(options.k)
                instancee.sub_aln(seq_sets)
    else:
        logging.info("Processing FASTA file from stdin")
        # stats = FastaStats().from_file(sys.stdin, options.minlen)
        # print(stats.pretty("stdin"))


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    print(HEADER)
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
