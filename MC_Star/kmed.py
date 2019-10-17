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
