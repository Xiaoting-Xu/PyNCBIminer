import os
import shutil
import time
import numpy as np
from Bio import SeqIO, SeqRecord, Seq
import traceback
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

import sys
sys.path.append(".")
from functional import create_folder, Timer
from nt_calculator import nt_Calculator

class Aligner:
    def __init__(self, in_file, out_path,
                 max_cluster_num=10, initial_cluster_num=3,
                 long_length_ratio=0.7, fragment_length_ratio=0.5,
                 DEBUG_MODE=False):
        """ Class Separator: do separate MSA
        ----------
        Parameters
        - in_file - the input fasta file (only one) in which sequences are to be aligend
        - out_path - destination folder of output, where log files and result will be written
        - max_cluster_num - maximum number of clusters. if exceeded, retain [min_long_record_ratio] as long sequence group
        - initial_cluster_num - the initial number of clusters.
        - min_long_record_ratio - the maximum ratio of long sequence group compared to original records
        - min_long_record_num - the minimum number of records in long_sequence_group
        - initial_PIT - the initial value of pairwise identity threshold, might decrease.
        - min_PIT - the minimum value of pairwise identity threshold.
        - step_length_PIT - the value (step length) each time pairwise identity will decrease by.
        """
        # logger
        self.log_path = os.path.join(out_path, "log")
        self.timer = Timer(os.path.join(self.log_path, "DEBUG_log.txt"))
        self.DEBUG_MODE = DEBUG_MODE
        
        # calculator
        self.nt_calculator = nt_Calculator()

        # file handles
        self.in_file = in_file
        self.out_path = out_path
        self.temp_folder = os.path.join(self.out_path, "temp")
        self.temp_long_seqs_filename = "long_sequence_group.fasta"
        self.temp_aligned_filename = "aligned_long_sequence_group.fasta"
        self.temp_unaligned_short_filename = "unaligned_short.fasta"
        self.temp_unaligned_other_filename = "unaligned_other.fasta"
        create_folder(self.temp_folder)

        # parameters on clustering
        self.MAX_CLUSTER_NUM = max_cluster_num
        self.INITIAL_CLUSTER_NUM = initial_cluster_num

        # parameters on ratio of long sequence group
        self.LONG_LENGTH_RATIO = long_length_ratio
        self.FRAGMENT_LENGTH_RATIO = fragment_length_ratio

        # store information of data and clustering result
        self.length_dct = {}
        self.cluster_result = []

        # just mention the variables used in calculation, each will be set before use
        self.cluster_num = 0
        self.record_num = 0
        self.PIT = 0
        self.compare_by_add = False

    def check_redundant_description(self):
        """ check if there are redundant description in the input file 
        -------
        Returns
        - redundant - True if there is redundancy, False if not
        """
        record_iter = SeqIO.parse(self.in_file, "fasta")
        record_description = [record.description for record in record_iter]
        set_record_description = list(set(record_description))

        # if no redundancy, then there should be no change after removing duplicates
        if len(record_description) == len(set_record_description):
            redundant = False
        else:
            redundant = True

        return redundant

    def extract_index_and_length(self):
        """ extract the index and length of each sequence
        -------
        Returns
        - length_dct - a dictionary storing the {description: length} of each sequence
        """
        length_dct = {}
        record_iter = SeqIO.parse(self.in_file, "fasta")
        for record in record_iter:
            description = record.description
            length = len(record.seq)
            length_dct[description] = length

        return length_dct

    def set_cluster_num(self, new_cluster_num):
        """
        Explicitly setting parameter "cluster_num"
        ----------
        Parameters
        - new_cluster_num - literally new value of cluster_num, must be smaller than MAX_CLUSTER_NUM
        """
        assert new_cluster_num <= self.MAX_CLUSTER_NUM, "Bad cluster num from <func> set_cluster_num"

        self.cluster_num = new_cluster_num

    def hierarchical_cluster_length(self):
        """ do hierarchical clustering on lengths
        -------
        Returns
        - clusters - list of int representing the clusters each length (sequence) belongs to
        """
        # data needed to be transformed into a 2-D numpy.array
        data = np.array(list(self.length_dct.values())).reshape((-1, 1))

        # Calculate distance matrix
        dist_mat = pdist(data)

        # Build hierarchical clustering object
        linkage_obj = linkage(dist_mat, method='single')

        # Get cluster labels based on number of clusters
        cluster_num = self.cluster_num
        clusters = fcluster(linkage_obj, t=cluster_num, criterion='maxclust')

        return clusters
    
    def hierarchical_cluster_distance(self, distance_matrix, max_cluster_num=None):
        if not max_cluster_num:
            max_cluster_num = len(distance_matrix)
        max_silhouette_score = 0
        final_cluster_result = []
        last_cluster_result = []
        decrease = 0

        for i in range(2, max_cluster_num):
            num_clusters = i
            clustering = AgglomerativeClustering(n_clusters=num_clusters, 
                                                 metric="precomputed", 
                                                 linkage="average") \
                                                .fit(distance_matrix)
                                                
            cluster_result = clustering.labels_
            score = silhouette_score(X=distance_matrix, labels=cluster_result, metric="precomputed")
            
            last_cluster_result = final_cluster_result
            if score >= max_silhouette_score:
                max_silhouette_score = score
                final_cluster_result = cluster_result
                decrease = 0
            elif decrease == 2:
                break
            else:
                decrease += 1
                
            if list(cluster_result).count(cluster_result[-1]) == 1:
                if i == 2:
                    final_cluster_result = np.ones(len(final_cluster_result))
                else:
                    final_cluster_result = last_cluster_result
                break
            
        if len(final_cluster_result) == 0:
            final_cluster_result = np.ones(1)
            
        return final_cluster_result

    def count_average_length(self):
        """ calculate average length of each cluster in self.cluster_result
        -------
        Returns
        - cluster_average_length - the length of each cluster ordered by cluster id
        """
        all_lengths = np.array(list(self.length_dct.values()))
        cluster_average_length = []

        for i in range(self.cluster_num):
            cluster_id = i + 1
            idx = np.where(self.cluster_result == cluster_id)[0]  # find idx of this cluster
            if len(idx) == 0:
                cluster_average_length.append(0)
            else:
                average_length = np.mean(all_lengths[idx])  # calculate average length
                cluster_average_length.append(average_length)

        return np.array(cluster_average_length)

    def set_PIT(self, new_PIT):
        """
        Explicitly setting parameter "PIT"
        ----------
        Parameters
        - new_cluster_num - literally new value of PIT, must be smaller than MAX_CLUSTER_NUM
        """
        assert new_PIT >= self.MIN_PIT, "Bad new_PIT from <func> set_PIT"
        self.PIT = round(new_PIT, 8)

    def align_and_write(self, longest_record_idx):
        """ do MSA among the long_sequence_group and write into disk aligned group
        ----------
        Parameters
        - longest_record_idx - the indices of long_sequence_group in length_dct
        """
        # substep 1: fetch long sequence group according to longest_record_idx
        record_list = list(SeqIO.parse(self.in_file, "fasta"))
        long_sequence_group = []
        for i in longest_record_idx:
            long_sequence_group.append(record_list[i])

        # substep 2: write into disk (temp) for multiple sequence alignment
        SeqIO.write(long_sequence_group,
                    os.path.join(self.temp_folder, self.temp_long_seqs_filename),
                    "fasta")

        # substep 3: perform multiple sequence alignment using --reorder
        mafft_in = os.path.join(self.temp_folder, self.temp_long_seqs_filename)
        mafft_out = os.path.join(self.temp_folder, self.temp_aligned_filename)
        
        if len(list(SeqIO.parse(mafft_in, "fasta"))) == 1:
            shutil.copyfile(mafft_in, mafft_out)
        else:
            command = f"mafft --retree 2 --thread -1 --reorder {mafft_in} > {mafft_out}"
            os.system(command)

    def set_in_path(self, in_file):
        self.in_file = in_file
        
    def set_out_path(self, out_path):
        self.out_path = out_path
        self.temp_folder = os.path.join(self.out_path, "temp")
        create_folder(self.temp_folder)
        self.log_path = os.path.join(out_path, "log.txt")

    def alofi(self):
        """ seperate long_sequence_group from sequences to be aligned using repeated clustering and MSA
            long_sequence_group aligend first and all others are added using "-addfragements"
            to imporve the quality of MSA
            
            * perhaps the only method of this class that needed calling outside
            
            the result of MSA and seperated long_sequence_group will be written in out_path
            a log file will also be created in the "log" folder in out_path
        """
        ## Exception situation: only one sequence in the input fasta file: just copy the file as is the alignment
        records = list(SeqIO.parse(self.in_file, "fasta"))
        out_path_for_exception = os.path.join(self.out_path, os.path.basename(self.in_file))
        if len(records) == 1:
            shutil.copyfile(self.in_file, out_path_for_exception)
            return 1
        
        ## STEP 1: Extract the index and length of each sequence, export a pair of [index, length] for each sequence.
        if self.check_redundant_description():
            error_message = "in file {self.in_file}: Redundant sequences detected, please filter your data first"
            print(error_message)
            error_message = "Execution abandoned."
            print(error_message)
            return

        self.length_dct = self.extract_index_and_length()
        self.record_num = len(self.length_dct.values())
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 1: load basic info of sequences and lengths")

        ## STEP 2: Initialize the number of clusters n to 3.
        self.set_cluster_num(self.INITIAL_CLUSTER_NUM)
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 2: initialize clustering on lengths")

        while True:
            ## STEP 3: Perform hierarchical clustering on the lengths and retain the clustering result with n classes.
            self.cluster_result = self.hierarchical_cluster_length()

            ## STEP 4: Select the cluster with the longest average length as the potential "long sequence group."
            avg_lengths = self.count_average_length()
            longest_cluster_idx = np.argmax(avg_lengths)

            ## STEP 5: If the min length of "long sequence group" < MIN_LONG_RECORD_RATIO, set n=n+1 and return step 3.
            longest_record_idx = np.where(np.array(self.cluster_result) == longest_cluster_idx + 1)[0]
            longest_record_num = len(longest_record_idx)

            message = f"Clustering into {self.cluster_num} clusters..."
            print(message)
            

            # if there are seq no longer than 70%, re-cluster into n+1 clusters
            all_lengths = np.array(list(self.length_dct.values()))
            max_length = np.max(all_lengths)
            length_threshold = self.LONG_LENGTH_RATIO*max_length

            if ((np.min(all_lengths[longest_record_idx]) <= self.LONG_LENGTH_RATIO*max_length) or
                    (longest_record_num <= 0)):  # this may happen when there are sequences of same length
                try:
                    self.set_cluster_num(self.cluster_num + 1)
                except:
                    longest_record_num = -1  # force to break and goto step 6.
                    break   
                continue
            else:
                break
            
        if self.DEBUG_MODE:
            self.timer.record("STEP 3-5: decide long sequence group by clustering")

        ## STEP 6: If only one or too many seqs in the "long sequence group," retain seqs longer than MIN_LONG_RECORD_RATIO.
        if longest_record_num <= 1:
            longest_record_idx = [i for i in range(len(all_lengths)) if all_lengths[i] > length_threshold]

        message = f"\nClustering finished, {len(longest_record_idx)} records in the long sequence group.\n"
        print(message)
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 6: decide long sequence group by counting lengths")

        ## STEP 7: Perform multiple sequence alignment on the "long sequence group" using mafft with --reorder.
        message = "Aligning long sequence group..."
        print(message)

        self.align_and_write(longest_record_idx)
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 7: align long sequence group")
        
        ## STEP 8: calculate consensus sequence and append to msa of long group 
        long_msa_records = list(SeqIO.parse(os.path.join(self.temp_folder,self.temp_aligned_filename), "fasta"))
        consensus_long = self.nt_calculator.get_consensus_sequence(long_msa_records)
        consensus_record = SeqRecord.SeqRecord(Seq.Seq(str(consensus_long)), description="consensus", id="", name="")
        long_msa_records = long_msa_records + [consensus_record]
        SeqIO.write(long_msa_records, os.path.join(self.temp_folder,"aligned_long_sequence_group_modified.fasta"), "fasta")
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 8: calculate consensus sequence of long sequence alignment")
        
        ## STEP 9: construct distance matrix on long sequence alignment
        distance_matrix = self.nt_calculator.\
            calculate_distance_matrix(in_path = os.path.join(self.temp_folder,"aligned_long_sequence_group_modified.fasta"),
                                      out_path = self.temp_folder,
                                      method="distance_calculation_only")
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 9: construct distance matrix on long sequence alignment")
        
        ## STEP 10: perform hierarchical clustering to search for the best match of consensus sequence
        cluster_result_consensus = self.hierarchical_cluster_distance(distance_matrix)
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 10: perform clustering on long sequence alignment")
        

        ## STEP 11: keep sequences together with the consensus sequence (last) as the long_consensus_group
        consensus_group_idx = cluster_result_consensus[-1] ## TODO : bug was here raised once [IndexError: list index out of range]
        consensus_group_records = []
        record_iter = SeqIO.parse(os.path.join(self.temp_folder, self.temp_long_seqs_filename),"fasta")
                
        count = 0
        for record in record_iter:
            if cluster_result_consensus[count] == consensus_group_idx:
                consensus_group_records.append(record)
            count += 1
        
        SeqIO.write(consensus_group_records, os.path.join(self.temp_folder,"long_consensus_group.fasta"), "fasta")
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 11: decide long consensus group")
        
        ## STEP 12: align the long consensus group
        mafft_in = os.path.join(self.temp_folder,"long_consensus_group.fasta")
        mafft_out = os.path.join(self.temp_folder, "aligned_long_consensus_group.fasta")
        if len(list(SeqIO.parse(mafft_in, "fasta"))) > 1:
            command = f"mafft --maxiterate 100 --thread -1 --reorder {mafft_in} > {mafft_out}"
            os.system(command)
        else:
            shutil.copyfile(mafft_in, mafft_out)
                
        if self.DEBUG_MODE:
            self.timer.record("STEP 12: align long consensus group")

        ## STEP 13: Use the current MSA result of the "long consensus group" as a reference to align fragments.
        message = "\nAligning by --add."
        print(message)

        
        # substep 1: long_sequence_group is written and other shorter ones also need to be written into disk
        record_iter = SeqIO.parse(self.in_file, "fasta")

        consensus_group_description = [record.description for record in consensus_group_records]
        unaligned = [record for record in record_iter if record.description not in consensus_group_description]
        unaligned_short = []
        unaligned_other = []
        for record in unaligned:
            if len(record) < self.FRAGMENT_LENGTH_RATIO*max_length:
                unaligned_short.append(record)
            else:
                unaligned_other.append(record)
        
        mafft_add_other = os.path.join(self.temp_folder, self.temp_unaligned_other_filename)
        mafft_add_short = os.path.join(self.temp_folder, self.temp_unaligned_short_filename)
        
        if len(unaligned_short) != 0:
            SeqIO.write(unaligned_short,
                        os.path.join(self.temp_folder, self.temp_unaligned_short_filename),
                        "fasta")
        if len(unaligned_other) != 0: 
            SeqIO.write(unaligned_other,
                        os.path.join(self.temp_folder, self.temp_unaligned_other_filename),
                        "fasta")

        # substep 2: do MSA, adding fragments to aligned_long_sequence_group.fasta
        mafft_in = os.path.join(self.temp_folder, "aligned_long_consensus_group.fasta")
        mafft_tmp = os.path.join(self.temp_folder, "partial_added.fasta")
        mafft_out = os.path.join(self.out_path, os.path.basename(self.in_file))
        
        if os.path.isfile(mafft_add_other):
            command = f"mafft --maxiterate 100 --thread -1 --reorder --add {mafft_add_other} {mafft_in} > {mafft_tmp}"
            os.system(command)

        if os.path.isfile(mafft_add_short):
            if os.path.isfile(mafft_tmp):
                command = f"mafft --maxiterate 100 --thread -1 --reorder --addfragments {mafft_add_short} {mafft_tmp} > {mafft_out}"
            else:
                command = f"mafft --maxiterate 100 --thread -1 --reorder --addfragments {mafft_add_short} {mafft_in} > {mafft_out}"
            os.system(command)

        if not os.path.isfile(mafft_out):
            if os.path.isfile(mafft_tmp):
                shutil.copyfile(mafft_tmp, mafft_out)
            else:
                shutil.copyfile(mafft_in,  mafft_out)
        
        message = "\n==========" \
            "\nMSA done by separately aligning and adding fragments" \
            f"\nSize of long consensus group:\t {len(consensus_group_records)}"
        print(message)
        
        if self.DEBUG_MODE:
            self.timer.record("STEP 13: align whole dataset by adding unaligned sequences")
    