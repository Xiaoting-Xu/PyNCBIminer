from itertools import product
from Bio import SeqIO
import numpy as np
import ctypes
import re
import os

class nt_Calculator:
    def __init__(self):
        self.base_equality_dict = self.create_base_equality_dict()
    
    @staticmethod
    def create_base_equality_dict():
        base_equality_dict = {    
"AA":1, "AT":0, "AC":0, "AG":0, "AR":0.5, "AK":0, "AM":0.5, "AS":0, "AW":0.5, "AB":0, "AD":0.3333333333333333, "AH":0.3333333333333333, "AV":0.3333333333333333, "AN":0.25, 
"TA":0, "TT":1, "TC":0, "TG":0, "TR":0, "TK":0.5, "TM":0, "TS":0, "TW":0.5, "TB":0.3333333333333333, "TD":0.3333333333333333, "TH":0.3333333333333333, "TV":0, "TN":0.25, 
"CA":0, "CT":0, "CC":1, "CG":0, "CR":0, "CK":0, "CM":0.5, "CS":0.5, "CW":0, "CB":0.3333333333333333, "CD":0, "CH":0.3333333333333333, "CV":0.3333333333333333, "CN":0.25, 
"GA":0, "GT":0, "GC":0, "GG":1, "GR":0.5, "GK":0.5, "GM":0, "GS":0.5, "GW":0, "GB":0.3333333333333333, "GD":0.3333333333333333, "GH":0, "GV":0.3333333333333333, "GN":0.25, 
"RA":0.5, "RT":0, "RC":0, "RG":0.5, "RR":0.5, "RK":0, "RM":0, "RS":0, "RW":0, "RB":0, "RD":0, "RH":0, "RV":0, "RN":0.25, 
"KA":0, "KT":0.5, "KC":0, "KG":0.5, "KR":0, "KK":0.5, "KM":0, "KS":0, "KW":0, "KB":0, "KD":0, "KH":0, "KV":0, "KN":0.25, 
"MA":0.5, "MT":0, "MC":0.5, "MG":0, "MR":0, "MK":0, "MM":0.5, "MS":0, "MW":0, "MB":0, "MD":0, "MH":0, "MV":0, "MN":0.25, 
"SA":0, "ST":0, "SC":0.5, "SG":0.5, "SR":0, "SK":0, "SM":0, "SS":0.5, "SW":0, "SB":0, "SD":0, "SH":0, "SV":0, "SN":0.25, 
"WA":0.5, "WT":0.5, "WC":0, "WG":0, "WR":0, "WK":0, "WM":0, "WS":0, "WW":0.5, "WB":0, "WD":0, "WH":0, "WV":0, "WN":0.25, 
"BA":0, "BT":0.3333333333333333, "BC":0.3333333333333333, "BG":0.3333333333333333, "BR":0, "BK":0, "BM":0, "BS":0, "BW":0, "BB":0.3333333333333333, "BD":0, "BH":0, "BV":0, "BN":0.25, 
"DA":0.3333333333333333, "DT":0.3333333333333333, "DC":0, "DG":0.3333333333333333, "DR":0, "DK":0, "DM":0, "DS":0, "DW":0, "DB":0, "DD":0.3333333333333333, "DH":0, "DV":0, "DN":0.25, 
"HA":0.3333333333333333, "HT":0.3333333333333333, "HC":0.3333333333333333, "HG":0, "HR":0, "HK":0, "HM":0, "HS":0, "HW":0, "HB":0, "HD":0, "HH":0.3333333333333333, "HV":0, "HN":0.25, 
"VA":0.3333333333333333, "VT":0, "VC":0.3333333333333333, "VG":0.3333333333333333, "VR":0, "VK":0, "VM":0, "VS":0, "VW":0, "VB":0, "VD":0, "VH":0, "VV":0.3333333333333333, "VN":0.25, 
"NA":0.25, "NT":0.25, "NC":0.25, "NG":0.25, "NR":0.25, "NK":0.25, "NM":0.25, "NS":0.25, "NW":0.25, "NB":0.25, "ND":0.25, "NH":0.25, "NV":0.25, "NN":0.25, 
"A-":0, "T-":0, "C-":0, "G-":0, "R-":0, "K-":0, "M-":0, "S-":0, "W-":0, "B-":0, "D-":0, "H-":0, "V-":0, "N-":0.25, 
"-A":0, "-T":0, "-C":0, "-G":0, "-R":0, "-K":0, "-M":0, "-S":0, "-W":0, "-B":0, "-D":0, "-H":0, "-V":0, "-N":0.25,
"AY":0, "TY":0.5, "CY":0.5, "GY":0, "RY":0, "YA":0, "YT":0.5, "YC":0.5, "YG":0, "YR":0, "YY":0.5, "YK":0, "YM":0, "YS":0, "YW":0, "YB":0, "YD":0, "YH":0, "YV":0, "YN":0.25, "Y-":0, "KY":0, "MY":0, "SY":0, "WY":0, "BY":0, "DY":0, "HY":0, "VY":0, "NY":0.25, "-Y":0, 
        }
        
        return base_equality_dict
    
    
    @staticmethod
    def identify_gap(seq_array, column):
        """ identify if the [column] column of [seq_array] contains any gap
        ----------
        Parameters
        - seq_array - the numpy.array form of supermatrix
        - column - the index of column to be identified
        -------
        Returns
        - contains_gap - whether this column contains any gap
        """
        contains_gap = "-" in seq_array[:, column]
        return contains_gap
    
    
    def identify_commen_ends(self, num_line, record_iter):
        """ identify the commen 3' end and 5' end of first [num_line] lines from a list of records
            seqs are shortened until meet a column with no gap
        ----------
        Parameters
        - num_line - identify the commen ends of first [num_line] lines from a list of records 
        - record_iter - the Fasta_iterator containing records to be identified
        -------
        Returns
        - five_commen_end - the left-most position where all [num_line] lines are NOT gap
        - three_commen_end - the right-most position where all [num_line] lines are NOT gap
        """
        ## STEP 1: transform data structure to numpy.array for convenience
        try:
            seq_array = np.array([record.seq for record in record_iter])
        except:
            seq_array = np.array([list(record) for record in record_iter])
        seq_array = seq_array[:num_line]
        
        ## STEP 2: identify commen ends at 5' end
        five_commen_end = 0
        while True:
            if not self.identify_gap(seq_array, five_commen_end):
                break
            five_commen_end += 1
            if five_commen_end == seq_array.shape[1]:
                return [None, None]
            
        ## STEP 3: identify commen ends at 3' end
        three_commen_end = seq_array.shape[1] - 1  # the last colomn
        while True:
            if not self.identify_gap(seq_array, three_commen_end):
                break
            three_commen_end -= 1
            
        return [five_commen_end, three_commen_end+1]
        
    
    @staticmethod    
    def create_kmer_dictionary(k):
        base_list = list("ACGT") # all possible bases
        kmer_iter = base_list
        kmer_iter = list(product(base_list, repeat=k))
        kmer = list(map(lambda x: "".join(x), kmer_iter)) # product to list
        kmer_index = list(range(len(kmer))) # index from 0 to 4^k-1
        kmer_dict = dict([[kmer[i], kmer_index[i]] for i in range(len(kmer))]) # merge
        return kmer_dict

    
    def calculate_PI(self, record1, record2, record1_ends=None, record2_ends=None):
        """ calculate pairwise identity of two aligned nucleotide sequences
        ----------
        Parameters
        - seq1, seq2 - two records to be compared
        -------
        Returns
        - identity - percentage of equivalent basesin the whole sequence
        """
        if not isinstance(record1, str):
            seq1 = str(record1.seq).upper()
        else:
            seq1 = record1.upper()
            
        if not isinstance(record2, str):
            seq2 = str(record2.seq).upper()
        else:
            seq2 = record2.upper()
            
        matches = 0
        length = len(seq1)
        
        if not record1_ends:
            for i in range(length):
                if seq1[i] != "-":
                    break
            for j in range(length-1, 0-1, -1):
                if seq1[j] != "-":
                    break
            record1_ends = [i, j]
            
        if not record2_ends:
            for i in range(length):
                if seq2[i] != "-":
                    break
            for j in range(length-1, 0-1, -1):
                if seq2[j] != "-":
                    break
            record2_ends = [i, j]
        
        commen_ends = [max(record1_ends[0],record2_ends[0]), min(record1_ends[1],record2_ends[1])]
        if commen_ends[0]>=commen_ends[1]:
            return 0
        seq1 = seq1[commen_ends[0]:commen_ends[1]]
        seq2 = seq2[commen_ends[0]:commen_ends[1]]
        length = len(seq1)

        for i in range(length):
            base1 = seq1[i]
            base2 = seq2[i]
            if base1 == "-" and base2 == "-":
                length -=1
            else:
                matches += self.base_equality_dict[base1+base2]

        identity = matches / length
            
        return identity


    def calculate_avg_PI(self, records, DEBUG_MODE=False):
        """ calculate pairwise identity among a group of sequences and get average
        ----------
        Parameters
        - records - a list of sequences to get average pairwise identity
        -------
        Returns
        - avg_PI - average pairwise identity of all sequences in the given group
        """
        if isinstance(records, str):
            if os.path.isfile(records):
                records = SeqIO.parse(records, "fasta")
            
        records = list(records)
        length = len(records)
        PIs = []
        for i in range(length):
            for j in range(i + 1, length):
                PI = self.calculate_PI(records[i], records[j])
                PIs.append(PI)
                if DEBUG_MODE:
                    print(f"record {i},record {j}: {PI}")
        avg_PI = np.mean(PIs)
        return avg_PI

    @staticmethod
    def get_position_in_alignment(record, position_in_raw, left_to_right=True):
        """ get the position of given number in alignment, also the number of base other than gaps
                for example, the raw sequence is AACCTGG, the position of T is 5 (from left to right)
                while in an alignment, the sequence may be like A-ACC-TG-G, then the corresponsing position should be 7
                so given record A-ACC-TG-G, position_in_raw 5, this method/function will return 7
        ----------
        Parameters
        - record - Bio.SeqRecord.SeqRecord: the aligned sequence (contains gaps)
        - position_in_raw - the position of base when gaps are excluded
        - left_to_right - if True, count the position from left to right, else from right to left
        -------
        Returns: int
        - position_in_alignment - the position of base when gaps are included.
        """
        direction = 1 if left_to_right else -1
        sequence = record.seq
        count = 0
        
        for i in range(len(record)):
            position_i = int(i*direction + 0.5*(direction-1)) # adjust direction
            if sequence[position_i] != "-":
                count += 1
                if count == position_in_raw:
                    position_in_alignment = i+1
                    break
        try:        
            return position_in_alignment
        except Exception as e:
            print("\n\n====================================")
            print(e)
            print(record)
            print(f"position_in_raw: {position_in_raw}")
            print("==================================\n\n")
            
            
    def get_consensus_sequence(self, record_iter):
        
        if isinstance(record_iter, str):
            if os.path.isfile(record_iter):
                record_iter = SeqIO.parse(record_iter, "fasta")
        
        seq_array = np.array([list(str(record.seq).upper()) for record in record_iter])
        seq_array = self.remove_gaps_in_ends(seq_array)
        
        consensus_sequence = ""
        for i in range(seq_array.shape[1]):
            dominate_base = self.get_dominate_base(seq_array[:,i])
            consensus_sequence += dominate_base
        return consensus_sequence


    @staticmethod
    def get_dominate_base(base_array):
        """ wobble base look up table:
            A --> adenosine           M --> A C (amino)
            C --> cytidine            S --> G C (strong)
            G --> guanine             W --> A T (weak)
            T --> thymidine           B --> G T C
            U --> uridine             D --> G A T
            R --> G A (purine)        H --> A C T
            Y --> T C (pyrimidine)    V --> G C A
            K --> G T (keto)          N --> A G C T (any)
        """
        base_count = {"A":0, "C":0, "G":0, "T":0, "-":0}
        base_dict =  {"A":[["A",1]], "C":[["C",1]], "G":[["G",1]], "T":[["T",1]],
                      "R":[["G",0.5],["A",0.5]], "Y":[["T",0.5],["C",0.5]], "K":[["G",0.5],["T",0.5]], 
                      "M":[["A",0.5],["C",0.5]], "S":[["G",0.5],["C",0.5]], "W":[["A",0.5],["T",0.5]],
                      "B":[["G",0.33333],["T",0.33333],["C",0.33333]],
                      "D":[["G",0.33333],["A",0.33333],["T",0.33333]],
                      "H":[["A",0.33333],["C",0.33333],["T",0.33333]],
                      "V":[["G",0.33333],["C",0.33333],["A",0.33333]],
                      "N":[["A",0.25],["C",0.25],["G",0.25],["T",0.25]],
                      "-":[["-",1]], " ":[]}
        
        base_summary = {"A":"A", "C":"C", "G":"G", "T":"T",
                        "AG":"R", "CT":"Y", "GT":"K", "AC":"M", "CG":"S", "AT":"W",
                        "CGT":"B", "AGT":"D", "ACT":"H", "ACG":"V", "ACGT":"N"}
        
        for base in base_array:
            score_list = base_dict[base]
            for score in score_list:
                base_count[score[0]] += score[1]
                
        max_occurrence = max(base_count.values())
        dominate_base_list = []
        for key, value in base_count.items():
            if value == max_occurrence:
                dominate_base_list.append(key)

        if "-" in dominate_base_list:
            if len(dominate_base_list) > 1:
                dominate_base = "N"
            else:
                dominate_base = "-"
        else:
            dominate_base = base_summary["".join(dominate_base_list)]
                
        return dominate_base
        
    
    def calculate_distance_matrix(self, in_path, out_path, method="alignment", len_kmer=None, fast=False):
        """ calculate pairwise distance of given records, by measuring different metric (only MSA available now)
        ----------
        Parameters
        - in_path - the input fatsa file of records to calculate pairwise distance.
        - out_path - destination folder of output, where log files and result will be written.
        - method - only "aligment" is available now, using MSA to calculate pairwise identity as distance matrix.
        - len_kmer - the length of kmer while using kmer filter to estimate the pairwise similarity, currently no use.
        -------
        Returns
        - distance_matrix - the matrix that stores the pairwise distance (only pairwise identity now)
        """
        ## TODO protection: if methods=="kmer" and len_kmer is not a valid positive integer
        ## TODO: if input is an alignment then just get distance matrix
        ## TODO: parllelize
        
        if method == "alignment":
            ## STEP 1: preparation
            try:
                basename = os.path.basename(in_path)
                out_filename = os.path.join(out_path, "distance_matrix_" + basename)
                if fast:
                    command = f"mafft --retree 1 {in_path} > {out_filename}"
                else:
                    command = f"mafft --auto {in_path} > {out_filename}"
                os.system(command)
                records = list(SeqIO.parse(out_filename, "fasta"))
            except:
                records = in_path
        else:
            records = list(SeqIO.parse(in_path, "fasta"))
        
        ## STEP 2: identify ends for each sequence
        all_ends = []  # [[seq1_five_end, seq1_three_end],[seq2_five_end, seq2_three_end]]
        for record in records:
            all_ends.append(self.identify_ends(record))
        
        ## STEP 3: calculate distance
        num_records = len(records)
        identity_matrix = np.ones((num_records, num_records))
        for i in range(num_records):
            for j in range(i+1, num_records):
                pairwise_identity = self.calculate_PI(records[i], records[j], all_ends[i], all_ends[j])
                identity_matrix[i,j] = pairwise_identity
                identity_matrix[j,i] = pairwise_identity
        
        distance_matrix = 1 - identity_matrix
        return distance_matrix

    @staticmethod
    def remove_gaps_in_ends(seq_array):
        for i in range(seq_array.shape[0]):
            for j in range(seq_array.shape[1]):
                if seq_array[i,j] == "-":
                    seq_array[i,j] = " "
                else:
                    break
        for i in range(seq_array.shape[0]):
            for j in range(seq_array.shape[1]-1, -1, -1):
                if seq_array[i,j] == "-":
                    seq_array[i,j] = " "
                else:
                    break
        return seq_array

    
    def calculate_sum_of_pairs_score(self, record_path, conservative=True):
        penalty = -2 if conservative else -1
        gap_open = 2
        gap_extend = 1
        
        if os.path.isfile(record_path):
            record_iter = SeqIO.parse(record_path, "fasta")
            records = np.array([list(str(record.seq).upper()) for record in record_iter])
        else:
            records = record_path
        
        sum_of_pairs = 0
        for col in range(records.shape[1]):
            column = records[:,col]
            for i in range(len(column)-1):
                base1 = column[i]
                if base1 == "-":
                    continue
                for j in range(i+1, len(column)):
                    if column[j] == "-":
                        continue
                    
                    base_relation = self.base_equality_dict[base1+column[j]]
                    sum_of_this_pair = penalty if base_relation==0 else base_relation
                    sum_of_pairs += sum_of_this_pair
        
        records = self.remove_gaps_in_ends(records)
        for i in range(records.shape[0]):
            last_base = None
            for j in range(records.shape[1]):
                if records[i,j] == "-":
                    if last_base:
                        sum_of_pairs -= gap_extend
                    else:
                        sum_of_pairs -= gap_open
                    last_base = "-"
                else:
                    last_base = None # not a gap anyway
                    
        return sum_of_pairs
    
    @staticmethod
    def identify_ends(record):
        if not isinstance(record, str):
            sequence = str(record.seq)
        else:
            sequence = record
        length = len(sequence)
        
        for i in range(0, length, 1):
            if sequence[i] != "-":
                five_end = i
                break
        
        for i in range(length-1,0-1,-1):
            if sequence[i] != "-":
                three_end = i
                break
            
        return [five_end, three_end]
    
    
    
    
    
    
    