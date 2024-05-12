from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import numpy as np
from itertools import product


class Sequence_indexer:
    def __init__(self, len_kmer):
        self.kmer_dict = self.create_kmer_dictionary(len_kmer)
        self.len_kmer = len_kmer
        self.num_total_kmer = 4**self.len_kmer
        self.wobbles = "RYKMSWBDHVN"
    
    @staticmethod
    def create_kmer_dictionary(k):
        base_list = list("ACGT") # all possible bases
        kmer_iter = base_list
        kmer_iter = list(product(base_list, repeat=k))
        kmer = list(map(lambda x: "".join(x), kmer_iter)) # product to list
        kmer_index = list(range(len(kmer))) # index from 0 to 4^k-1
        kmer_dict = dict([[kmer[i], kmer_index[i]] for i in range(len(kmer))]) # merge
        return kmer_dict
    
    def clean_record(self, record, keep=False):
        """ if there are wobble bases in a record, replace them or rid them.
        ----------
        Parameters
        - record - the record to remove continuous wobbles from
        - keep - if False, keep no wobble bases; if true, keep only one of countinuous ones.
        -------
        Returns
        - record - the record whose continuous wobbles have been replaced with one of that.
        """
        seq = str(record.seq).upper()
        if keep:  # not yet ready
            for base in "BDHVN":
                clean_seq = re.sub(f"{base}+", base, seq)
                record.seq = clean_seq
        else:
            for base in self.wobbles:
                record.seq = Seq(str(record.seq).replace(base, ""))
                
        record.seq = record.seq.upper()
        
        return record

    @staticmethod    
    def explain_kmer(kmer):
        """ if there are wobble in the sequence, explain them to multiple basic ones 
        ----------
        Parameters
        - sequence - the sequence to be explained with only ATCG other that wobbles.
        -------
        Returns
        - explained_sequences - a list of sequences using only ATCG.
        """
        base_dictionary = {"A": "A",
                           "C": "C",
                           "T": "T",
                           "G": "G",
                           "R": "GA",
                           "Y": "TC",
                           "K": "GT",
                           "M": "AC",
                           "S": "GC",
                           "W": "AT",
                           "B": "GTC",
                           "D": "GAT",
                           "H": "ACT",
                           "V": "GCA",
                           "N": "AGCT"}
        base_list = list(kmer)
        explained_base_list = list(map(lambda x: base_dictionary[x], base_list))
        explained_sequences = list(map(lambda x: "".join(x), product(*explained_base_list)))
        return explained_sequences
            
    def split_sequence_into_kmer(self, record):
        """ split the complete sequence into overlapping 5mer (step length = 1)
        ----------
        Parameters
        - record - the record to split
        -------
        Returns
        - kmer_list - the list of kmers split from original record 
        """
        length = len(record)
        kmer_list = []
        for i in range(length-(self.len_kmer-1)):
            kmer_list.append(record.seq[i:i+self.len_kmer])
        return kmer_list
    
    @staticmethod
    def is_wobble_included(seq):
        """ decide whether there are wobbles in the given sequence 
        ----------
        Parameters
        - seq - the sequence to check
        -------
        Returns
        - is_included - True if there are wobbles in the sequence, False if not
        """
        if any(base in seq for base in "RYKMSWBDHVN"):
            return True
        else:
            return False
    
    def explain_kmer_list(self, kmer_list):
        """ explain every kmer in the given list, expressing them with only ATCG 
        ----------
        Parameters
        - kmer_list - the list of kmers
        -------
        Returns
        - explained_kmer_list - the list of kmers expressed using only ATCG
        """
        explained_list = []
        for kmer in kmer_list:
            if not self.is_wobble_included(kmer):
                explained_list.append(kmer)
            else:
                explained_list += self.explain_kmer(kmer)
        return explained_list
                
    def kmer_to_index(self, kmer_list):
        """ express kmers using integers through hash table (dictionary in Python)
        ----------
        Parameters
        - kmer_list - the list of kmers
        -------
        Returns
        - index_array - the list of indices from the kmer_list according to kmer_dictionary,
            where each number in each postition represents the number of that kmer 
            i.e. [0,0,1,2] represents 2 AAAAA, 1 AAAAT and 2 AAAAC.
        """
        index_array = list(map(lambda x: self.kmer_dict[x], kmer_list))
        index_array.sort()
        return index_array
    
    def index_record(self, record):
        """ index a sequence of a record into kmer (index) list 
        ----------
        Parameters
        - record - the record to index
        -------
        Returns
        - index_array - the index result
        """
        record = self.clean_record(record)
        kmer_list = self.split_sequence_into_kmer(record)
        kmer_list = self.explain_kmer_list(kmer_list)
        index_array = self.kmer_to_index(kmer_list)
        return index_array
    
    def set_kmer(self, k):
        self.len_kmer = k
    
if __name__ == "__main__":
    """
    spliter = Sequence_indexer(len_kmer=7)
    test_record = SeqRecord("AATTCNCGGAATTCCGBN")
    kmer_list = spliter.split_sequence_into_kmer(test_record)
    kmer_list = spliter.explain_kmer_list(kmer_list)
    index_array = spliter.kmer_to_index(kmer_list)
    """
    """
    spliter = Sequence_indexer(len_kmer=7)
    test_record = SeqRecord("AATTCNCGGAATTCCGBN")
    test_record = spliter.clean_record(test_record)
    kmer_list = spliter.split_sequence_into_kmer(test_record)
    index_array = spliter.kmer_to_index(kmer_list)
    """
    
    spliter = Sequence_indexer(len_kmer=7)
    test_record = SeqRecord("AATTNAATTATATYACCGCGCNNCNCCGAGACNC")
    index_array = spliter.index_record(test_record)
    