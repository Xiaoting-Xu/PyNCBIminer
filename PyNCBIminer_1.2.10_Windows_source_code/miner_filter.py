import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from operator import add
from functools import reduce
import pandas as pd
import shutil
import re
import traceback
import time

from functional import check_inpath_validity, check_outpath_validity, get_file_handles, create_folder, get_checked_path, get_fasta
from message_logger import MessageLogger
from nt_calculator import nt_Calculator


class Miner_filter:
    """ CLass Miner_filter - filter retrieved seqs, including three main functions:
        1. <func> control_extension: control the extension of all seqs, trim if necessary.
        2. <func> combine_species: combine specific records to corresponding species (subsp., var., f.)
        3. <func> remove_exceptional_records: remove specific records (sp., cf., aff., x, and short ones)
        4. <func> reduce_dataset: select best marker for each taxon.
        
    Public:
        <func> Miner_filter - construction method
        <func> control_extension: control the extension of all seqs, trim if necessary.
        <func> combine_species: combine specific records to corresponding species (subsp., var., f.)
        <func> remove_exceptional_records: remove specific records (sp., cf., aff., x, and short ones)
        <func> reduce_dataset: select best marker for each taxon.
        <func> tnrs_name_correction: to correct names of organisms based on tnrs (or rTNRS)
        <func> get_consensus_dict: "get method" to get taxa consensus dict
        
    Private:
        <attr> in_path - input path of this class, usually the output (working directory) of "retrived sequences"
        <attr> out_path - output path of this class, usually the output (working directory) of "retrived sequences"
        <attr> tmp_path - destination folder of temporary files
        <attr> log_path - destination folder of log files
        <attr> logger - an instance of MessageLogger, to show or write log messages
        <attr> num_query - the number of sequences in the query
        <attr> quality_control_max_size_subset - the maximum size when spliting subsets
        <attr> <class> nt_calculator - some functions for pairwise identity calculation, etc.
        <attr> taxa_consensus_dict - a dictionary storing consensus sequence for each taxon
        
        <func> get_input_filename - get the most valid (suitable) filename as input fasta file
        
        ================================== for <func> control_extension ===========================================
        <func>[static] get_upper_taxonomic_unit - get the upper taxonomic group of the given taxonomy
        <func> split_by_genus - split the fasta file according to genus of the records
        <func> split_by_length - split the fastas according to relative length of the records
        <func> split_large_subset - split the fastas according to number of the records
        <func> align_subset - align the subsets split previously
        <func> remove_erroneous_extension - remove extension if the extension part is too gappy in the MSA
        
        ================================== for <func> reduce_dataset ==============================================
        <func>[static] get_length_without_wobble - get the number of ATCGs in the sequence
        <func>[static] check_record_coverage - check if given segment contains too many gaps (NOT IN USE)
        <func> calculate_consensus_dict - calculate consensus sequence for each taxon and store to a dictionary
        <func> count_num_query - count the number of sequences in the first query
        <func> align_long_seq - align sequences from the blast_result_long (chosen ones) using --add (NOT IN USE)
        <func> evaluate_seq - evaluate each sequence for representative selecton
        <func> select_seqs_to_keep - select sequences from each taxon according to 3 criteria and write into tsv
        <func> save_selected_seqs - save only the selected sequences into fasta according to the saved csv     
        """
    
    def __init__(self, in_path, out_path, DEBUG_MODE=False):
        """ construction and method
        ----------
        Parameters
        - in_path - input path of this class, usually the output (working directory) of "retrived sequences"
        - out_path - output path of this class, usually the output (working directory) of "retrived sequences"
        """
        self.__in_path = in_path  # usually the output folder (working directory) of the "retrived sequences"
        self.__out_path = out_path
        self.__tmp_path = os.path.join(out_path, "tmp_files")
        self.__log_path = out_path
        self.__logger = MessageLogger(self.__log_path)
        self.__num_query = 0  # the number of sequences in the query, initiated as 0
        self.__quality_control_max_size_subset = 0
        self.__nt_calculator = nt_Calculator()
        self.__taxa_consensus_dict = {} # {"Magnolia coco": "AATTCCGG", "taxon 2": "AATCGCCTT", ...}
        self.DEBUG_MODE = DEBUG_MODE
        
        files = ["blast_results_non_duplicate.fasta",
                 "blast_results_exception_removed.fasta",
                 "blast_results_checked_seq_info_modified.txt",
                 "blast_result_kept.txt",
                 "blast_results_filtered.fasta"]
        
        file_time = False
        backup_folder = os.path.join(out_path,"results","history_backup")
        create_folder(backup_folder)
        for file in os.listdir(os.path.join(out_path,"results")):
            if file in files:
                if not file_time:
                    file_time = time.localtime(os.stat(os.path.join(out_path,"results",file)).st_mtime)
                    file_time = time.strftime("%Y-%m-%d.%H.%M'.%S''", file_time)
                    backup_folder = os.path.join(backup_folder, file_time)
                    create_folder(backup_folder)
                shutil.move(os.path.join(out_path,"results",file), os.path.join(backup_folder, file))
        
        tmp_files = ["consensus_calculation",
                     "blast_result_kept.txt",
                     "blast_result_long.fasta",
                     "blast_result_long.txt"]
        for path in os.listdir(self.__tmp_path):
            if path in tmp_files:
                path = os.path.join(self.__tmp_path, path)
                try:
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                except:
                    pass
                
    def tnrs_name_correction(self):
        """ Main author: Yuxuan Wang, Yang Yi
            to correct names of organisms based on tnrs (or rTNRS), original file will be back up as seq_info_ori.txt
            dependencies required: folder 'TNRS_dep' should be placed in the same path with this script (miner_filter.py)
        """
        if "R_HOME" not in os.environ:
            os.environ["R_HOME"] = "./TNRS_dep"
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        
        pandas2ri.activate()
        tnrs = importr("TNRS")
        blast_info_path = os.path.join(self.__in_path,"results","blast_results_checked_seq_info.txt")
        df = pd.read_table(blast_info_path)
        df_main_info = df[['accession', 'organism']]
        result = tnrs.TNRS(df_main_info, sources="wcvp", classification="wfo", mode="resolve", matches="best", skip_internet_check=True)
        with localconverter(ro.default_converter + pandas2ri.converter):
            pdf = ro.conversion.rpy2py(result)
        pdf_selected = pdf.iloc[:, [1,32,33,34,37,40]]
        pdf_selected.to_csv(os.path.join(self.__tmp_path,'temp_TNRS.txt'), index=False, sep='\t')
        dfc = pd.merge(df, pdf_selected, left_on='organism', right_on='Name_submitted', how='left')
        dfc.rename(columns={'organism': 'organism_ori', 'Accepted_name': 'organism'}, inplace=True)
        
        out_path = os.path.join(self.__out_path,"results","blast_results_checked_seq_info.txt")
        if os.path.isfile(out_path):
            shutil.move(out_path, 
                        out_path.replace("blast_results_checked_seq_info.txt","blast_results_checked_seq_info_ori.txt"))
        dfc.to_csv(out_path, index=False, sep='\t')

        
    def remove_exceptional_records(self, in_path=None, out_path=None,
                                   sp=True, cf=True, aff=True, x=True, 
                                   length_threshold=0, ignore_gap=True):
        """ remove specific records (sp., cf., aff., x, and short ones)
        ----------
        Parameters
        - in_path - input path of fasta file (for independent calling), or None to remain the same as this class
        - out_path - output path of fasta file (for independent calling), or None to remain the same as this class
        - sp  - if True, records with " sp. "  will be removed
        - cf  - if True, records with " cf. "  will be removed
        - aff - if True, records with " aff. " will be removed
        - x   - if True, records with " x "    will be removed
        - length_threshold - sequences shorter than this value will be removed
        - ignore_gap - if True, gaps will be ignored when counting lengths of sequences
            for example:
                for record.seq = "AATT-CGGAA", length_threshold = 10, ignore_gap = True,
                    this length will be recorded as 9, thus lower than threshold, will be removed
                while for record.seq = "AATT-CGGAA", length_threshold = 10, ignore_gap = False,
                    this length will be recorded as 10, thus not lower than threshold, will be kept.
        """
        ## STEP 1: get paths ready, independent calling allowed
        if in_path:
            this_in_path = in_path
        else:
            curr_step = "blast_results_exception_removed.fasta"
            this_in_path = os.path.join(self.__in_path, "results", self.__get_input_filename(curr_step))
            
        if out_path:
            this_out_path = out_path
        else:
            this_out_path = os.path.join(self.__out_path, "results")
        
        df_path = os.path.join(self.__in_path, "results", self.__get_info_csv())
        df_records_info = pd.read_csv(df_path, sep="\t")
        out_file_name = "blast_results_exception_removed.fasta"
        
        ## STEP 2: set special tokens to be removed
        remove_list = []
        if sp:
            remove_list.append("_sp._")
            remove_list.append("_sp_")
        if cf:
            remove_list.append("_cf._")
            remove_list.append("_cf_")
        if aff:
            remove_list.append("_aff._")
            remove_list.append("_aff_")
        if x:
            remove_list.append("_x._")
            remove_list.append("_x_")
            
        ## STEP 3: keep seqs other than sequences with specified token or lower than length threshold
        record_iter = SeqIO.parse(this_in_path, "fasta")
        filtered_records = []
        for record in record_iter:
            if ignore_gap:
                record.seq = Seq(str(record.seq).replace("-", ""))
            
            accession = record.description.split("|")[0].split(":")[0]
            # print(df_records_info.loc[df_records_info["accession"]==accession]["organism"])
            frame = list(df_records_info.loc[df_records_info["accession"]==accession]["organism"])
            if len(frame) == 1:
                organism = frame[0]
            else:
                continue
            
            organism = organism.replace(" ", "_")
            if any([banned_word in organism for banned_word in remove_list]):
                continue
            
            if len(record.seq) < length_threshold:
                continue
            filtered_records.append(record)
        
        out_fasta = os.path.join(this_out_path, out_file_name)
        SeqIO.write(filtered_records, out_fasta, "fasta")
        
        df_records_info.dropna(subset=["organism"], inplace=True)
        df_records_info.to_csv(df_path, index=False, sep="\t")
        
    
    def combine_species(self, subsp=True, var=True, f=True):
        """ combine specific records(subsp, var, f) into their species 
                method: by renaming the column ['organism'] in table blast_results_checked_seq_info.txt
                example: if subsp == True, treat V._rafinesquianum_var._affine and V._rafinesquianum equally in next steps
        ----------
        Parameters
        - subsp - if True, records with " subsp. "  will be combined to their species
        - var   - if True, records with " var. "    will be combined to their species
        - f     - if True, records with " f. "      will be combined to their species
        """
        ## STEP 1: load blast_results_checked_seq_info.txt
        df_records_info = pd.read_csv(os.path.join(self.__in_path, "results", "blast_results_checked_seq_info.txt"), sep="\t")
        
        ## STEP 2: set special tokens to be removed
        combine_list = []
        if subsp:
            combine_list.append(" subsp. ")
            combine_list.append(" subsp ")
            combine_list.append(" ssp. ")
            combine_list.append(" ssp ")
        if var:
            combine_list.append(" var. ")
            combine_list.append(" var ")
        if f:
            combine_list.append(" f. ")
            combine_list.append(" f ")
                    
        ## STEP 3: modify the info csv (blast_results_checked_seq_info.txt)
        for i in df_records_info.index:
            for word in combine_list:
                if not isinstance(df_records_info["organism"][i], str):
                    df_records_info.drop(i, axis=0, inplace=True)
                    break
                if word in df_records_info["organism"][i]:
                    organism = df_records_info["organism"][i]
                    organism = organism[:organism.index(word)]
                    df_records_info.loc[i,"organism"] = organism
        
        csv_out_path = os.path.join(self.__in_path, "results", "blast_results_checked_seq_info_modified.txt")
        df_records_info.to_csv(csv_out_path, sep="\t", index=False)

    def reduce_dataset(self,
                       name_correction=False,  # for name correction using tnrs (rtnrs)
                       consensus_value=True,
                       subsp=True, var=True, f=True,  # for species combination
                       sp=True, cf=True, aff=True, x=True, length_threshold=0, ignore_gap=True,  # for exception removal
                       max_insertion_length=20, max_insertion_num=1  # for deletion of large inserted fragment
                       ):
        """ to reduce the dataset by select the best representative sequence for each taxon
        ----------
        Parameters
        - sp  - if True, records with " sp. "  will be removed
        - cf  - if True, records with " cf. "  will be removed
        - aff - if True, records with " aff. " will be removed
        - x   - if True, records with " x "    will be removed
        - length_threshold - sequences shorter than this value will be removed
        - ignore_gap - if True, gaps will be ignored when counting lengths of sequences
        """
        if name_correction:
            self.tnrs_name_correction()
        print("Filtering... turn to tmp_files/consensus_calculation for approximate progress")
        
        df = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()), sep="\t")        
        for row_index, row in df.iterrows():
            organism = row["organism"]
            if isinstance(organism, str):
                if "'" in organism:
                    df.loc[row_index, "organism"] = organism.replace("'","")
        df.to_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()), sep="\t", index=False)
        
        self.__count_consensus_value = consensus_value
        
        self.combine_species(subsp=subsp, var=var, f=f)
        self.remove_exceptional_records(sp=sp, cf=cf, aff=aff, x=x, 
                                        length_threshold=length_threshold, ignore_gap=ignore_gap)
        self.remove_duplicate()
        if not ((max_insertion_length<=0 or max_insertion_num<=0) and not self.__count_consensus_value):
            self.__calculate_consensus_dict(length_threshold=max_insertion_length, taxa_threshold=max_insertion_num)
        self.__evaluate_seq()  # may need further discussion
        self.__select_seqs_to_keep()
        self.__save_selected_seqs()
        
    def remove_duplicate(self):
        """ remove redundant records (records with exactly same voucher name and sequence are defined redundancy)
        """
        ## STEP 1: load blast_results_checked_seq_info.txt and fasta file
        curr_step = "blast_results_non_duplicate.fasta"
        df_records_info = pd.read_csv(os.path.join(self.__in_path, "results", "blast_results_checked_seq_info.txt"), sep="\t")
        record_iter = list(SeqIO.parse(os.path.join(self.__in_path, "results", self.__get_input_filename(curr_step)),"fasta"))
        df_records_info["specimen_voucher"] = df_records_info["specimen_voucher"].fillna("unknown")
        record_iter = [record for record in record_iter if record.description.split("|")[0].split(":")[0] in list(df_records_info["accession"])]
        
        record_iter = sorted(list(record_iter), 
                             key=lambda record:list(df_records_info
                                                    .loc[df_records_info["accession"]==record.description.split("|")[0].split(":")[0]]
                                                    ["specimen_voucher"])[0])
        vouchers = [list(df_records_info.loc[df_records_info["accession"]==record.description.split("|")[0].split(":")[0]]["specimen_voucher"])[0]
                    for record in record_iter]
        vouchers = ["      " if voucher == "unknown" else voucher for voucher in vouchers]
        
        organisms = [list(df_records_info.loc[df_records_info["accession"]==record.description.split("|")[0].split(":")[0]]["organism"])[0]
                    for record in record_iter]
        
        ## STEP 2: traverse all the records and record non duplicate ones, those with vouchers first
        
        df = pd.DataFrame({"accession":[record.description.split("|")[0].split(":")[0] for record in record_iter],
                           "seq":[record.seq for record in record_iter],
                           "voucher":vouchers,
                           "organism":organisms})
        
        df = df.groupby("organism") \
               .apply(
                   lambda x: x.sort_values(by=["seq","voucher"], ascending=[True, False]) \
                              .drop_duplicates(subset=["seq","voucher"], keep='first')
               )
               
        accession_numbers = list(df["accession"])
        
        ## STEP 3: save records according to list accession_numbers
        records = [record for record in record_iter if record.description.split("|")[0].split(":")[0] in accession_numbers]
        SeqIO.write(records, os.path.join(self.__in_path, "results", curr_step),"fasta")
        
        
    def control_extension(self, length_ratio=0.6, max_subset_size=200, gappyness_threshold=0.5):
        """ control the extension of all seqs, trim if necessary
        ----------
        Parameters
        - length_ratio - records longer than this ratio will be decided as longer seqs while spliting by length
        - max_size - fasta contains more than [max_size] records will be split into smaller ones
        - gappyness_threshold - if extension with gappyness more than this number will be removed/trimmed
        """
        create_folder(os.path.join(self.__out_path, "tmp_files/extension_control"))
        self.__split_by_genus()
        self.__split_by_length(length_ratio=length_ratio)
        self.__split_large_subset(max_size=max_subset_size)
        self.__align_subset()
        self.__remove_erroneous_extension(gappyness_threshold=gappyness_threshold)
        
    def get_consensus_dict(self):
        """ get method: get taxa consensus dict
        -------
        Returns
        - self.__taxa_consensus_dict - the dictionary of consensus sequence of each taxon
        """
        return self.__taxa_consensus_dict
    
    def __get_input_filename(self, curr_step=None):
        """ get the most valid (suitable) filename as input fasta file 
        ----------
        Parameters
        - curr_step - output filename of current step, for limitation of searching space
        -------
        Returns
        - matching_filename - the most valid (suitable) filename as input fasta file
        """
        path = os.path.join(self.__in_path, "results")
        existing_files = os.listdir(path)
        files = ["blast_results_checked.fasta",
                 "blast_results_controlled.fasta",
                 "blast_results_exception_removed.fasta",
                 "blast_results_non_duplicate.fasta"]
        if not curr_step:
            curr_step = len(files)-1
        else:
            curr_step = files.index(curr_step)
            
        file_range = max(len(files)-1, curr_step-1)
        for i in range(file_range,-1,-1):
            if files[i] in existing_files:
                matching_filename = files[i]
                return matching_filename
    
    
    ## ===========================================================================================================
    ## ================================== for <func> reduce_dataset ==============================================    
    def __get_info_csv(self):
        """ check if there are modified seq info table in the path, if not, load original one
        """
        if os.path.isfile(os.path.join(self.__in_path, "results", "blast_results_checked_seq_info_modified.txt")):
            return "blast_results_checked_seq_info_modified.txt"
        else:
            return "blast_results_checked_seq_info.txt"

    @staticmethod
    def __get_length_without_wobble(record):
        """ get the number of ATCGs in the sequence (only ATCG are counted, RYKMSWBDHVN and gap '-' are ignored)
        ----------
        Parameters
        - record - the record to count length
        -------
        Returns
        - seq_length - the length of the sequence (not counting wobbles like WSPYN...)
        """
        sequence = record.seq.upper()
        for wobble in "RYKMSWBDHVN":
            sequence = str(sequence).replace(wobble, "")
        seq_length = len(sequence)
        return seq_length
    
    @staticmethod
    def __check_record_coverage(record, check_start, check_end, missing_ratio=0.2):
        """ check if base of given segment contains no more than [missing_ratio] gaps
        ----------
        Parameters
        - record - the record (Bio.SeqRecord.SeqRecord) to check if double ends are covered
        - check_start - the start (5') position of the segment
        - check_end - the end (3') position of the segment
        - missing_ratio - the maximum ratio of gaps in the segment, 
            record exceeding it will be marked 0, else 1 (means double ends are covered)
            Default: 0.2
        -------
        Returns
        - coverage - whether the record meets the standard to be considered "covering double ends"
        """
        segment = record.seq[check_start:check_end]
        coverage = False
        
        if segment.count("-") <= missing_ratio*len(segment):
            coverage = True
            
        return coverage
    
    def __calculate_consensus_dict(self, length_threshold=20, taxa_threshold=1):
        """ calculate consensus sequence for each taxon and store to a dictionary
        ----------
        Parameters
        - length_threshold - for <func> remove_long_insertion:  insertion longer than this threshold will be removed
        - taxa_threshold - for <func> remove_long_insertion: insertion in at most [taxa_threshold] taxa will be removed
        """
        in_path = os.path.join(self.__in_path, "results")
        tmp_path= os.path.join(self.__tmp_path, "consensus_calculation")
        create_folder(tmp_path)
        
        df_records_info = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()), sep="\t")
        record_iter = SeqIO.parse(os.path.join(in_path, self.__get_input_filename()), "fasta")
        records_grouped = {} # {"Magnolia coco":[SeqRecord1, SeqRecord2], "taxon 2": [SeqRecord1], ...}
        records_consensus = {} # {"Magnolia coco": "AATTCCGG", "taxon 2": "AATCGCCTT", ...}
        
        for record in record_iter:
            accession = record.description.split("|")[0].split(":")[0]
            organism = list(df_records_info.loc[df_records_info["accession"]==accession]["organism"])[0]
            organism = organism.replace(" ", "_")
            records_grouped.setdefault(organism, [])
            records_grouped[organism].append(record)
            
        for taxon in records_grouped.keys():
            records_path = os.path.join(tmp_path, f"{taxon}.fasta")
            msa_path = os.path.join(tmp_path, f"{taxon}_msa.fasta")
            
            records = records_grouped[taxon]
            if len(records) <= 3:
                continue
            SeqIO.write(records, records_path, "fasta")
            
            command = f"mafft --auto --thread -1 --reorder {records_path} > {msa_path}"
            os.system(command)
            
            self.__remove_long_insertion(taxon, length_threshold, taxa_threshold)
        
        for taxon in records_grouped.keys():   
            records = records_grouped[taxon]
            if len(records) <= 3:
                continue
            msa_path = os.path.join(tmp_path, f"{taxon}_msa.fasta")
            consensus_sequence = self.__nt_calculator.get_consensus_sequence(SeqIO.parse(msa_path, "fasta"))
            records_consensus[taxon] = consensus_sequence
            
        self.__taxa_consensus_dict = records_consensus
    
    def __count_num_query(self, query_number=1):
        """ count the number of sequences in the specified query
        ----------
        Parameters
        - query_number - the number index of query to count
        """
        initial_msa = os.path.join(self.__in_path, "parameters", "ref_msa", f"msa_queries_{query_number}.fasta")
        record_iter = SeqIO.parse(initial_msa, "fasta")
        self.__num_query = len(list(record_iter))
        
        # write log file
        msg = f"in <func> count_num_query:\n  {self.__num_query} sequences identified as the reference"
        self.__logger.write_message(msg)

    def __align_long_seq(self):
        """ align the sequences from the blast_result_long (chosen ones) using --add
        """
        blast_result = os.path.join(self.__tmp_path, "blast_result_long.fasta")
        msa_query = os.path.join(self.__in_path, "parameters", "ref_msa", "msa_queries_1.fasta")
        msa_blast_result = os.path.join(self.__tmp_path, "msa_blast_result_long.fasta")
        command = f"mafft --add {blast_result} {msa_query} > {msa_blast_result}"
        os.system(command)
        
        # write log file
        msg = f'in <func> align_long_seq:\n  Long sequences are aligned with MAFFT using command "{command}"'
        self.__logger.write_message(msg)

    def __evaluate_seq(self):
        """ get [max_num] most longest sequence of each taxon and save
        ----------
        Parameters
        - int max_num - the number of seqs that are chosen to consider
        """
        ## STEP 1: load blast result
        blast_result = os.path.join(self.__out_path, "results", self.__get_input_filename())
        record_iter = SeqIO.parse(blast_result, "fasta")
        df_records_info = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()), sep="\t")

        ## STEP 6: store information of records chosen
        # substep 1: get info (that will be used for filtering) from original blast result text file
        blast_result_txt = os.path.join(self.__in_path, "results/blast_results.txt")
        df_blast_info = pd.read_csv(blast_result_txt, sep="\t", 
                                    usecols=["subject_acc.ver", "sum_hits_score", "Source"])
        
        seq_info_txt = os.path.join(self.__in_path, "results/blast_results_checked_seq_info.txt")
        df_seq_info = pd.read_csv(seq_info_txt, sep="\t",
                                  usecols=["accession","date","journal","specimen_voucher"])
        df_seq_info.rename(columns={'accession': 'subject_acc.ver'}, inplace=True)
        
        df_blast_info = pd.merge(df_blast_info, df_seq_info, on="subject_acc.ver", how="left")
        
        # substep 2: create a file consisting of acc_num and taxon_name
        record_taxon_info = []
        record_iter = sorted(list(SeqIO.parse(blast_result, "fasta")), 
                             key=lambda record:list(df_records_info
                                                    .loc[df_records_info["accession"]==record.description.split("|")[0].split(":")[0]]
                                                    ["organism"])[0])
        curr_taxon = ""
        consensus_msa_path = os.path.join(self.__in_path, "tmp_files/consensus_calculation")
        
        for record in record_iter:
            subject_acc_ver = str(record.description.split("|")[0].split(":")[0])
            organism = list(df_records_info.loc[df_records_info["accession"]==subject_acc_ver]["organism"])[0]
            organism = organism.replace(" ", "_")
            record_length = self.__get_length_without_wobble(record) # length without wobble
            
            if organism in self.__taxa_consensus_dict.keys() and self.__count_consensus_value:
                if curr_taxon != organism:
                    aligned_record_iter = list(SeqIO.parse(os.path.join(consensus_msa_path,f"{organism}_msa.fasta"), "fasta"))
                    curr_taxon = organism
                aligned_record = [record for record in aligned_record_iter 
                                  if record.description.split("|")[0].split(":")[0] == subject_acc_ver][0]
                consensus_value = self.__nt_calculator.calculate_PI(aligned_record,
                                                                    self.__taxa_consensus_dict[organism],
                                                                    [0,len(aligned_record)],
                                                                    [0,len(aligned_record)])
            
            else:
                consensus_value = -1
                
            record_taxon_info.append([subject_acc_ver,organism,record_length,consensus_value])

        df_taxon_name = pd.DataFrame(record_taxon_info, columns=["subject_acc.ver", "taxon_name", "record_length", "consensus_value"])

        # substep 3: merge two tables (and change the order of columns) 
        df_blast_info = df_blast_info.merge(df_taxon_name)
        column_order = ["subject_acc.ver", "taxon_name", "date","journal","specimen_voucher",
                        "consensus_value", "sum_hits_score", "Source", "record_length"]
        df_blast_info = df_blast_info[column_order]  # change the order of columns
        df_blast_info.to_csv(os.path.join(self.__tmp_path, "blast_result_long.txt"),
                             index=False, sep="\t")
        
    
    def __select_seqs_to_keep(self, voucher_info_required=1, publish_required=2, newest=3):
        """ select sequences from each taxon according to all criteria and write into tsv
            only default is implemented in current version:
            ones with highest consensus value first, 
            then ones with voucher info,
            then ones with published info,
            then ones that is updated most recently,
            then ones that are longer and has higher sum_hits_score, and are found during earilier blast iertations
        """
        ## STEP 1: load input parameter
        blast_result_long_txt = os.path.join(self.__tmp_path, "blast_result_long.txt")
        df = pd.read_csv(blast_result_long_txt, sep="\t")
        df["specimen_voucher"] = df["specimen_voucher"].fillna("unknown")
        creteria = {"date":newest, "journal":publish_required,"specimen_voucher":voucher_info_required}
        creteria = dict(sorted(creteria.items(), key=lambda item: item[1])) # sort by dict.values()
        for key, value in creteria.items():
            if value == False:
                del creteria[key]
        creteria = list(creteria.keys())
        
        if "specimen_voucher" in creteria:
            for i in df.index:
                if df.loc[i,"specimen_voucher"].lower() == "unknown":
                    df.loc[i,"specimen_voucher"] = "       "
                    
        if "journal" in creteria:
            unpublished = ["unknown", "unpublished", "published only in database"]
            for i in df.index:
                journal = df.loc[i,"journal"].lower()
                for word in unpublished:
                    if word in journal:
                        df.loc[i,""] = "       "
                        break
        
        if "date" in creteria:
            df["date_f"] = pd.to_datetime(df["date"],dayfirst=True)
                   
        
        ## STEP 2: for each sequence, get rank of non_gap_length, sum_hits_score, source
        """
        column_order = ["subject_acc.ver", "taxon_name", "date", "journal", "specimen_voucher",
                        "consensus_value", "sum_hits_score", "Source", "record_length",
                        "Rank_sum_hits_score", "Rank_source", "Rank_record_length", "Sum_of_rank"]"""
        
        df["Rank_sum_hits_score"] = df["sum_hits_score"].rank(ascending=False, method="average")
        df["Rank_source"] = df["Source"].rank(ascending=True, method="average")
        df["Rank_record_length"] = df["record_length"].rank(ascending=False, method="average")
        df["Sum_of_rank"] = df[["Rank_sum_hits_score", "Rank_source", "Rank_record_length"]].sum(axis=1)

        ## STEP 3: keep only one sequence for every taxon, the kept sequence has the largest sum of rank
        keeping_frame_list = []
        name_list = list(set(list(df["taxon_name"])))
        for organism in name_list:
            df_organism = df.loc[df["taxon_name"]==organism]
            
            ## substep 1: try searching for a seq with the highest consensus value
            max_consensus_value = max(df_organism["consensus_value"])
            df_organism = df_organism.loc[df_organism["consensus_value"]==max_consensus_value]
            
            if df_organism.shape[0] == 1:
                keeping_frame_list.append(df_organism)
                continue
            
            ## substep 2: if there are multiple seqs with highest consensus value, select according to "creteria"
            df_organism = self.__filter_on_voucher_info(df_organism)
            if df_organism.shape[0] == 1:
                keeping_frame_list.append(df_organism)
                continue
            
            df_organism = self.__filter_on_journal_info(df_organism)
            if df_organism.shape[0] == 1:
                keeping_frame_list.append(df_organism)
                continue
            
            df_organism = self.__filter_on_date_info(df_organism)
            if df_organism.shape[0] == 1:
                keeping_frame_list.append(df_organism)
                continue
            
            ## substep 3: if there are still more than one, select according to ranks
            df_organism.sort_values(by="Sum_of_rank", ascending=True, inplace=True)
            df_organism.drop_duplicates(subset="taxon_name", keep="first", inplace=True)
            keeping_frame_list.append(df_organism)

        df = pd.concat(keeping_frame_list)        
        df.drop("date_f", axis=1, inplace=True)
        
        ## STEP 4: save to blast_result_kept.txt
        blast_result_kept_txt = os.path.join(self.__tmp_path, "blast_result_kept.txt")
        df.to_csv(blast_result_kept_txt, sep="\t", index=False)
        blast_result_kept_txt = os.path.join(self.__out_path, "results", "blast_result_kept.txt")
        df.to_csv(blast_result_kept_txt, sep="\t", index=False)
        
    def __remove_long_insertion(self, taxon, length_threshold=20, taxa_threshold=1):
        """ if there are large insertion in an MSA, remove the insertion fragment
        ----------
        Parameters
        - taxon - the taxon name to remove insertion
        - length_threshold - insertion longer than this threshold will be removed
        - taxa_threshold - if insertion in at most [taxa_threshold] taxa, then remove, else ignore
        """
        consensus_calculation_folder = os.path.join(self.__tmp_path, "consensus_calculation")
        for file in os.listdir(consensus_calculation_folder):
            if file.endswith("_msa.fasta"):        
                self.__nt_calculator.remove_minor_large_insertion(os.path.join(consensus_calculation_folder,file),
                                                                  length_threshold=length_threshold,
                                                                  taxa_threshold=taxa_threshold,
                                                                  keep_tmp=self.DEBUG_MODE)
        
    @staticmethod
    def __filter_on_voucher_info(df):
        """ if all records have no voucher info, do nothing, else delete ones without voucher info
        ----------
        Parameters
        - df - the dataframe of an organism that is to be filtered
        -------
        Returns
        - df - the filtered dataframe
        """
        voucher_set = list(set(list(df["specimen_voucher"])))
        if voucher_set == ["       "]:
            return df
        else:
            df = df.loc[df["specimen_voucher"]!="       "]
            return df
        
    @staticmethod
    def __filter_on_journal_info(df):
        """ if all records have no publication info, do nothing, else delete ones without publication info
        ----------
        Parameters
        - df - the dataframe of an organism that is to be filtered
        -------
        Returns
        - df - the filtered dataframe
        """
        journal_set = list(set(list(df["journal"])))
        if journal_set == ["       "]:
            return df
        else:
            df = df.loc[df["journal"]!="       "]
            return df
        
    @staticmethod
    def __filter_on_date_info(df):
        """ sort dataframe by date info
        ----------
        Parameters
        - df - the dataframe of an organism that is to be filtered
        -------
        Returns
        - df - the filtered dataframe
        """
        df = df.loc[df["date_f"]==max(df["date_f"])]
        return df
    
    def __save_selected_seqs(self):
        """ save only the selected sequences into fasta according to the saved csv
        """
        # load related files
        blast_result = os.path.join(self.__out_path, "results", self.__get_input_filename())
        record_iter = SeqIO.parse(blast_result, "fasta")
        df_records_info = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()), sep="\t")
        
        blast_result_long_txt = os.path.join(self.__tmp_path, "blast_result_kept.txt")
        df = pd.read_csv(blast_result_long_txt, sep="\t", usecols = ["subject_acc.ver"])
        keeping_list = list(df["subject_acc.ver"])
        
        # select and write
        name_error_log = os.path.join(self.__out_path, "name_error.txt")
        if os.path.isfile(name_error_log):
            os.remove(name_error_log)
        
        keeping_records = []
        for record in record_iter:
            acc = record.description.split("|")[0].split(":")[0]
            if acc in keeping_list:
                description_list = record.description.split("|")
                accession = description_list[0].split(":")[0]
                organism = list(df_records_info.loc[df_records_info["accession"]==accession]["organism"])[0]
                organism = organism.replace(" ", "_")
                
                description_list[1] = organism
                record.description = "|".join(description_list)
                record.id = record.description.split(" ")[0]
                record.name = record.id
                
                try:
                    record.description.encode("gbk")    
                    keeping_records.append(record)
                except:
                    f_err = open(name_error_log, "a")
                    msg = f"{accession}: 'gbk' can't encode description of {accession}"
                    f_err.write(msg)
                    f_err.close()
                
        self.keeping_records = keeping_records
        
        filtered_records = os.path.join(self.__out_path, "results", "blast_results_filtered.fasta")
        SeqIO.write(keeping_records, filtered_records, "fasta")
        
        # write log file
        msg = "in <func> save_selected_seqs:\n  Most qualified sequence for each taxon is saved to 'blast_results_filtered.fasta'"
        self.__logger.write_message(msg)
        
    ## ===========================================================================================================    
    ## ================================== for <func> control_extension ===========================================
    @staticmethod
    def __get_upper_taxonomic_unit(this_taxonomy, info_list):
        """ get the upper taxonomic group of the given taxonomy, for example given Magnolia, return Magnoliaceae
        ----------
        Parameters
        - this_taxonomy - the taxonomy whose upper group is to be returned
        - info_list - the information list of taxonomy, in table blast_results_checked_seq_info.txt
        -------
        Returns
        - upper_taxonomy - the name of upper taxonomic group.
        """
        for info in info_list:
            groups = info.split("|")
            for i in range(1, len(groups)):
                if this_taxonomy == groups[-i]:
                    for j in range(i+1, len(groups)):
                        if " " not in groups[-j]:
                            upper_taxonomy = info.split("|")[-j]
                            return upper_taxonomy
                        
    def __split_by_genus(self):
        """ split the fasta file according to genus of the records
        """
        ## STEP 1: set and prepare folders
        in_path = os.path.join(self.__in_path, "results/blast_results_checked.fasta")
        out_path = os.path.join(self.__out_path, "tmp_files/extension_control/split_by_genus")
        create_folder(out_path)
        
        ## STEP 2: parse the blast_results_checked.fasta and save into dictionary by genus
        taxonomy_df = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()),
                                  usecols=["accession", "taxonomy"],
                                  sep="\t")
        taxonomy_dict = dict(zip(taxonomy_df["accession"],taxonomy_df['taxonomy']))
        record_iter = SeqIO.parse(in_path, "fasta")
        genus_dict = {}
        for record in record_iter:
            accession_number = record.description.split("|")[0].split(":")[0]
            genus = [unit for unit in taxonomy_dict[accession_number].split("|") if " " not in unit][-1]
            genus_dict.setdefault(genus, [])
            genus_dict[genus].append(record)
            
        ## STEP 3: write the records (each genus a individual fasta file)
        for genus, species_list in genus_dict.items():
            SeqIO.write(species_list, os.path.join(out_path, f"{genus}.fasta"), "fasta")
            
            
    def __split_by_length(self, length_ratio=0.6):
        """ split the fastas (previously split by genus) according to relative length of the records
        ----------
        Parameters
        - length_ratio - records longer than this ratio will be decided as "longer" sequences
        """
        ## STEP 1: set and prepare folders
        in_path = os.path.join(self.__in_path, "tmp_files/extension_control/split_by_genus")
        out_path = os.path.join(self.__out_path, "tmp_files/extension_control/split_by_length")
        create_folder(out_path)
        
        ## STEP 2: for each fasta file in the in_path, split the dataset according to length ratio
        for file in os.listdir(in_path):
            file_abs_path = os.path.join(in_path, file)
            genus = os.path.splitext(file)[0]
            
            ## substep 1: get the max length
            length_list = []
            record_iter = SeqIO.parse(file_abs_path, "fasta")
            for record in record_iter:
                length_list.append(len(record.seq))
            max_length = max(length_list)
            length_threshold = max_length * length_ratio
            
            ## substep 2: split the dataset, longer than ratio and no longer than the ratio
            longer_records = []
            shorter_records = []
            record_iter = SeqIO.parse(file_abs_path, "fasta")
            for record in record_iter:
                length = len(record.seq)
                if length > length_threshold:
                    longer_records.append(record)
                else:
                    shorter_records.append(record)
                    
            ## substep 3: sava the current split result
            SeqIO.write(longer_records, os.path.join(out_path, f"{genus}_longer.fasta"), "fasta")
            if len(shorter_records) > 0:
                SeqIO.write(shorter_records, os.path.join(out_path, f"{genus}_shorter.fasta"), "fasta")
            
            
    def __split_large_subset(self, max_size=200):
        """ split the fastas (previously split by genus and length) according to number of the records
        ----------
        Parameters
        - max_size - fasta contains more than [max_size] records will be split into smaller ones
        """
        ## STEP 1: set and prepare folders
        self.__quality_control_max_size_subset = max_size
        in_path = os.path.join(self.__in_path, "tmp_files/extension_control/split_by_length")
        out_path = os.path.join(self.__out_path, f"tmp_files/extension_control/split_max_{max_size}")
        create_folder(out_path)
        
        ## STEP 2: split large fasta files
        for file in os.listdir(in_path):
            file_abs_path = os.path.join(in_path, file)
            records = list(SeqIO.parse(file_abs_path, "fasta"))
            total_size = len(records)
            
            if total_size > max_size:
                num_subset = total_size // max_size + 1
                sub_size = int(total_size / num_subset)
                
                for i in range(num_subset):
                    sub_records = records[i*sub_size : (i+1)*sub_size]
                    sub_filename = os.path.splitext(file)[0] + f"_{i}.fasta"
                    SeqIO.write(sub_records, os.path.join(out_path, sub_filename), "fasta")
                
            else:
                shutil.copyfile(file_abs_path, os.path.join(out_path, file))
                
    def __align_subset(self, add_threshold=5):
        """ align the subsets
        ----------
        Parameters
        - add_threshold - files contain seqs less than this number will use --add (refer to another MSA)
        """
        ## STEP 1: set and prepare folders
        in_path = os.path.join(self.__in_path,  f"tmp_files/extension_control/split_max_{self.__quality_control_max_size_subset}")
        out_path = os.path.join(self.__out_path, "tmp_files/extension_control/subset_MSA")
        create_folder(out_path)
        
        ## STEP 2: decide the largest fasta file (contains the biggest number of sequences)
        file_size_dict = {}
        for file in os.listdir(in_path):
            record_iter = SeqIO.parse(os.path.join(in_path, file), "fasta")
            iter_length = len(list(record_iter))
            file_size_dict[file] = iter_length
        file_size_dict = dict(sorted(file_size_dict.items(), key=lambda x:x[1], reverse=True)) # sort by length
        
        ## STEP 2: perform multiple sequence alignment
        file_waiting_list = []
        for file in os.listdir(in_path):
            file_abs_path = f"{in_path}/{file}"
            file_out_path = f"{out_path}/{os.path.splitext(file)[0]}_MSA.fasta"
            
            ## substep 1: if the fasta file contains more than [add_threshold] seqs, exec command using --auto
            if file_size_dict[file] > add_threshold:
                command = f"mafft --auto {file_abs_path} > {file_out_path}"
                os.system(command)
            
            ## substep 2: if file contains no more than [add_threshold] seqs, store these for further --add
            else:
                file_waiting_list.append(file)
                
        ## substep 3: get taxonomic information
        df_taxonomy = pd.read_csv(os.path.join(self.__in_path, "results", self.__get_info_csv()),
                                  sep="\t")
        taxonomy_genus_and_above = list(set(list(map(lambda x:x[:x.rindex("|")], df_taxonomy["taxonomy"]))))
        
        ## substep 4: for each file in waiting list, find the most related taxonomic group to --add to MSA
        for file in file_waiting_list:
            file_abs_path = f"{in_path}/{file}"
            file_out_path = f"{out_path}/{os.path.splitext(file)[0]}_MSA.fasta"
            this_taxonomy = file.split("_")[0]
            upper_unit = self.__get_upper_taxonomic_unit(this_taxonomy, taxonomy_genus_and_above)
            for ref_file in file_size_dict:
                # if the genus of the candidate file is the same as this_taxonomy
                if (self.__get_upper_taxonomic_unit(ref_file.split("_")[0], taxonomy_genus_and_above) == upper_unit
                    and ref_file != file):
                    reference = ref_file
                    file_ref_path = f"{out_path}/{ref_file.replace('.fasta','_MSA.fasta')}"
                    break
            
            try:
                command = f"mafft --add {file_abs_path} {file_ref_path} > {file_out_path}"
            except:
                warning_msg = "In file {file}: there may be error in extension check" \
                              "because no other genus from the same family can be used as reference."
                self.__logger.collect_warning(warning_msg)
                reference = list(file_size_dict.keys())[0]
                file_ref_path = f"{out_path}/{reference.replace('.fasta','_MSA.fasta')}"
                command = f"mafft --add {file_abs_path} {file_ref_path} > {file_out_path}"
            
            os.system(command)
            file_name = os.path.join(out_path, command.split(">")[-1].strip())
            record_iter = SeqIO.parse(file_name, "fasta")
            
            count = 0
            records = []
            for record in record_iter:
                count += 1
                if count >= file_size_dict[reference]:
                    break
            for record in record_iter:
                records.append(record)
            SeqIO.write(records, file_name, "fasta")
            del file_ref_path
            
        self.__logger.write_warning()
             
    def __remove_erroneous_extension(self, gappyness_threshold=0.5):
        """ remove extension if the extension part is too gappy in the MSA (so-called errorneous ones)
        ----------
        Parameters
        - gappyness_threshold - if extension with gappyness more than this number will be removed/trimmed"""
        self.__logger.write_message("Into removal.")
        record_ids = []
        
        ## STEP 1: set and prepare folders
        in_path = os.path.join(self.__in_path,  "tmp_files/extension_control/subset_MSA")
        tmp_path = os.path.join(self.__out_path, "tmp_files/extension_control")
        out_path = os.path.join(self.__out_path, "results")
        create_folder(tmp_path)
        
        ## STEP 2: load the information table results/blast_results.txt
        blast_result_path = os.path.join(self.__in_path, "results/blast_results.txt")
        df_blast_result = pd.read_csv(blast_result_path, usecols=["subject_acc.ver","s_start","s_end"], sep="\t")
        
        ## STEP 3: get the position of extension part in MSA (according to results/blast_results.txt)
        df_modification = pd.DataFrame(columns=["subject_acc.ver", "new_start", "new_end", "s_new_start", "s_new_end"])
        for file in os.listdir(in_path):
            self.__logger.write_message(f"Performing removal on {file}.")
            file_abs_path = os.path.join(in_path, file)
            record_iter = SeqIO.parse(file_abs_path, "fasta")
                
            for record in record_iter:
                accession_number = record.description.split(":")[0]
                info_line = df_blast_result.loc[df_blast_result["subject_acc.ver"]==accession_number]
                
                # explanation of some variables:
                # before extension:  s_start  ---> TTCCGG <---  s_end
                # after  extension:  r_start --> AATTCCGGAA <-- r_end
                # after  alignment: as_start ---> T-TCCG-G <--- as_end
                # after  alignment: ar_start -> AAT-TCCG-GAA <- ar_end
                
                s_start, s_end = info_line.iloc[0][["s_start","s_end"]]
                s_start, s_end = min([s_start, s_end]), max([s_start, s_end])
                r_start, r_end = record.description.split("|")[0].split(":")[1].split("_")[0].split("-") # r for blast_'R'esult
                r_start, r_end = min([r_start, r_end]), max([r_start, r_end])
                s_start, s_end, r_start, r_end = list(map(int, [s_start, s_end, r_start, r_end]))
                
                as_start = self.__nt_calculator.get_position_in_alignment(record, max(s_start-(r_start-1),1), left_to_right=True)
                ar_start = self.__nt_calculator.get_position_in_alignment(record, max(r_start-(r_start-1),1), left_to_right=True)
                as_end   = self.__nt_calculator.get_position_in_alignment(record, max(r_end-(s_end-1),1), left_to_right=False)
                ar_end   = self.__nt_calculator.get_position_in_alignment(record, max(r_end-(r_end-1),1), left_to_right=False)
                
                ## STEP 4: check gappyness and remove ones that are too gappy
                ## substep 1: check gappyness
                if ar_start == as_start:
                    gappyness_start = 0
                else:
                    gappyness_start = 1 - (s_start-r_start) / (as_start-ar_start)
                    
                if ar_end == as_end:
                    gappyness_end = 0
                else:
                    gappyness_end = 1 - (r_end-s_end) / (as_end-ar_end)
                
                new_start, new_end = None, None
                if gappyness_start > gappyness_threshold:
                    new_start = s_start-r_start+1
                    s_new_start = s_start
                if gappyness_end > gappyness_threshold:
                    new_end = s_end-r_start+1
                    s_new_end = s_end
                if any([new_start, new_end]):
                    record_ids.append(accession_number)
                    if new_start == None:
                        new_start = 1
                        s_new_start = r_start
                    if new_end == None:
                        new_end = r_end-r_start+1
                        s_new_end = r_end
                        
                    position_modification = [accession_number, new_start, new_end, s_new_start, s_new_end]
                    df_modification.loc[len(df_modification.index)] = position_modification
        
        df_modification.to_csv(os.path.join(tmp_path, "modification.txt"), index=False, sep="\t")
        
        ## STEP 5: save the raw seqs (trimmed but not aligned) and the alignment (empty columns are removed)
        blast_result_checked = os.path.join(self.__in_path, "results/blast_results_checked.fasta")
        existing_accession = []
        new_records = []
        records_to_trim = list(df_modification["subject_acc.ver"])
        record_iter = SeqIO.parse(blast_result_checked, "fasta")
        for record in record_iter:
            accession_number = record.description.split(":")[0]
            if accession_number in existing_accession:
                continue
            else:
                existing_accession.append(accession_number)
                
            if accession_number in records_to_trim:
                info_line = df_modification.loc[df_modification["subject_acc.ver"]==accession_number]
                start, end, s_start, s_end = info_line.iloc[0][["new_start","new_end","s_new_start","s_new_end"]]
                record.seq = record.seq[start-1:end]
                record.description = re.sub(r"\d+-\d+", f"{s_start}-{s_end}", record.description, count=1)
                record.id = record.description.split(" ")[0]
                record.name = record.id
                new_records.append(record)
            else:
                new_records.append(record)
                
        self.__logger.write_message("Extension control finished.")
        new_records.sort(key=lambda x:x.description)
        SeqIO.write(new_records, os.path.join(out_path, "blast_results_controlled.fasta"), "fasta")
          
    
    
if __name__ == "__main__":

    """ 延伸控制使用示例如下 
    第一步，创建一个类实例，要求指定输入文件夹和输出文件夹，要求直接为Sequence Retrieving中的Working Directory，
        可以有其他指定参数，均有默认值，我觉得应该可以设置成可变的，或许组会前和师姐或者和老师讨论下；
    第二步，执行控制延伸的函数 control_extension()。
    
    输出文件：
        多数中间文件在 tmp_files/extension_control 中，每个文件夹是一个子步骤；
        最终需要裁剪掉的延伸部分在 tmp_files/modification.txt 中；
        裁剪结果在results/blast_results_controlled.fasta 中 """
    """
    my_filter = Miner_filter("C:/Users/Yy/学习资料/Project/pyNCBIminer/Windows/simple_test/test_retrieving/ITS", 
                             "C:/Users/Yy/学习资料/Project/pyNCBIminer/Windows/simple_test/test_retrieving/ITS")
    my_filter.control_extension()"""
    
    def debug_combine_species(path):
        try:
            my_filter = Miner_filter(path, path)
            my_filter.combine_species()
        except:
            with open((f"{path}/BUG_log.txt"),"w") as f:
                f.write(str(traceback.format_exc()))
    
    def debug_reduce_dataset(path):
        try:
            my_filter = Miner_filter(path, path, DEBUG_MODE=True)
            my_filter.reduce_dataset(consensus_value=True,
                                     subsp=True, var=True, f=True,  # for species combination
                                     sp=True, cf=True, aff=True, x=True, length_threshold=20, ignore_gap=True)
            return my_filter
        except:
            with open((f"{path}/BUG_log.txt"),"w") as f:
                f.write(str(traceback.format_exc()))
            return None
        
    def debug_control_and_reduce(path):
        try:
            my_filter = Miner_filter(path, path, DEBUG_MODE=True)
            my_filter.control_extension()
            my_filter.reduce_dataset(consensus_value=True,
                                     subsp=True, var=True, f=True,  # for species combination
                                     sp=True, cf=True, aff=True, x=True, length_threshold=20, ignore_gap=True)
            return my_filter
        except:
            with open((f"{path}/BUG_log.txt"),"w") as f:
                f.write(str(traceback.format_exc()))
            return None
                
    path1 = "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Magnoliaceae/GAI1"
    path2 = "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Asterales/rbcL_20230412"   # 13881 seqs    
    path3 = "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Asterales/ITS_20230412"  # 40000+ seqs
    path4 = "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Dipsacales/rbcL_20180820_Dipsacales"
    bug_path = r"C:\Users\Yy\学习资料\Project\pyNCBIminer\测试分析\bugs\2024.3.20\_psbA-trnH"
    # my_filter = debug_reduce_dataset(bug_path)
    
    path = bug_path
    debug_control_and_reduce(bug_path)

    
    """
    my_filter = Miner_filter("单独调用时，类的输入文件夹无所谓，可以初始化为空，但输出要求与实际输出文件夹一致",
                             "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Magnoliaceae/rbcL")
    my_filter.remove_exceptional_records("C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Magnoliaceae/rbcL/results/blast_results_checked.fasta",
                                         "C:/Users/Yy/学习资料/Project/pyNCBIminer/开发数据/Magnoliaceae/rbcL")
    """