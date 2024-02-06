# *-* coding:utf-8 *-*
# @Time:2024/2/3 15:17
# @Author:Ruijing Cheng
# @File:my_filter.py
# @Software:PyCharm

import os
from pathlib import Path
from miner_filter import Miner_filter
from datetime import datetime
import shutil
import pandas as pd


def rename_results(wd):
    file_list = os.listdir(Path(wd) / Path("results"))
    if "blast_results_controlled.fasta" in file_list and "blast_results_checked.fasta" in file_list:
        os.makedirs(Path(wd) / Path("results") / Path("not_controlled"))
        shutil.copy(Path(wd) / Path("results") / Path("blast_results_checked.fasta"),
                    Path(wd) / Path("results") / Path("not_controlled") / Path("blast_results_checked.fasta"))
        os.remove(Path(wd) / Path("results") / Path("blast_results_checked.fasta"))
        os.rename(Path(wd) / Path("results") / Path("blast_results_controlled.fasta"),
                  Path(wd) / Path("results") / Path("blast_results_checked.fasta"))
        print("Copied blast_result_checked.fasta the not_controlled folder. ")
        print("Renamed blast_result_controlled.fasta as blast_results_checked.fasta. ")


def combine_keep_records(wd_list):
    print("Combining keep records...")
    combined_records = None
    col_list = ["taxon_name"]
    for wd in wd_list:

        name = wd.name
        col_list.append(name)
        try:
            df = pd.read_table(Path(wd) / Path("results") / Path("blast_result_kept.txt"), sep="\t")
            if combined_records is None:
                combined_records = df[["taxon_name", "subject_acc.ver"]]
            else:
                combined_records = pd.merge(combined_records, df[["taxon_name", "subject_acc.ver"]], how="outer",
                                            on="taxon_name")
            combined_records.columns = col_list
        except FileNotFoundError:
            print("%s has not been reduced." % name)

        # species = species | set(df["taxon name"])
        # print(Path(file).stem, end="\t")
        # print(df.shape[0])
    combined_records = combined_records.fillna("-")
    combined_records.to_csv(Path(wd).parent / Path("combined_records.txt"), index=False, sep="\t")
    print("Combined records save in %s" % Path(wd).parent)


def call_miner_filter(in_path, out_path, action, len_shresh, name_correction=False):
    """
    Call miner_filter, and modify input and output file names, using one thread.
    :param in_path: working directory of one marker or the parent directory of multiple working directories
    :param action: 1-control extension, 2-reduce dataset, 3-control extension, then reduce dataset
    :param len_shresh:
    :param max_num:
    :return:
    """
    # check in_path
    dir_list = os.listdir(in_path)
    dir_list = [x for x in dir_list if os.path.isdir(Path(in_path)/Path(x))]
    wd_list = []
    if 'results' in dir_list and "tmp_files" in dir_list:
        wd_list = [in_path]
    else:
        for directory in dir_list:
            sub_dir_list = os.listdir(Path(in_path)/Path(directory))
            if 'results' in sub_dir_list and "tmp_files" in sub_dir_list:
                wd_list.append(Path(in_path)/Path(directory))
    if len(wd_list) == 0:
        print("The input path is not correct.")
        print("Please provide working directory of one marker or the parent directory of multiple working directories.")

    # check out_path

    if action == 1:  # control extension
        for wd in wd_list:
            print("Control extension: %s" % wd)
            t0 = datetime.now()
            my_miner_filter = Miner_filter(wd, wd)
            my_miner_filter.control_extension(length_ratio=0.6, max_subset_size=200, gappyness_threshold=0.5)
            t1 = datetime.now()
            print("Running time: %s seconds" % (t1 - t0))
    elif action == 2:  # reduce dataset
        for wd in wd_list:
            rename_results(wd)
            print("Reduce dataset: %s" % wd)
            t0 = datetime.now()
            my_miner_filter = Miner_filter(wd, wd)
            my_miner_filter.reduce_dataset(name_correction=name_correction,  # for name correction using trns (rtrns)
                                           subsp=True, var=True, f=True,  # for species combination
                                           sp=True, cf=True, aff=True, x=True, length_threshold=len_shresh, ignore_gap=True,
                                           # for exception removal
                                           )
            t1 = datetime.now()
            print("Running time: %s seconds" % (t1 - t0))
        combine_keep_records(wd_list)
    elif action == 3:  # control extension, then reduce dataset
        for wd in wd_list:
            my_miner_filter = Miner_filter(wd, wd)
            print("Control extension: %s" % wd)
            t0 = datetime.now()
            my_miner_filter.control_extension(length_ratio=0.6, max_subset_size=200, gappyness_threshold=0.5)
            t1 = datetime.now()
            print("Running time: %s seconds" % (t1 - t0))

            rename_results(wd)

            print("Reduce dataset: %s" % wd)
            my_miner_filter.reduce_dataset(name_correction=False,  # for name correction using trns (rtrns)
                                           subsp=True, var=True, f=True,  # for species combination
                                           sp=True, cf=True, aff=True, x=True, length_threshold=len_shresh,
                                           ignore_gap=True,
                                           # for exception removal
                                           )
            t2 = datetime.now()
            print("Running time: %s seconds" % (t2 - t1))
        combine_keep_records(wd_list)
