# *-* coding:utf-8 *-*
# @Time:2024/2/2 15:11
# @Author:Ruijing Cheng
# @File:blast_results_extend.py
# @Software:PyCharm

import os
from pathlib import Path
import pandas as pd
from math import floor
from Bio import SeqIO
from tools import get_query_accession


def add_all_queries2(wd):
    print("Aligning all reference sequences to calculate the missing length on the left and right side...")
    queries_file_list = os.listdir(Path(wd) / Path("parameters") / Path("ref_seq"))
    # todo: think about only 1 round of blast
    if len(queries_file_list) == 0:
        print("Could not find query sequences in the parameters folder.")
        return

    # for queries_file in queries_file_list:
    #     fw = open(Path(wd)/Path("parameters")/Path("ref_msa")/Path(queries_file), "w")
    #     for record in SeqIO.parse(Path(wd)/Path("parameters")/Path("ref_seq")/Path(queries_file), "fasta"):
    #         fw.write((">%d|"+record.description+"\n") % int(Path(queries_file).stem.split("_")[-1]))
    #         fw.write(str(record.seq)+"\n")
    #     fw.close()
    if not os.path.exists(
            Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta")) or os.path.getsize(
        Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta")) == 0:
        os.system(
            r"mafft --localpair --maxiterate 1000 %s > %s" %
            (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta"),
             Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta")))

    if len(queries_file_list) > 1:
        ref_msa_file = "msa_queries_1_to_%d.fasta" % len(queries_file_list)
        ref_msa_path = Path(wd) / Path("parameters") / Path("ref_msa") / Path(ref_msa_file)
        if not os.path.exists(ref_msa_path) or os.path.getsize(ref_msa_path) == 0:

            for n in range(2, len(queries_file_list) + 1):
                if n == 2:
                    mafft_cmd = r"mafft --multipair --addfragments %s %s > %s" % \
                                (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_2.fasta"),
                                 Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta"),
                                 Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1_to_2.fasta"))
                else:
                    mafft_cmd = r"mafft --multipair --addfragments %s %s > %s" % \
                                (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_%d.fasta" % n),
                                 Path(wd) / Path("parameters") / Path("ref_msa") / Path(
                                     "msa_queries_1_to_%d.fasta" % (n - 1)),
                                 Path(wd) / Path("parameters") / Path("ref_msa") / Path(
                                     "msa_queries_1_to_%d.fasta" % n))
                os.system(mafft_cmd)


    else:
        ref_msa_file = "msa_queries_1.fasta"
    # todo:删除序列
    return ref_msa_file


def extend_hits(df, maxlen, qreflen, missing_left, missing_right):
    # df = df.copy()
    df[["sum_hits_score", "q_start", "q_end", "s_start", "s_end"]] = \
        df[["sum_hits_score", "q_start", "q_end", "s_start", "s_end"]].apply(pd.to_numeric)
    df["hit_len"] = abs(df["s_end"] - df["s_start"]) + 1
    df["max_extension"] = maxlen - df["hit_len"]

    # auto_extension is used to add the unmatched two ends between query and subject to the subject
    df["auto_extension_left"] = df["q_start"] - 1
    df["auto_extension_right"] = qreflen - df["q_end"]

    # extra_extension extends the subject further if a query is a partial sequence of the target genetic marker
    df["extra_extension_left"] = missing_left
    df["extra_extension_right"] = missing_right

    # todo: use 1.2 * max_ref_len as maxlen?
    # todo: is int faster than float?
    # todo: extra 和 auto 的计算都是基于query，可以进行向量操作
    # todo: subject正向和反向的情况筛选之后也可以批量进行？
    if maxlen - qreflen > 0:  # maxlen > qreflen
        df["extension_left"] = df["auto_extension_left"] + df["extra_extension_left"]
        df["extension_right"] = df["auto_extension_right"] + df["extra_extension_right"]
        adjust_list = ((df["extension_left"] + df["extension_right"]) > 0) & (df["extension_left"] + df["extension_right"] > df["max_extension"])
        p_ext_l = df.loc[adjust_list, "extension_left"]/(df.loc[adjust_list, "extension_left"]+df.loc[adjust_list, "extension_right"])
        p_ext_r = 1 - p_ext_l
        df.loc[adjust_list, "extension_left"] = (df.loc[adjust_list, "max_extension"] * p_ext_l).apply(floor)
        df.loc[adjust_list, "extension_right"] = (df.loc[adjust_list, "max_extension"] * p_ext_r).apply(floor)
    else:  # maxlen < qreflen
        df["extension_left"] = df["auto_extension_left"]
        df["extension_right"] = df["auto_extension_right"]
    # todo: check the results of extension for reverse complement sequences
    df["q_extstart"] = df["q_start"] - df["extension_left"]
    df["q_extend"] = df["q_end"] + df["extension_right"]
    # positive strand
    df.loc[df["s_strand"], "s_extstart"] = df.loc[df["s_strand"], "s_start"] - df.loc[df["s_strand"], "extension_left"]
    df.loc[df["s_strand"], "s_extend"] = df.loc[df["s_strand"], "s_end"] + df.loc[df["s_strand"], "extension_right"]
    # negative strand
    # # s_end <------ s_start **** s_end - extension_left <----- s_start + extension_right
    # df.loc[~df["s_strand"], "s_extstart"] = df.loc[~df["s_strand"], "s_start"] + df.loc[~df["s_strand"], "extension_right"]
    # df.loc[~df["s_strand"], "s_extend"] = df.loc[~df["s_strand"], "s_end"] - df.loc[~df["s_strand"], "extension_left"]
    # s_end <------ s_start **** s_end - extension_right <----- s_start + extension_left
    df.loc[~df["s_strand"], "s_extstart"] = df.loc[~df["s_strand"], "s_start"] + df.loc[~df["s_strand"], "extension_left"]
    df.loc[~df["s_strand"], "s_extend"] = df.loc[~df["s_strand"], "s_end"] - df.loc[~df["s_strand"], "extension_right"]

    return df


def calculate_missing_length(wd, ref_msa_file):
    print("Calculating the missing length on the left and right side of reference sequences...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(Path(wd) / Path("parameters") / Path("ref_msa") / Path(ref_msa_file), "fasta"),
                             key_function=get_query_accession)
    all_queries_info = pd.read_table(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), sep="\t")
    all_queries_info.index = all_queries_info["ID"]
    all_queries_info["Missing_left"] = 0
    all_queries_info["Missing_right"] = 0
    for key in seq_dict.keys():
        seq = seq_dict[key].seq.upper()
        left = [seq.find("A"), seq.find("T"), seq.find("C"), seq.find("G")]
        right = [seq.rfind("A"), seq.rfind("T"), seq.rfind("C"), seq.rfind("G")]
        all_queries_info.loc[key, "Missing_left"] = min(left)
        all_queries_info.loc[key, "Missing_right"] = len(seq) - max(right) - 1
    all_queries_info.to_csv(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), sep="\t", index=False)


def blast_results_extend_main(wd, max_len):
    ref_msa_file = add_all_queries2(wd)
    calculate_missing_length(wd, ref_msa_file)
    # todo: "query_acc.ver"改成ID？
    blast_results = pd.read_table(Path(wd) / Path("results") / Path("blast_results.txt"), sep="\t", engine="python")
    # if "s_extstart" in blast_results.columns:
    #     return
    ref_info = pd.read_table(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), sep="\t", engine="python")
    extended_blast_results = []
    for blast_round in blast_results["Source"].value_counts().index:
        tmp_df = blast_results[blast_results["Source"] == blast_round].copy()
        sum_table = pd.read_table(
            Path(wd) / Path("tmp_files") / Path("BLAST_%d" % blast_round) / Path("blast_summary.txt"),
            sep='\t', engine='python')
        # extend hits
        extended_hit_tables = []
        for (name, group) in tmp_df.groupby(tmp_df["query_acc.ver"]):
            if name.startswith("Query"):
                # todo: use query accession for query_acc.ver, parse hit-tables also need to be changed
                ref_id = sum_table[sum_table["Query"] == name].iloc[0]["ID"]
                ref_id = ref_id.split("|")[0].split(":")[0]
                group["query_acc.ver"] = ref_id
            else:
                ref_id = name
            # print(ref_id)

            qreflen = ref_info[ref_info["ID"] == ref_id].iloc[0]["Sequence_length"]
            missing_left = ref_info[ref_info["ID"] == ref_id].iloc[0]["Missing_left"]
            missing_right = ref_info[ref_info["ID"] == ref_id].iloc[0]["Missing_right"]
            print("Extending sequences found by %s..." % ref_id)
            extended_hit_tables.append(extend_hits(group, max_len, qreflen, missing_left, missing_right))
        extended_hit_tables = pd.concat(extended_hit_tables)
        extended_blast_results.append(extended_hit_tables)
    extended_blast_results = pd.concat(extended_blast_results)
    extended_blast_results.to_csv(Path(wd) / Path("results") / Path("blast_results.txt"), index=False, sep="\t")




