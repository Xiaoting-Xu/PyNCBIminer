# *-* coding:utf-8 *-*
# @Time:2024/2/2 15:00
# @Author:Ruijing Cheng
# @File:hits_parse_join_select.py
# @Software:PyCharm

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime


def parse_xml_by_re(wd, in_file):
    fr = open(Path(wd) / Path(in_file), "r")
    xml = fr.read()
    fr.close()
    find_query = re.compile('<BlastOutput_query-ID>(.*?)</BlastOutput_query-ID>')
    find_hit = re.compile('<Hit>(.*?)</Hit>', re.S)
    find_hit_id = re.compile('<Hit_id>(.*?)</Hit_id>')
    find_hit_def = re.compile('<Hit_def>(.*?)</Hit_def>')  # &gt; multiple def
    find_hsp = re.compile('<Hsp>(.*?)</Hsp>', re.S)
    find_identity = re.compile('<Hsp_identity>(.*?)</Hsp_identity>')
    find_align_len = re.compile('<Hsp_align-len>(.*?)</Hsp_align-len>')
    find_gaps = re.compile('<Hsp_gaps>(.*?)</Hsp_gaps>')
    find_q_start = re.compile('<Hsp_query-from>(.*?)</Hsp_query-from>')
    find_q_end = re.compile('<Hsp_query-to>(.*?)</Hsp_query-to>')
    find_s_start = re.compile('<Hsp_hit-from>(.*?)</Hsp_hit-from>')
    find_s_end = re.compile('<Hsp_hit-to>(.*?)</Hsp_hit-to>')
    find_evalue = re.compile('<Hsp_evalue>(.*?)</Hsp_evalue>')
    find_bit_score = re.compile('<Hsp_bit-score>(.*?)</Hsp_bit-score>')

    out_file_name = in_file.split("_")[0] + "_HitTable.txt"
    fw = open(Path(wd) / Path(out_file_name), "w")
    fw.write(
        "query_acc.ver\tsubject_acc.ver\thit_id\thit_def\t%_identity\talignment_length\tgap_opens\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score\n")
    query = re.findall(find_query, xml)[0]
    hits = re.findall(find_hit, xml)
    # print("Hits number: %d" % len(hits))
    for hit in hits:
        hit_id = re.findall(find_hit_id, hit)[0]
        hit_def = re.findall(find_hit_def, hit)[0].split("&gt;")
        if len(hit_id.split("|")[-1]) in [1, 2]:
            subject = hit_id.split("|")[-2] + "_" + hit_id.split("|")[-1]
        else:
            subject = hit_id.split("|")[-2]
        columns0 = [subject + "\t" + hit_id + "\t" + hit_def[0]]
        if len(hit_def) > 1:
            i = 1
            while i < len(hit_def):
                if len(hit_def[i].split(" ", 1)) == 2:
                    multi_id = hit_def[i].split(" ", 1)[0]
                    multi_def = hit_def[i].split(" ", 1)[1]
                    # print(multi_id+" "+multi_def)
                    if len(multi_id.split("|")) == 5:
                        columns0.append(multi_id.split("|")[-2] + "\t" + multi_id + "\t" + multi_def)
                i += 1
                # print(i)
        hsps = re.findall(find_hsp, hit)
        for hsp in hsps:
            identity = re.findall(find_identity, hsp)[0]
            align_len = re.findall(find_align_len, hsp)[0]
            gaps = re.findall(find_gaps, hsp)[0]
            q_start = re.findall(find_q_start, hsp)[0]
            q_end = re.findall(find_q_end, hsp)[0]
            s_start = re.findall(find_s_start, hsp)[0]
            s_end = re.findall(find_s_end, hsp)[0]
            evalue = re.findall(find_evalue, hsp)[0]
            bit_score = re.findall(find_bit_score, hsp)[0]
            columns1 = identity + "\t" + align_len + "\t" + gaps + "\t" + q_start + "\t" + q_end + "\t" + s_start + "\t" + s_end + "\t" + evalue + "\t" + bit_score
            for columns in columns0:
                fw.write(query + "\t" + columns + "\t" + columns1 + "\n")
    fw.close()
    os.remove(Path(wd) / Path(in_file))
    print("Results saved in %s" % out_file_name)


def join_group(group, allowed_length):
    # hit with the highest bit score is kept as the start of merging
    group = group.sort_values(by=["bit_score", "alignment_length", "%_identity"], ascending=[False, False, False])
    # delete hits with different direction from the highest bit-score hit
    group = group[group["q_strand"] == group.iloc[0]["q_strand"]]
    group = group[group["s_strand"] == group.iloc[0]["s_strand"]]
    # delete hits with questionable relative position of query fragment and subject fragment
    if len(group) > 1:
        group = group[
            (group["q_start"] - group.iloc[0]["q_start"] > 0) == (group["q_end"] - group.iloc[0]["q_end"] > 0)]
        group = group[
            (group["s_start"] - group.iloc[0]["s_start"] > 0) == (group["s_end"] - group.iloc[0]["s_end"] > 0)]
        group = group[
            (group["q_start"] - group.iloc[0]["q_start"] < 0) == (group["q_end"] - group.iloc[0]["q_end"] < 0)]
        group = group[
            (group["s_start"] - group.iloc[0]["s_start"] < 0) == (group["s_end"] - group.iloc[0]["s_end"] < 0)]
    # todo: query 和 subject 位置错位
    if len(group) > 1:
        group = group[
            (group["s_start"] - group.iloc[0]["s_start"] > 0) == (group["q_start"] - group.iloc[0]["q_start"] > 0)]
        group = group[
            (group["s_start"] - group.iloc[0]["s_start"] < 0) == (group["q_start"] - group.iloc[0]["q_start"] < 0)]

    if group.iloc[0]["s_strand"]:
        s_min = group.iloc[0]["s_start"]
        s_max = group.iloc[0]["s_end"]
        q_min = group.iloc[0]["q_start"]
        q_max = group.iloc[0]["q_end"]
        j = 1
        while s_max - s_min + 1 <= allowed_length and j < len(group):
            if group.iloc[j]["s_start"] < s_min:  # and group.iloc[j]["q_start"] < q_min:
                s_minnew = group.iloc[j]["s_start"]
                q_minnew = group.iloc[j]["q_start"]
            else:
                s_minnew = s_min
                q_minnew = q_min
            if group.iloc[j]["s_end"] > s_max:  # and group.iloc[j]["q_end"] > q_max:
                s_maxnew = group.iloc[j]["s_end"]
                q_maxnew = group.iloc[j]["q_end"]
            else:
                s_maxnew = s_max
                q_maxnew = q_max
            if s_maxnew - s_minnew + 1 <= allowed_length:  # todo: joined length will not exceed max length?
                s_min = s_minnew
                s_max = s_maxnew
                q_min = q_minnew
                q_max = q_maxnew
            j = j + 1
        s_start = s_min
        s_end = s_max
        q_start = q_min
        q_end = q_max

    else:
        s_max = group.iloc[0]["s_start"]
        s_min = group.iloc[0]["s_end"]
        q_max = group.iloc[0]["q_end"]
        q_min = group.iloc[0]["q_start"]
        j = 1
        while s_max - s_min + 1 <= allowed_length and j < len(group) - 1:
            if group.iloc[j]["s_start"] > s_max:  # and group.iloc[j]["q_start"] > q_max:
                s_maxnew = group.iloc[j]["s_start"]
                q_maxnew = group.iloc[j]["q_end"]
            else:
                s_maxnew = s_max
                q_maxnew = q_max
            if group.iloc[j]["s_end"] < s_min:  # and group.iloc[j]["q_end"] < q_min:
                s_minnew = group.iloc[j]["s_end"]
                q_minnew = group.iloc[j]["q_start"]
            else:
                s_minnew = s_min
                q_minnew = q_min
            if s_maxnew - s_minnew + 1 <= allowed_length:
                s_min = s_minnew
                s_max = s_maxnew
                q_min = q_minnew
                q_max = q_maxnew
            j = j + 1
        s_start = s_max
        s_end = s_min
        q_start = q_min
        q_end = q_max
    hits_num = len(group)
    sum_hits_alignlen = sum(group["alignment_length"])
    sum_hits_score = sum(group["bit_score"])
    q_strand = group.iloc[0]["q_strand"]
    s_strand = group.iloc[0]["s_strand"]
    return s_start, s_end, q_start, q_end, hits_num, sum_hits_alignlen, sum_hits_score, q_strand, s_strand


def join_hits(df, maxlen, qreflen):
    """
    merge hits of the same accession number
    :param df:
    :param maxlen:
    :param qreflen:
    :return:
    """
    allowed_length = qreflen * 1.05
    if allowed_length > maxlen:
        allowed_length = maxlen
    df[["alignment_length", "bit_score", "q_start", "q_end", "s_start", "s_end"]] = \
        df[["alignment_length", "bit_score", "q_start", "q_end", "s_start", "s_end"]].apply(pd.to_numeric)
    df["q_strand"] = (df["q_end"] - df["q_start"]) > 0
    df["s_strand"] = (df["s_end"] - df["s_start"]) > 0
    df["hits_num"] = 1
    df["sum_hits_alignlen"] = df["alignment_length"]
    df["sum_hits_score"] = df["bit_score"]
    groups = df.groupby(df["subject_acc.ver"])
    print("Number of hits: %d, number of accessions: %d" % (df.shape[0], len(groups)))

    len_groups = len(groups)
    arr = np.zeros((len_groups, df.shape[1]), dtype=str)
    df1 = pd.DataFrame(arr, columns=df.columns, dtype=str)

    # TODO: try apply
    for i, (name, group) in enumerate(groups):
        df1.iloc[i] = group.iloc[0]  # faster than append
        if len(group) > 1:
            s_start, s_end, q_start, q_end, hits_num, sum_hits_alignlen, sum_hits_score, q_strand, s_strand = \
                join_group(group, allowed_length)
            df1.iloc[i]["s_start"] = s_start
            df1.iloc[i]["s_end"] = s_end
            df1.iloc[i]["q_start"] = q_start
            df1.iloc[i]["q_end"] = q_end
            df1.iloc[i]["hits_num"] = hits_num
            df1.iloc[i]["sum_hits_alignlen"] = sum_hits_alignlen
            df1.iloc[i]["sum_hits_score"] = sum_hits_score
            df1.iloc[i]["q_strand"] = q_strand
            df1.iloc[i]["s_strand"] = s_strand
    return df1


def select_hits(df):
    df[["sum_hits_alignlen", "sum_hits_score"]] = df[["sum_hits_alignlen", "sum_hits_score"]].apply(pd.to_numeric)
    groups = df.groupby(df["subject_acc.ver"])
    mat = np.zeros((len(groups), df.shape[1]), dtype=str)
    df1 = pd.DataFrame(mat, columns=df.columns, dtype=str)

    i = 0
    for name, group in groups:
        if len(group) == 1:
            df1.iloc[i] = group.iloc[0]  # faster than append
        else:
            group = group.sort_values(by=["sum_hits_score", "sum_hits_alignlen"], ascending=[False, False])
            df1.iloc[i] = group.iloc[0]
        i += 1
    return df1


def hits_parse_join_select_main(wd, max_len):
    # parse xml files
    file_list = os.listdir(wd)
    for file in file_list:
        if file.split("_")[-1] == "XML.txt":
            try:
                print("Parsing %s..." % file, end="")
                time0 = datetime.now()
                parse_xml_by_re(wd, file)  # save results as HitTable and delete xml files
                time1 = datetime.now()
                print("Running time: %s Seconds" % (time1 - time0))
            except Exception as result:
                print(result)
    # join hits with the same accession in each query
    file_list = os.listdir(wd)
    hit_tables = [x.split("_")[0] for x in file_list if x.split("_")[-1] == "HitTable.txt"]
    joined_hit_tables = [x.split("_")[0] for x in file_list if x.split("_")[-1] == "joined.txt"]
    not_joined_hit_tables = set(hit_tables) - set(joined_hit_tables)
    if "blast_summary.txt" in file_list:
        sum_table = pd.read_table(Path(wd) / Path("blast_summary.txt"), sep='\t', engine='python')
    else:
        print("Cound not find blast_summary.txt")
        return
    for hit_table in not_joined_hit_tables:
        file = hit_table + "_HitTable.txt"
        try:
            print("Joining %s..." % file, end="")
            time0 = datetime.now()
            hit_table = pd.read_table(Path(wd) / Path(file), sep='\t', engine='python')
            query_ref_len = sum_table[sum_table["RID"] == file.split("_")[0]].iloc[0]["Sequence_length"]
            hit_table_merged = join_hits(hit_table, max_len, query_ref_len)
            hit_table_merged.to_csv(Path(wd) / Path(file.split("_")[0] + "_joined.txt"), index=False, sep="\t")
            time1 = datetime.now()
            print("Running time: %s Seconds" % (time1 - time0))
        except Exception as result:
            print(result)
    # select the highest bit-score hits from multiple queries
    file_list = os.listdir(wd)
    if "hits_selected.txt" not in file_list:
        time0 = datetime.now()
        print("Reading joined HitTables...", end="")
        hit_tables = []
        for file in file_list:
            if file.split("_")[-1] == "joined.txt":
                hit_table = pd.read_table(Path(wd) / Path(file), sep='\t', engine='python')
                hit_tables.append(hit_table)
        hit_tables = pd.concat(hit_tables)  # hit_tables.shape
        print("Selecting hits...", end="")
        hits_selected = select_hits(hit_tables)
        hits_selected.to_csv(Path(wd) / Path("hits_selected.txt"), index=False, sep="\t")
        time1 = datetime.now()
        print("Results saved in hits_selected.txt")
        print("Running time: %s Seconds" % (time1 - time0))

    else:
        hits_selected = pd.read_table(Path(wd) / Path("hits_selected.txt"), sep='\t', engine='python')
    return hits_selected