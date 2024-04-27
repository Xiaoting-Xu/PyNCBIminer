# *-* coding:utf-8 *-*
# @Time:2024/2/2 15:08
# @Author:Ruijing Cheng
# @File:select_new_queries.py
# @Software:PyCharm

import os
import markov_clustering as mc
import networkx as nx
import pandas as pd
import numpy as np
from Bio import AlignIO
from pathlib import Path
from datetime import datetime
from scipy.sparse import csr_matrix
from Bio import SeqIO
from tools import print_line, get_query_accession
from seq_check_download import seq_check_download_main


def cluster_queries(wd, ref_list=None):
    # todo: do this again if all selected sequences have low quality?
    df = pd.read_table(Path(wd) / Path("hits_selected.txt"), sep='\t', engine='python')
    # df.index = df["subject_acc.ver"]

    # delete reference sequences that have been used in previous search
    if ref_list is not None:
        df = df[~df["subject_acc.ver"].isin(ref_list)]

    df["qcluster"] = df["query_acc.ver"]
    sum_table = pd.read_table(Path(wd) / Path("blast_summary.txt"), sep='\t', engine='python')
    groups = df.groupby(df["query_acc.ver"])
    index_list = []
    for name, group in groups:
        try:
            print_line()
            print("%s, %d sequences" % (name, len(group)))
            if len(group) < 3:
                print("Selecting one sequence found by this query randomly ...")
                index_list.append(np.random.choice(group.index, 1, replace=False)[0])
            else:
                query_ref_len = sum_table[sum_table["Query"] == name].iloc[0]["Sequence_length"]

                # todo: remove the sequences with too high or too low bit-score?
                print("remove the sequences with too high or too low bit-score")
                high_score = group["sum_hits_score"].quantile(0.75)
                low_score = group["sum_hits_score"].quantile(0.25)
                group = group[group["sum_hits_score"] <= high_score]
                group = group[group["sum_hits_score"] >= low_score]  # 228 -> 173

                if len(group) > 1000:
                    print("Selecting 1000 sequences randomly...")
                    group = group.loc[np.random.choice(group.index, 1000, replace=False)]
                time0 = datetime.now()
                print("Extracting positions...", end="")
                positions = {i: (group.iloc[i]["q_start"], group.iloc[i]["q_end"]) for i in range(len(group))}
                time1 = datetime.now()
                print("running time: %s Seconds" % (time1 - time0))
                print("Generating network...", end="")
                network = nx.random_geometric_graph(group.shape[0], radius=0.1 * query_ref_len, pos=positions)  # get 56 clusters
                time2 = datetime.now()
                print("running time: %s Seconds" % (time2 - time1))
                print("Converting to matrix...", end="")
                matrix = nx.to_scipy_sparse_array(network)
                time3 = datetime.now()
                print("running time: %s Seconds" % (time3 - time2))
                print("Running MCL...", end="")
                result = mc.run_mcl(matrix)
                time4 = datetime.now()
                print("running time: %s Seconds" % (time4 - time3))
                print("Getting clusters...", end="")
                # todo: what if get no cluster
                clusters = mc.get_clusters(result)
                time5 = datetime.now()
                print("running time: %s Seconds" % (time5 - time4))
                print("Total running time: %s Seconds" % (time5 - time0))
                # block=True
                print("get %d clusters" % len(clusters))
                # mc.draw_graph(matrix, clusters, pos=positions, node_size=10, with_labels=False, edge_color="silver")

                # # the cluster results are reordered in group, the indices are different from the original df
                # print("Selecting one sequence randomly from each cluster...")
                # for (i, cluster) in enumerate(clusters):
                #     indices = group.iloc[list(cluster)].index
                #     # print(i, cluster)  # tuple
                #     # print(indices)
                #
                #     df.loc[indices, "qcluster"] = i
                #     index_list.append(random.choice(indices, 1, replace=False)[0])
                # df.to_csv(Path(wd) / Path("hits_selected_clusters.csv"), index=False, sep=",")

                # print("Selecting the sequence with the highest bit-score from each cluster...")
                print("Selecting the sequence with the longest align length from each cluster...")
                for (i, cluster) in enumerate(clusters):
                    print("cluster %d, %d sequences" % (i, len(cluster)))
                    indices = group.iloc[list(cluster)].index
                    df.loc[indices, "qcluster"] = i
                    # todo: select one seq with highest bit-score?
                    index_list.append(df.loc[indices].sort_values(by=["sum_hits_alignlen", "sum_hits_score"],
                                                                  ascending=[False, True]).index[0])

                df.to_csv(Path(wd) / Path("hits_selected_clusters.txt"), index=False, sep="\t")
        except Exception as result:
            print(result)
            # print("Selecting one sequence found by this query randomly ...")
            # index_list.append(np.random.choice(group.index, 1, replace=False)[0])

    df1 = df.loc[index_list]
    df1 = df1.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=False)
    df1 = df1.sort_values(by="subject_acc.ver")
    # sort table
    df1.to_csv(Path(wd) / Path("hits_clustered.txt"), index=False, sep="\t")
    # Error: The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().
    print("Selected %d sequences" % len(index_list))


def filter_seq(wd, max_length, table="hits_clustered.txt", fasta_file="hits_clustered.fasta", ref_seq_list=None):
    def get_accession(record):
        parts = record.description.split(":")
        # assert len(parts) == 2
        return parts[0]

    df = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')
    df.index = df["subject_acc.ver"]
    file_list = os.listdir(wd)

    print("Removing errorneous sequences...")
    if ("erroneous_" + fasta_file in file_list) and ("erroneous_" + table not in file_list):
        # the second condition is to avoid interruption
        seq_dict = SeqIO.to_dict(SeqIO.parse(Path(wd) / Path("erroneous_" + fasta_file), "fasta"),
                                 key_function=get_accession)
        df_err = df.loc[seq_dict.keys()].copy()
        df_err.to_csv(Path(wd) / Path("erroneous_" + table), index=False, sep="\t")
        df = df.drop(seq_dict.keys())
        df.to_csv(Path(wd) / Path(table), index=False, sep="\t")

    print("Filtering low-quality sequences...")
    df["organism"] = ""
    df["seq_len"] = 0
    index_list = []
    n = 0
    if not os.path.exists(Path(wd) / Path("hits_clustered.fasta")):
        return 0
    if os.path.getsize(Path(wd) / Path("hits_clustered.fasta")) == 0:
        return 0
    scluster = SeqIO.parse(Path(wd) / Path("hits_clustered.fasta"), "fasta")
    with open(Path(wd) / Path("hits_clustered_filtered.fasta"), "w") as fw:
        for record in scluster:
            index = record.description.split(":")[0]
            # print(index)
            # print(record.description.split("|")[1].replace("_", " "))
            df.loc[index, "organism"] = record.description.split("|")[1].replace("_", " ")
            seq = record.seq.lower()
            seq_len = seq.count("a") + seq.count("t") + seq.count("c") + seq.count("g")
            df.loc[index, "seq_len"] = seq_len
            if seq_len < 0.2 * max_length:
                print("Removed short sequence: " + index)
                n += 1
            elif seq_len < 0.995 * len(record.seq):
                print("Removed sequence with too many ambiguous bases: " + index)
                n += 1
            # # in theory, the subject will be the same as query
            elif (ref_seq_list is not None) and (record.seq in ref_seq_list):
                print("Removed duplicate sequence: " + index)
                n += 1
            else:
                # print(record.description)
                SeqIO.write(record, fw, "fasta")
                index_list.append(index)
    # df.to_csv(Path(wd) / Path("hits_clustered.txt"), index=False, sep="\t")
    df1 = df.loc[index_list]
    df1.to_csv(Path(wd) / Path("hits_clustered_filtered.txt"), index=False, sep="\t")
    print("Deleted %d sequences." % n)
    print("Filtered sequences saved in hits_clustered_filtered.fasta")
    # todo: what if no sequence left after filtering? Stop iteration.
    return len(df1)


def p_distance(wd, in_file):
    # alignment = AlignIO.read(Path(wd)/Path("msa_hits_clustered.fasta"), "fasta")
    alignment = AlignIO.read(Path(wd) / Path(in_file), "fasta")

    r = len(alignment)  # number of rows
    c = len(alignment[0])  # number of columns
    arr = np.zeros((r, r), dtype=float)
    df = pd.DataFrame(arr, dtype=float)
    for i in range(r):
        for j in range(i, r):
            different = 0
            # gap = 0
            for k in range(c):
                # if alignment[i][k] == "-" and alignment[j][k] == "-":
                #     gap += 1
                if alignment[i][k] != alignment[j][k]:
                    different += 1
            # df.iloc[i, j] = different/(c-gap)
            df.iloc[i, j] = different / c
            df.iloc[j, i] = df.iloc[i, j]
    # df.to_csv(Path(wd)/Path("Sequences_distance.csv"), index=True, sep=",")

    return df


def my_mcl(wd, df, table):
    # Nodes are considered adjacent if the distance between them is <= 0.3 units
    matrix = np.array(df)
    matrix[matrix == 0] = 1
    adjacent = (matrix <= 0.3)
    matrix[adjacent] = 1
    matrix[~adjacent] = 0
    matrix = csr_matrix(matrix)
    result = mc.run_mcl(matrix)  # run MCL with default parameters
    clusters = mc.get_clusters(result)  # 4 clusters, 11 clusters
    # mc.draw_graph(matrix, clusters, node_size=10, with_labels=True, edge_color="silver")

    df1 = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')
    df1["scluster"] = -1
    index_list = []
    for (i, cluster) in enumerate(clusters):
        df1.loc[list(cluster), "scluster"] = i
        # index_list.append(np.random.choice(cluster, 1, replace=False)[0])
        index_list.append(df1.loc[list(cluster)].sort_values(by="seq_len", ascending=False).index[0])
    df1.to_csv(Path(wd) / Path(table), index=False, sep="\t")
    # for i in range(len(clusters)):
    #     # print("select one sequence from cluster %d" % i)
    #     cluster = df1.iloc[list(clusters[i])]
    #     cluster = cluster.sort_values(by=["sum_hits_alignlen", "sum_hits_score"], ascending=[False, False])
    #     # cluster = cluster.sort_values(by="real_seq_len", ascending=False)
    #     index_list.append(cluster.index[0])

    # todo: restrict number of new reference sequences
    print("get %d clusters" % len(index_list))
    df2 = df1.iloc[index_list]
    df2 = df2.sort_values(by="subject_acc.ver")
    return df2


def cluster_sequences_main(wd, fasta_file=r"hits_clustered_filtered.fasta"):
    print("Clustering sequences...")
    # todo: if hits_clustered_filtered.fasta only has two sequence
    if os.path.getsize(Path(wd) / Path(fasta_file)) > 0:
        n_seq = 0
        for record in SeqIO.parse(Path(wd) / Path(fasta_file), "fasta"):
            n_seq += 1
        if n_seq < 5:
            seq_clustered = pd.read_table(Path(wd) / Path("hits_clustered_filtered.txt"), sep="\t")
            seq_clustered.to_csv(Path(wd) / Path("sequences_clustered.txt"), index=False, sep="\t")
            print("Only %d sequences left after filtering. No need to do MCL step2." % n_seq)
            return seq_clustered
        mafft_cmd = "mafft --localpair --maxiterate 1000 %s > %s" % (
            Path(wd) / Path(fasta_file), Path(wd) / Path("msa_" + fasta_file))
        os.system(mafft_cmd)

    if os.path.getsize(Path(wd) / Path("msa_" + fasta_file)) > 0:
        seq_distance = p_distance(wd, "msa_" + fasta_file)
        # todo: what if no return
        seq_distance.to_csv(Path(wd) / Path(Path("msa_" + fasta_file).stem + "_distance.txt"), index=True, sep="\t")
        table = r"hits_clustered_filtered.txt"
        seq_clustered = my_mcl(wd, seq_distance, table)
        seq_clustered.to_csv(Path(wd) / Path("sequences_clustered.txt"), index=False, sep="\t")
        return seq_clustered


def cluster_sequences(wd, fasta_file=r"hits_clustered_filtered.fasta"):
    print("Clustering sequences...")
    # todo: if hits_clustered_filtered.fasta only has two sequence
    if os.path.getsize(Path(wd) / Path(fasta_file)) > 0:
        n_seq = 0
        for record in SeqIO.parse(Path(wd) / Path(fasta_file), "fasta"):
            n_seq += 1
        if n_seq < 5:
            seq_clustered = pd.read_table(Path(wd) / Path("hits_clustered_filtered.txt"), sep="\t")
            seq_clustered.to_csv(Path(wd) / Path("sequences_clustered.txt"), index=False, sep="\t")
            print("Only %d sequences left after filtering. No need to do MCL step2." % n_seq)
            return seq_clustered
        mafft_cmd = "mafft --localpair --maxiterate 1000 %s > %s" % (
            Path(wd) / Path(fasta_file), Path(wd) / Path("msa_" + fasta_file))
        os.system(mafft_cmd)
    else:
        return None

    if os.path.getsize(Path(wd) / Path("msa_" + fasta_file)) > 0:
        seq_distance = p_distance(wd, "msa_" + fasta_file)
        # todo: what if no return
        seq_distance.to_csv(Path(wd) / Path(Path("msa_" + fasta_file).stem + "_distance.txt"), index=True, sep="\t")
        table = r"hits_clustered_filtered.txt"
        # seq_clustered = my_mcl(wd, seq_distance, table)  # my_mcl(wd, df, table)
        try:
            matrix = np.array(seq_distance)
            matrix[matrix == 0] = 1
            adjacent = (matrix <= 0.3)
            matrix[adjacent] = 1
            matrix[~adjacent] = 0
            matrix = csr_matrix(matrix)
            result = mc.run_mcl(matrix)  # run MCL with default parameters
            clusters = mc.get_clusters(result)  # 4 clusters, 11 clusters
            # mc.draw_graph(matrix, clusters, node_size=10, with_labels=True, edge_color="silver")

            df1 = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')
            df1["scluster"] = -1
            index_list = []
            for (i, cluster) in enumerate(clusters):
                df1.loc[list(cluster), "scluster"] = i
                # index_list.append(np.random.choice(cluster, 1, replace=False)[0])
                index_list.append(df1.loc[list(cluster)].sort_values(by="seq_len", ascending=False).index[0])
            df1.to_csv(Path(wd) / Path(table), index=False, sep="\t")

            # todo: if only get one cluster? or if clustering failed?
            print("get %d clusters" % len(index_list))
            df2 = df1.iloc[index_list]
            df2 = df2.sort_values(by="subject_acc.ver")
            df2.to_csv(Path(wd) / Path("sequences_clustered.txt"), index=False, sep="\t")
            return df2
        except Exception as result:
            print("MCL step2 clustering according to sequence distance failed.")
            print("Select one sequence randomly from MCL step1.")
            df1 = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')
            df2 = df1.loc[np.random.choice(df1.index, 1, replace=False)[0]].to_frame().T
            df2.to_csv(Path(wd) / Path("sequences_clustered.txt"), index=False, sep="\t")
            print(result)
            return df2


def select_new_queries(wd, tmp_wd, blast_round, ref_number):
    seq_clustered = pd.read_table(Path(tmp_wd) / Path(r"sequences_clustered.txt"), sep='\t', engine='python')
    seq_clustered.index = seq_clustered["subject_acc.ver"]
    # file_list = os.listdir(Path(wd) / Path("parameters"))
    """
    if "all_new_queries_info.txt" in file_list:
        all_new_queries = pd.read_table(Path(wd) / Path("parameters") / Path(r"all_new_queries_info.txt"), sep='\t',
                                        engine='python')
    #     all_new_queries.index = all_new_queries["subject_acc.ver"]
    #     indices = set(seq_clustered.index) - set(all_new_queries.index)
    #     if len(indices) < 1:
    #         print("Cannot find more new references. Stop iteration.")
    #         return 0
    #     seq_clustered = seq_clustered.loc[indices]
    else:
        all_new_queries = None
    """

    if seq_clustered.shape[0] > ref_number:
        print("Selecting %d sequences randomly..." % ref_number)
        index_list = np.random.choice(range(seq_clustered.shape[0]), ref_number, replace=False)
        seq_clustered = seq_clustered.iloc[index_list]
        # print("Selecting %d longest sequences..." % ref_number)
        # seq_clustered = seq_clustered.sort_values(by="seq_len", ascending=False)
        # seq_clustered = seq_clustered.iloc[range(ref_number)]

    seq_clustered["blast_round"] = blast_round
    # seq_clustered = seq_clustered.sort_values(by="subject_acc.ver")
    # ValueError: 'subject_acc.ver' is both an index level and a column label, which is ambiguous.
    """
    if all_queries is None:
        all_queries = seq_clustered
    else:
        all_queries = all_queries.append(seq_clustered)
    """

    # def get_accession(record):
    #     parts = record.description.split(":")
    #     # assert len(parts) == 2
    #     return parts[0]

    fw = open(Path(tmp_wd) / Path("new_queries.fasta"), "w")
    fw.close()
    fasta_file = r"hits_clustered_filtered.fasta"
    seq_dict = SeqIO.to_dict(SeqIO.parse(Path(tmp_wd) / Path(fasta_file), "fasta"), key_function=get_query_accession)
    for index in seq_clustered.index:
        record = seq_dict[index]
        with open(Path(tmp_wd) / Path("new_queries.fasta"), "a") as fw:
            SeqIO.write(record, fw, "fasta")
    print("Selected %d new references. " % seq_clustered.shape[0])

    seq_clustered.to_csv(Path(tmp_wd) / Path("new_queries_info.txt"), index=False,
                         sep="\t")  # write info table after new seq added to fas file
    # all_new_queries.to_csv(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), index=False, sep="\t")

    return seq_clustered.shape[0]


def select_new_queries_main(wd, tmp_wd, key_annotations, exclude_sources, entrez_email, max_length, blast_round, ref_number):
    file_list = os.listdir(tmp_wd)
    # select new reference sequences from different clusters
    if "hits_clustered.txt" not in file_list:
        print("Selecting new references...")
        file_list = os.listdir(Path(wd) / Path("parameters"))
        # if os.path.exists(Path(wd) / Path("parameters") / Path("all_new_queries_info.txt")):
        #     all_new_queries = pd.read_table(Path(wd) / Path("parameters") / Path(r"all_new_queries_info.txt"),
        #                                     sep='\t', engine='python')
        #     all_new_queries_list = all_new_queries["subject_acc.ver"]
        #     cluster_queries(wd=tmp_wd, ref_list=all_new_queries_list)
        if os.path.exists(Path(wd) / Path("parameters") / Path("all_queries_info.txt")):
            all_queries = pd.read_table(Path(wd) / Path("parameters") / Path(r"all_queries_info.txt"),
                                        sep='\t', engine='python')
            all_queries_list = all_queries["ID"]
            cluster_queries(wd=tmp_wd, ref_list=all_queries_list)
            # todo: what to do if no more new reference sequnces could be found
        else:
            cluster_queries(wd=tmp_wd)
    else:
        qcluster_table = pd.read_table(Path(tmp_wd) / Path("hits_clustered.txt"), sep="\t", engine="python")
        print("Get %d clusters according to query start and query end." % qcluster_table.shape[0])
        del qcluster_table
    # todo: hits_clustered.txt main contain repeated rows???

    if "hits_clustered_filtered.txt" not in file_list:
        seq_check_download_main(wd=tmp_wd, acc_file=r"hits_clustered.txt",
                                out_file=r"hits_clustered.fasta",
                                key_annotations=key_annotations, exclude_sources=exclude_sources,
                                entrez_email=entrez_email)
        if os.path.exists(Path(wd) / Path("parameters") / Path("ref_seq")):
            ref_file_list = os.listdir(Path(wd) / Path("parameters") / Path("ref_seq"))
            ref_seq_list = []
            for file in ref_file_list:
                for record in SeqIO.parse(Path(wd) / Path("parameters") / Path("ref_seq") / Path(file), "fasta"):
                    ref_seq_list.append(record.seq)
            n_hits_clustered_filtered = filter_seq(wd=tmp_wd, max_length=max_length, ref_seq_list=ref_seq_list)
            # print("Number of hits_clustered_filtered %d" % n_hits_clustered_filtered)
            # todo: sometimes the clustering will fail when sequences number is too small.
            if n_hits_clustered_filtered < 1:
                print("Cannot find more new reference, stop iteration. ")
                return None
    # if hits_clustered.fasta only contain no more than 5 sequences, then no need to do MCL step2
    if "sequences_clustered.txt" not in file_list:
        scluster_table = cluster_sequences(wd=tmp_wd)
    else:
        scluster_table = pd.read_table(Path(tmp_wd) / Path("sequences_clustered.txt"), sep="\t", engine="python")

    if scluster_table is None:
        print("No sequence passed filtering.")
        print("Cannot find more new reference, stop iteration. ")
        return None

    elif scluster_table.shape[0] > 0:
        if "scluster" in scluster_table.columns:
            print("Get %d clusters according to sequence distance." % scluster_table.shape[0])
        else:
            print("Using all of the %d sequences selected in MCL step1." % scluster_table.shape[0])
        del scluster_table

    else:
        print("Cannot find more new reference, stop iteration. ")
        del scluster_table
        return None

    blast_round += 1
    if "new_queries_info.txt" not in file_list:
        new_quereis_num = select_new_queries(wd, tmp_wd, blast_round, ref_number)
        # todo: new ref seqs need to be more than 2???
        if new_quereis_num < 1:
            print("Cannot find more new reference. Stop iteration.")
            return None
    else:
        new_queries = pd.read_table(Path(tmp_wd) / Path(r"new_queries_info.txt"), sep='\t', engine='python')
        new_quereis_num = new_queries.shape[0]
        print("Selected %d new queries." % new_quereis_num)

    if not os.path.exists(Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_%d.fasta" % blast_round)):
        os.system("copy %s %s" % (Path(tmp_wd) / Path("new_queries.fasta"),
                                  Path(wd) / Path("parameters") / Path("ref_seq") / Path(
                                      "queries_%d.fasta" % blast_round)))

    # extend hits After BLAST iteration.
    return new_quereis_num
