# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:44
# @Author:Ruijing Cheng
# @File:iterated_blast.py
# @Software:PyCharm

import pandas as pd
import numpy as np
from pathlib import Path
import os
from math import ceil
from Bio import SeqIO

from tools import print_line, get_query_accession
from my_entrez import format_entrez_query
from seq_check_download import check_annotation, seq_check_download_main
from blast_put_get import blast_put_get_main
from hits_parse_join_select import hits_parse_join_select_main
from select_new_queries import select_new_queries_main
from blast_results_extend import blast_results_extend_main


def combine_iterated_blast(wd, blast_round, tmp_df, key_annotations, exclude_sources):
    if blast_round == 1:
        tmp_df.to_csv(Path(wd) / Path("results") / Path("blast_results.txt"), index=False, sep="\t")
        blast_count_new = tmp_df.shape[0]
        blast_count = blast_count_new
        blast_count_new1 = blast_count_new
        print("Results saved in blast_results.txt")
    else:
        print("Reading results from multiple BLAST...")
        df = pd.read_table(Path(wd) / Path("results") / Path("blast_results.txt"), sep="\t", engine="python")
        # if blast_round in df["Source"]: # always True, fail to combine new sequences
        if blast_round in set(df["Source"].values):
            print("Results of BLAST round %d already in blast_results.txt" % blast_round)
            new_df = df[df["Source"] == blast_round]
            blast_count_new = len(new_df)
            blast_count = len(df)
        else:
            # df.index = df["subject_acc.ver"]
            # tmp_df.index = tmp_df["subject_acc.ver"]
            acc_set = set(df["subject_acc.ver"])
            tmp_acc_set = set(tmp_df["subject_acc.ver"])
            # common_acc_set = acc_set & tmp_acc_set
            new_acc_set = tmp_acc_set - acc_set
            blast_count_new = len(new_acc_set)
            blast_count = len(acc_set) + blast_count_new
            # if blast_count_new == 0:
            #     source_list = df["Source"].value_counts()
            #     if blast_round in source_list.index:
            #         blast_count_new = source_list.loc[blast_round]

            # # compare sequences with the same accession, kep the one with highest bit-score
            # for acc in common_acc_set:
            #     if df.loc[acc]["sum_hits_score"] < tmp_df.loc[acc]["sum_hits_score"]:
            #         df.loc[acc] = tmp_df.loc[acc]

            new_df = tmp_df[tmp_df["subject_acc.ver"].isin(new_acc_set)]

            df = pd.concat([df, new_df])
            df.to_csv(Path(wd) / Path("results") / Path("blast_results.txt"), index=False, sep="\t")
        print("Combined results saved in BLAST_results.txt")
        # check definition of all new sequences

        print("Checking definition of all new sequences...")
        blast_count_new1 = blast_count_new
        for index in new_df.index:
            hit_def = new_df.loc[index, "hit_def"]
            hit_def = [hit_def.replace(" ", "")]
            key_annotations = [x.replace(" ", "").lower() for x in key_annotations]
            exclude_sources = [x.replace(" ", "").lower() for x in exclude_sources]
            if not check_annotation(hit_def, key_annotations, exclude_sources):
                blast_count_new1 -= 1

    return blast_count_new, blast_count_new1, blast_count


def iterated_blast_main(wd, organisms, count,
                        expect, gapcosts, word_size,
                        nucl_reward, nucl_penalty, max_length,
                        key_annotations, exclude_sources, ref_number, date_from, date_to, entrez_email,
                        blast_round=1):
    blast_round = 1  # to correct error in combining blast results.
    print("Start BLAST iteration...")
    entrez_query = format_entrez_query(organisms=organisms, date_from=date_from, date_to=date_to)
    alignments = ceil(count * 1.05)
    last_new = 9999
    # read previous BLAST results if blast_results.txt exists in results folder
    if os.path.exists(Path(wd) / Path("results") / Path("blast_results.txt")):
        print("Reading previous BLAST results...")
        blast_results = pd.read_table(Path(wd) / Path("results") / Path("blast_results.txt"), sep='\t', engine='python')
        source_list = blast_results["Source"].value_counts()
        total_seq_num = 0
        for i in range(blast_round - 1):
            if i + 1 in source_list.index:
                seq_num = source_list.loc[i + 1]
            else:
                seq_num = 0
            total_seq_num += seq_num
            print("BLAST round %d: %d new sequences, %d sequences in total." % (i + 1, seq_num, total_seq_num))
            last_new = seq_num
        if blast_round > 1:
            print("Find %d new sequences in the last round." % last_new)

        if blast_round > 2:
            last_df = blast_results[blast_results["Source"] == blast_round - 1]
            for index in last_df.index:
                hit_def = last_df.loc[index, "hit_def"]
                hit_def = [hit_def.lower().replace(" ", "")]
                key_annotations = [x.replace(" ", "").lower() for x in key_annotations]
                exclude_sources = [x.replace(" ", "").lower() for x in exclude_sources]
                if not check_annotation(hit_def, key_annotations, exclude_sources):
                    last_new -= 1
            print("%d sequences passed definition checking." % last_new)
            del last_df
        del blast_results

    while True:
        print_line("*")
        print("BLAST round %d" % blast_round)
        if blast_round == 1:
            if not os.path.exists(Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta")):
                description_list = []
                seq_list = []
                m = 0
                n = 0
                # inital_queries.fasta may contain duplicate sequences, queries_1.fasta only saves unique sequences.
                with open(Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta"), "w") as fw:
                    for record in SeqIO.parse(Path(wd) / Path("parameters") / Path("initial_queries.fasta"), "fasta"):
                        m += 1
                        if record.description in description_list:
                            pass
                        elif record.seq in seq_list:
                            pass
                        else:
                            description_list.append(record.description)
                            seq_list.append(record.seq)
                            SeqIO.write(record, fw, "fasta")
                            n += 1
                if m > n:
                    print("Removed %d duplicate queries" % (m - n))
                if n < 2:
                    print("At least 2 initial queries are required")
                    return
            query_path = Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta")

            """
            query_path = Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta")
            if not os.path.exists(Path(query_path)) or os.path.getsize(Path(query_path)) == 0:
                print("Aligning initial queries...")
                # os.system(
                #     "mafft --auto %s > %s" % (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta"),
                #                               Path(wd) / Path("parameters") / Path("ref_msa") / Path(
                #                                   "msa_queries_1.fasta")))

                # using L-INS-i (Very slow; recommended for <200 sequences with one conserved domain and long gaps)

                os.system(
                    "mafft --localpair --maxiterate 1000 %s > %s" %
                    (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_1.fasta"), query_path))
            if os.path.getsize(query_path) == 0:
                print(
                    "Failed to align inital queries, "
                    "please check if MAFFT is correctly installed.")
                return
            """
        else:
            query_path = Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_%d.fasta" % blast_round)
            """
            query_path = Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_%d.fasta" % blast_round)
            if not os.path.exists(Path(query_path)) or os.path.getsize(Path(query_path)) == 0:
                print("Aligning new queries...")
                if os.path.getsize(Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta")) == 0:
                    print("The file size of msa_queries_1.fasta is zero, please check the initial queries.")
                    return
                align_new_query(wd, blast_round)
            if os.path.getsize(Path(query_path)) == 0:
                print("Failed to align new queries, please check the queries_%d.fasta file" % blast_round)
                return
            """

        # save queries information
        queries = SeqIO.to_dict(SeqIO.parse(Path(query_path), "fasta"), key_function=get_query_accession)
        column_list = ["ID", "Description", "Sequence", "Sequence_length", "blast_round"]
        sum_mat = np.zeros((len(queries), len(column_list)), dtype=str)
        sum_table = pd.DataFrame(sum_mat, columns=column_list, dtype=str)

        for (i, key) in enumerate(queries.keys()):
            sum_table.loc[i, "ID"] = key
            sum_table.loc[i, "Description"] = queries[key].description
            sum_table.loc[i, "Sequence"] = str(queries[key].seq.upper())
            sum_table.loc[i, "Sequence_length"] = len(sum_table.loc[i, "Sequence"])
            sum_table.loc[i, "blast_round"] = blast_round
        if not os.path.exists(Path(wd) / Path("parameters") / Path("all_queries_info.txt")):
            sum_table.to_csv(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), index=False, sep="\t")
            print("Queries information is saved in all_queries_info.txt.")
        else:
            all_queries = pd.read_table(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), sep="\t")
            if blast_round not in set(all_queries["blast_round"]):
                all_queries = pd.concat([all_queries, sum_table])
                all_queries.to_csv(Path(wd) / Path("parameters") / Path("all_queries_info.txt"), index=False, sep="\t")
                print("Queries information is updated in all_queries_info.txt.")
            else:
                print("Queries information is already saved in all_queries_info.txt.")

        folder = "BLAST_%d" % blast_round
        tmp_wd = Path(wd) / Path("tmp_files") / Path(folder)
        if not os.path.exists(tmp_wd):
            os.makedirs(tmp_wd)

        # put BLAST request and get BLAST results
        blast_put_get_main(wd=tmp_wd, queries_path=query_path, entrez_query=entrez_query,
                           alignments=alignments, expect=expect, gapcosts=gapcosts, word_size=word_size,
                           nucl_reward=nucl_reward, nucl_penalty=nucl_penalty)
        # join and extend hits
        file_list = os.listdir(tmp_wd)
        # todo: do not extend in each round of BLAST.
        # tmp_results = hits_join_extend_main(wd=tmp_wd, max_len=max_length)
        tmp_results = hits_parse_join_select_main(wd=tmp_wd, max_len=max_length)
        tmp_results["Source"] = blast_round
        # combine results of iterated BLAST
        blast_count_new, blast_count_new1, blast_count = combine_iterated_blast(wd=wd, blast_round=blast_round,
                                                                                tmp_df=tmp_results,
                                                                                key_annotations=key_annotations,
                                                                                exclude_sources=exclude_sources)
        print("Find %d new sequences in round %d of BLAST." % (blast_count_new, blast_round))
        print("%d sequences passed definition checking." % blast_count_new1)
        print("Find %d sequences in total." % blast_count)

        # check if the iteration should stop
        # if blast_count > count * 1.1:
        if blast_count > count * 1.2:
            print("Total sequences number is greater than Entrez search results number. Stop BLAST iteration.")
            break
        if blast_count_new1 < 3 and last_new < 3:
            print("Cannot find more than 2 correct new sequences. Stop BLAST iteration.")
            break
        last_new = blast_count_new1

        # select new reference sequences
        new_quereis_num = select_new_queries_main(wd, tmp_wd, key_annotations, exclude_sources, entrez_email, max_length, blast_round,ref_number)
        if new_quereis_num is None:
            break
        blast_round += 1

        """
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
                    break
        # if hits_clustered.fasta only contain no more than 5 sequences, then no need to do MCL step2
        if "sequences_clustered.txt" not in file_list:
            scluster_table = cluster_sequences(wd=tmp_wd)
        else:
            scluster_table = pd.read_table(Path(tmp_wd) / Path("sequences_clustered.txt"), sep="\t", engine="python")

        if scluster_table is None:
            print("No sequence passed filtering.")
            print("Cannot find more new reference, stop iteration. ")
            break

        elif scluster_table.shape[0] > 0:
            if "scluster" in scluster_table.columns:
                print("Get %d clusters according to sequence distance." % scluster_table.shape[0])
            else:
                print("Using all of the %d sequences selected in MCL step1." % scluster_table.shape[0])
            del scluster_table

        else:
            print("Cannot find more new reference, stop iteration. ")
            del scluster_table
            break

        blast_round += 1
        if "new_queries_info.txt" not in file_list:
            new_quereis_num = select_new_queries(wd, tmp_wd, blast_round, ref_number)
            # todo: new ref seqs need to be more than 2???
            if new_quereis_num < 1:
                print("Cannot find more new reference. Stop iteration.")
                break
        else:
            new_queries = pd.read_table(Path(tmp_wd) / Path(r"new_queries_info.txt"), sep='\t', engine='python')
            print("Selected %d new queries." % new_queries.shape[0])

        if not os.path.exists(Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_%d.fasta" % blast_round)):
            os.system("copy %s %s" % (Path(tmp_wd) / Path("new_queries.fasta"),
                                      Path(wd) / Path("parameters") / Path("ref_seq") / Path(
                                          "queries_%d.fasta" % blast_round)))

    # extend hits After BLAST iteration.
    ref_msa_file = add_all_queries2(wd)
    blast_results_extend_main(wd, max_length, ref_msa_file)
    """

    blast_results_extend_main(wd, max_length)

    seq_check_download_main(wd=Path(wd) / Path("results"), acc_file=r"blast_results.txt",
                            out_file=r"blast_results_checked.fasta",
                            key_annotations=key_annotations, exclude_sources=exclude_sources, entrez_email=entrez_email,
                            extend=True)
    print("Stop thread.")