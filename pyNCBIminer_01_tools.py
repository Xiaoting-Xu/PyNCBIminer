# *-* coding:utf-8 *-*
# @Time:2022/5/6 15:15
# @Author:Ruijing Cheng
# @File:pyNCBIminer_01_tools.py
# @Software:PyCharm


import os
import re
import time
import urllib.error
import markov_clustering as mc
# Backend QtAgg is interactive backend. Turning interactive mode on.
import networkx as nx
from math import floor, ceil
import pandas as pd
import numpy as np
from Bio import AlignIO, Entrez  # , SeqIO
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request, urlretrieve
from pathlib import Path
from datetime import datetime
from scipy.sparse import csr_matrix
import func_timeout.exceptions
from func_timeout import func_set_timeout
# from fake_useragent import UserAgent
import shutil
import sys
from Bio import SeqIO, SeqRecord
import zipfile


def print_line(character="#"):
    print(character * 50)


def schedule(blocknum, blocksize, totalsize):
    per = 100.0 * blocknum * blocksize / totalsize
    if per > 100:
        per = 100
    if per % 10 < (100.0 * blocksize / totalsize):
        # print("%d %d %d" % (blocknum, blocksize, totalsize))
        print("%.2f.%%" % per, end=" ")


def install_mafft():
    """
    download mafft to the same directory as pyncbiminer and extract the zip file
    https://mafft.cbrc.jp/alignment/software/mafft-7.520-win64-signed.zip
    path to trimal after extraction: ./mafft/mafft-win
    :return:
    """
    # todo: this url needs to be updated or dynamically maintained.
    mafft_url = r"https://mafft.cbrc.jp/alignment/software/mafft-7.520-win64-signed.zip"
    root_path = os.getcwd()
    mafft_dir = Path(root_path) / Path(r"./mafft")
    if not os.path.exists(mafft_dir):
        os.makedirs(mafft_dir)
    file_path = Path(mafft_dir) / Path(r"mafft-7.520-win64-signed.zip")
    print("Downloading MAFFT...")
    file_path, _ = urlretrieve(mafft_url, file_path, schedule)
    print("Extracting files...")
    with zipfile.ZipFile(file_path, "r") as zip:
        zip.extractall(mafft_dir)
    print("Installation finished.")


def install_trimal():
    """
    download trimal to the same directory as pyncbiminer and extract the zip file
    http://trimal.cgenomics.org/_media/trimal.v1.2rev59.zip
    path to trimal after extraction: ./trimal/trimAl/bin
    :return:
    """
    # todo: this url needs to be updated or dynamically maintained.
    trimal_url = r"http://trimal.cgenomics.org/_media/trimal.v1.2rev59.zip"
    root_path = os.getcwd()
    trimal_dir = Path(root_path) / Path(r"./trimal")
    if not os.path.exists(trimal_dir):
        os.makedirs(trimal_dir)
    file_path = Path(trimal_dir) / Path(r"trimal.v1.2rev59.zip")
    print("Downloading trimAl...")
    file_path, _ = urlretrieve(trimal_url, file_path, schedule)

    print("Extracting files...")
    with zipfile.ZipFile(file_path, "r") as zip:
        zip.extractall(trimal_dir)
    print("Installation finished.")


def format_entrez_query(organisms, entrez_qualifier="", date_from="", date_to=""):
    """
    connects organism, entrez qualifier and publication date
    :param organisms: a list of organisms
    :param entrez_qualifier: the specific entrez search range
    :param date_from: start of publication date span
    :param date_to: end of publication date span
    :return: formated entrez query
    """
    if len(organisms) > 0:
        entrez_query = '%s[Organism]' % organisms[0]
        if len(organisms) > 1:
            for i in range(1, len(organisms)):
                entrez_query = entrez_query + ' OR %s[Organism]' % organisms[i]
        if entrez_qualifier.upper().startswith("NOT") or entrez_qualifier.upper().startswith("OR"):
            entrez_query = entrez_query + entrez_qualifier
        else:
            entrez_query = ('(%s) AND (%s)' % (entrez_query, entrez_qualifier))
    else:
        entrez_query = entrez_qualifier
    if len(date_to) > 0:
        if len(date_from) > 0:
            date_span = '"%s"[Publication Date] : "%s"[Publication Date]' % (date_from, date_to)
        else:
            date_span = '("1900"[Publication Date] : "%s"[Publication Date])' % date_to
        if len(entrez_query) > 0:
            entrez_query += ' AND ' + date_span
        else:
            entrez_query += date_span
    return entrez_query


@func_set_timeout(60)
def entrez_count(entrez_query, entrez_email):
    """
    send entrez query to NCBI and get the entrez search results count
    :param entrez_query: the entrez query
    :param entrez_email: the user's email address
    :return: entrez search results count
    """
    if len(entrez_email) == 0:
        print("Warning: email address is not specified.")
        print("To make use of NCBI's E-utilities, NCBI requires you to specify your email address with each request. ")
        print("In case of excessive usage of the E-utilities, "
              "NCBI will attempt to contact a user at the email address provided "
              "before blocking access to the E-utilities.")
    Entrez.email = entrez_email
    handle = Entrez.esearch(db="nucleotide", term=entrez_query)
    record = Entrez.read(handle)
    print("Entrez search results count: %s" % record["Count"])
    return int(record["Count"])


def get_query_accession(record):
    parts = record.description.split(" ")[0].split("|")[0].split(":")
    # assert len(parts) == 2
    return parts[0]


def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is probably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = handle.read().decode()
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i + len("RID ="): j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i + len("RTOE ="): j].strip()

    if not rid and not rtoe:
        # Can we reliably extract the error message from the HTML page?
        # e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        # or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        # This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i + len('<div class="error msInf">'):].strip()
            msg = msg.split("</div>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i + len('<p class="error">'):].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        # Generic search based on the way the error messages start:
        i = s.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = s[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError("Error message from NCBI: %s" % msg)
        # We didn't recognise the error layout :(
        # print(s)
        raise ValueError(
            "No RID and no RTOE found in the 'please wait' page, "
            "there was probably an error in your request but we "
            "could not extract a helpful error message."
        )
    elif not rid:
        # Can this happen?
        raise ValueError(
            "No RID found in the 'please wait' page. (although RTOE = %r)" % rtoe
        )
    elif not rtoe:
        # Can this happen?
        raise ValueError(
            "No RTOE found in the 'please wait' page. (although RID = %r)" % rid
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            "A non-integer RTOE found in the 'please wait' page, %r" % rtoe
        ) from None


@func_set_timeout(60)
def put_blast_requests(url_base, message, header):
    request = Request(url_base, message, headers=header)
    handle = urlopen(request)
    print("Parsing BLAST ref page...")
    rid, rtoe = _parse_qblast_ref_page(handle)  # get rid and rtoe
    return rid, rtoe


def put_blast(wd, queries_path=None,
              program="blastn",
              database="nt",
              url_base="https://blast.ncbi.nlm.nih.gov/Blast.cgi",
              entrez_query="(none)",
              expect=10.0,
              gapcosts="2 1",
              alignments=1000,
              word_size=7,
              nucl_reward=1,
              nucl_penalty=-1,
              table="blast_summary.txt"):
    """
    Format the "Put" command, send search requests to NCBI, get RID and RTOE
    RTOE is probably 'Request Time of Execution' and RID would be 'Request Identifier'
    """
    if os.path.exists(Path(wd) / Path(table)):
        sum_table = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')
        sum_table.loc[sum_table["RID"].isna(), "RID"] = ""
    else:
        # queries = AlignIO.read(Path(queries_path), "fasta")
        # sum_mat = np.zeros((len(queries), 14), dtype=str)
        # sum_table = pd.DataFrame(sum_mat, columns=["ID", "Description", "Sequence", "Sequence_length",
        #                                            "Missing_left", "Missing_right", "Message_put", "Time_put", "RID",
        #                                            "RTOE", "Message_get", "Time_get", "Status", "Query"], dtype=str)

        queries = SeqIO.to_dict(SeqIO.parse(Path(queries_path), "fasta"), key_function=get_query_accession)
        column_list = ["ID", "Description", "Sequence", "Sequence_length",
                       "Message_put", "Time_put", "RID",
                       "RTOE", "Message_get", "Time_get", "Status", "Query"]
        sum_mat = np.zeros((len(queries), len(column_list)), dtype=str)
        sum_table = pd.DataFrame(sum_mat, columns=column_list, dtype=str)

        for (i, key) in enumerate(queries.keys()):
            sum_table.loc[i, "ID"] = key
            sum_table.loc[i, "Description"] = queries[key].description
            sum_table.loc[i, "Sequence"] = str(queries[key].seq.upper())
            sum_table.loc[i, "Sequence_length"] = len(sum_table.loc[i, "Sequence"])
            # print(key)
            # print(queries[key].description)
            # print(str(queries[key].seq.upper()))

        # for i in range(len(queries)):
        #     sum_table.loc[i, "ID"] = queries[i].id
        #     sum_table.loc[i, "Description"] = queries[i].description
        #     seq = str(queries[i].seq.upper())
        #     sum_table.loc[i, "Sequence"] = seq.replace("-", "")
        #     sum_table.loc[i, "Sequence_length"] = len(sum_table.loc[i, "Sequence"])
        # left = [seq.find("A"), seq.find("T"), seq.find("C"), seq.find("G")]
        # right = [seq.rfind("A"), seq.rfind("T"), seq.rfind("C"), seq.rfind("G")]
        # sum_table.loc[i, "Missing_left"] = min(left)
        # sum_table.loc[i, "Missing_right"] = len(seq) - max(right) - 1

        # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
        # new website: https://ncbi.github.io/blast-cloud/dev/api.html (2023/06/12)

        for index in sum_table.index:
            sequence = sum_table.loc[index, "Sequence"]
            parameters = [
                ("HITLIST_SIZE", alignments),
                # ("ALIGNMENTS", alignments),
                ("DATABASE", database),
                ("ENTREZ_QUERY", entrez_query),
                ("EXPECT", expect),
                ("GAPCOSTS", gapcosts),
                ("NUCL_PENALTY", nucl_penalty),
                ("NUCL_REWARD", nucl_reward),
                ("PROGRAM", program),
                ("QUERY", sequence),
                ("WORD_SIZE", word_size),
                ("CMD", "Put")
            ]
            query = [x for x in parameters if x[1] is not None]
            message = urlencode(query)
            sum_table.loc[index, "Message_put"] = message
        sum_table.to_csv(Path(wd) / Path(table), index=False, sep="\t")  # save parameters
        print("Information of initial query saved in %s" % table)
    # queries_info = sum_table[["ID", "Description", "Sequence", "Sequence_length"]].copy()
    # queries_info["blast_round"] = blast_round

    # Note the NCBI do not currently impose a rate limit here,
    # other than the request not to make say 50 queries at once using multiple threads.
    # ua = UserAgent()
    # header = {"User-Agent": ua.random}
    indices = sum_table[sum_table["RID"] == ""].index
    while len(indices) > 0:
        for n in range(len(indices)):
            try:
                index = indices[n]
                message = sum_table.loc[index, "Message_put"]
                message = message.encode()
                print("Try submitting Query %s..." % sum_table.loc[index, "ID"])
                # request = Request(url_base, message, {"User-Agent": "BiopythonClient"})
                # # request = Request(url_base, message, headers=header)
                # handle = urlopen(request) # get rid and rtoe
                # rid, rtoe = _parse_qblast_ref_page(handle)
                rid, rtoe = put_blast_requests(url_base, message, {"User-Agent": "BiopythonClient"})
                sum_table.loc[index, "Time_put"] = time.time()
                sum_table.loc[index, "RID"] = rid
                sum_table.loc[index, "RTOE"] = rtoe
                sum_table.to_csv(Path(wd) / Path(table), index=False, sep="\t")
                print("Query %s submitted, RID = %s, RTOE = %s" % (sum_table.loc[index, "ID"], rid, rtoe))
            except func_timeout.exceptions.FunctionTimedOut:
                print("Time out, try submitting Query %s again..." % sum_table.loc[index, "ID"])
            except Exception as result:
                print(result)
        indices = sum_table[sum_table["RID"] == ""].index
    print("All queries submitted!")


@func_set_timeout(600)
def get_blast_results(url_base, message, header):
    request = Request(url_base, message, headers=header)
    handle = urlopen(request)  # time-consuming
    print("Decoding results...")
    results = handle.read().decode()  # tabular is fast to decode, but other formats take more time
    return results


def get_blast(wd,
              url_base="https://blast.ncbi.nlm.nih.gov/Blast.cgi",
              alignments=1000,
              format_type="XML",
              table="blast_summary.txt"
              ):
    """
    Format the "Get" command, get the formatted results from put_blast
    :param wd: the directory where blast_summary.txt file is saved
    :param format_type: report type, XML is the default.
    :param alignments:
    :param url_base:
    :param table:
    :return:
    """

    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007
    # new website: https://ncbi.github.io/blast-cloud/dev/api.html (2023/06/12)

    sum_table = pd.read_table(Path(wd) / Path(table), sep='\t', engine='python')

    for index in sum_table.index:
        rid = sum_table.loc[index, "RID"]
        parameters = [
            ("HITLIST_SIZE", alignments),
            # ("ALIGNMENTS", alignments),
            ("FORMAT_TYPE", format_type),
            ("RID", rid),
            ("CMD", "Get")
        ]
        query = [x for x in parameters if x[1] is not None]
        message = urlencode(query)
        sum_table.loc[index, "Message_get"] = message

    # Poll NCBI until the results are ready.
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
    # 1. Do not contact the server more often than once every 10 seconds.
    # 2. Do not poll for any single RID more often than once a minute.
    # 3. Use the URL parameter email and tool, so that the NCBI
    #    can contact you if there is a problem.
    # 4. Run scripts weekends or between 9 pm and 5 am Eastern time
    #    on weekdays if more than 50 searches will be submitted.
    # --
    # Could start with a 10s delay, but expect most short queries
    # will take longer thus at least 70s with delay. Therefore,
    # start with 20s delay, thereafter once a minute.

    previous = 0
    delay = 20  # seconds
    indices = sum_table[sum_table["Status"] != "Successful"].index
    while len(indices) > 0:
        for n in range(len(indices)):
            index = indices[n]
            # If the waiting time has exceeded 24 hours, try submitting query again.
            time_put = sum_table.loc[index, "Time_put"]
            if time.time() - time_put > 86400:  # 24 h
                sum_table.loc[index, "RID"] = ""
                sum_table.loc[index, "RTOE"] = ""
                sum_table.loc[index, "Message_get"] = ""
                sum_table.loc[index, "Time_get"] = ""
                sum_table.loc[index, "Status"] = ""
                sum_table.loc[index, "Query"] = ""
                sum_table.to_csv(Path(wd) / Path(table), index=False, sep="\t")
                print("The waiting time has exceeded 24 hours, try submitting query %s again..." % sum_table.loc[
                    index, "ID"])
                return False
            try:
                current = time.time()
                wait = previous + delay - current
                if wait > 0:
                    print("Waiting for %.2f seconds..." % wait)
                    time.sleep(wait)
                    previous = current + wait
                else:
                    previous = current

                # delay by at least 60 seconds only if running the request against the public NCBI API
                # ua = UserAgent()
                # header = {"User-Agent": ua.random}

                message = sum_table.loc[index, "Message_get"]
                message = message.encode()
                print("Contacting the server to get results for %s..." % sum_table.loc[index, "RID"])
                # results = get_blast_results(url_base, message, header)
                results = get_blast_results(url_base, message, {"User-Agent": "BiopythonClient"})

                # Can see an "\n\n" page while results are in progress,
                # if so just wait a bit longer...
                if results == "\n\n":
                    print("Results are not ready yet!")  # todo: print this information
                    continue
                if format_type != "XML":
                    # todo: XML results don't have the Status tag when finished
                    if "Status=" not in results:
                        print("Status is not found in results!")
                        continue
                    i = results.index("Status=")
                    j = results.index("\n", i)
                    status = results[i + len("Status="): j].strip()
                    if status.upper() == "READY":
                        print("Writing results...")
                        with open(Path(wd) / Path(sum_table.loc[index, "RID"] + "_Tabular.txt"), "w") as fw:
                            fw.write(results)
                        print("Results saved in %s_Tabular.txt" % sum_table.loc[index, "RID"])
                        sum_table.loc[index, "Status"] = "Successful"
                        sum_table.loc[index, "Time_get"] = time.time()
                        sum_table.to_csv(Path(wd) / Path("blast_summary.txt"), index=False, sep="\t")
                        # todo: sometimes the results are not completely downloaded
                        # for example: IncompleteRead(306872 bytes read)
                    else:
                        print("Results are not ready yet!")
                else:
                    if results.find("</BlastOutput>") == -1:
                        continue
                        # IncompleteRead

                    else:
                        if results.find("<Hit>") == -1:
                            sum_table.loc[index, "RID"] = ""
                            sum_table.loc[index, "RTOE"] = ""
                            sum_table.loc[index, "Message_get"] = ""
                            sum_table.loc[index, "Time_get"] = ""
                            sum_table.loc[index, "Status"] = ""
                            sum_table.loc[index, "Query"] = ""
                            sum_table.to_csv(Path(wd) / Path(table), index=False, sep="\t")
                            print("The results are invalid, try submitting query %s again..." % sum_table.loc[
                                index, "ID"])
                            return False
                        with open(Path(wd) / Path(sum_table.loc[index, "RID"] + "_XML.txt"), "w") as fw:
                            fw.write(results)
                            # todo: IncompleteRead(1347250 bytes read)
                        print("Results saved in %s_XML.txt" % sum_table.loc[index, "RID"])
                        sum_table.loc[index, "Status"] = "Successful"
                        sum_table.loc[index, "Time_get"] = time.time()
                        find_query = re.compile('<BlastOutput_query-ID>(.*?)</BlastOutput_query-ID>')
                        query = re.findall(find_query, results)[0]
                        sum_table.loc[index, "Query"] = query
                        sum_table.to_csv(Path(wd) / Path(table), index=False, sep="\t")

            except func_timeout.exceptions.FunctionTimedOut:
                print(
                    "Time out, try contacting the server to get results for %s again..." % sum_table.loc[index, "RID"])
            except Exception as result:
                print("Error: %s, try contacting the server to get results for %s again..." %
                      (result, sum_table.loc[index, "RID"]))
        indices = sum_table[sum_table["Status"] != "Successful"].index
        if len(indices) == 1:
            # Wasn't a quick return, must wait at least a minute
            delay = 80
        elif len(indices) == 2:
            delay = 40
    print("All BLAST results saved!")
    return True


def blast_put_get_main(wd, queries_path, entrez_query, alignments, expect, gapcosts, word_size, nucl_reward,
                       nucl_penalty):
    while True:
        put_blast(wd=wd,
                  queries_path=queries_path,
                  entrez_query=entrez_query,
                  alignments=alignments,
                  expect=expect,
                  gapcosts=gapcosts,
                  word_size=word_size,
                  nucl_reward=nucl_reward,
                  nucl_penalty=nucl_penalty
                  )
        finished = get_blast(wd=wd, alignments=alignments)
        if finished:
            break
        else:
            pass


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


# def hits_join_extend_main(wd, max_len):
#     # parse xml files
#     file_list = os.listdir(wd)
#     for file in file_list:
#         if file.split("_")[-1] == "XML.txt":
#             try:
#                 print("Parsing %s..." % file, end="")
#                 time0 = datetime.now()
#                 parse_xml_by_re(wd, file)  # save results as HitTable and delete xml files
#                 time1 = datetime.now()
#                 print("Running time: %s Seconds" % (time1 - time0))
#             except Exception as result:
#                 print(result)
#     # join hits
#     file_list = os.listdir(wd)
#     hit_tables = [x.split("_")[0] for x in file_list if x.split("_")[-1] == "HitTable.txt"]
#     extended_hit_tables = [x.split("_")[0] for x in file_list if x.split("_")[-1] == "joined.txt"]
#     not_extended_hit_tables = set(hit_tables) - set(extended_hit_tables)
#     if "blast_summary.txt" in file_list:
#         sum_table = pd.read_table(Path(wd) / Path("blast_summary.txt"), sep='\t', engine='python')
#     else:
#         print("Cound not find blast_summary.txt")
#         return
#     for hit_table in not_extended_hit_tables:
#         file = hit_table + "_HitTable.txt"
#         try:
#             print("Joining %s..." % file, end="")
#             time0 = datetime.now()
#             hit_table = pd.read_table(Path(wd) / Path(file), sep='\t', engine='python')
#             query_ref_len = sum_table[sum_table["RID"] == file.split("_")[0]].iloc[0]["Sequence_length"]
#             hit_table_merged = join_hits(hit_table, max_len, query_ref_len)
#             hit_table_merged.to_csv(Path(wd) / Path(file.split("_")[0] + "_joined.txt"), index=False, sep="\t")
#             time1 = datetime.now()
#             print("Running time: %s Seconds" % (time1 - time0))
#         except Exception as result:
#             print(result)
#     # select hits
#     file_list = os.listdir(wd)
#     if "hits_selected.txt" not in file_list:
#         time0 = datetime.now()
#         print("Reading joined HitTables...", end="")
#         hit_tables = []
#         for file in file_list:
#             if file.split("_")[-1] == "joined.txt":
#                 hit_table = pd.read_table(Path(wd) / Path(file), sep='\t', engine='python')
#                 hit_tables.append(hit_table)
#         hit_tables = pd.concat(hit_tables)  # hit_tables.shape
#         print("Selecting hits...", end="")
#         hits_selected = select_hits(hit_tables)
#         hits_selected.to_csv(Path(wd) / Path("hits_selected.txt"), index=False, sep="\t")
#         time1 = datetime.now()
#         print("Results saved in hits_selected.txt")
#         print("Running time: %s Seconds" % (time1 - time0))
#
#     else:
#         hits_selected = pd.read_table(Path(wd) / Path("hits_selected.txt"), sep='\t', engine='python')
#     # extend hits
#     if "hit_len" not in hits_selected.columns:
#         time0 = datetime.now()
#         print("Extending selected hits...", end="")
#         extended_hit_tables = []
#         for (name, group) in hits_selected.groupby(hits_selected["query_acc.ver"]):
#             query_ref_len = sum_table[sum_table["Query"] == name].iloc[0]["Sequence_length"]
#             missing_left = sum_table[sum_table["Query"] == name].iloc[0]["Missing_left"]
#             missing_right = sum_table[sum_table["Query"] == name].iloc[0]["Missing_right"]
#             extended_hit_tables.append(extend_hits(group, max_len, query_ref_len, missing_left, missing_right))
#         extended_hit_tables = pd.concat(extended_hit_tables)
#         extended_hit_tables.to_csv(Path(wd) / Path("hits_selected.txt"), index=False, sep="\t")
#         time1 = datetime.now()
#         print("Results saved in hits_selected.txt")
#         print("Running time: %s Seconds" % (time1 - time0))
#         return extended_hit_tables
#     else:
#         return hits_selected


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
    # # extend hits
    # if "hit_len" not in hits_selected.columns:
    #     time0 = datetime.now()
    #     print("Extending selected hits...", end="")
    #     extended_hit_tables = []
    #     for (name, group) in hits_selected.groupby(hits_selected["query_acc.ver"]):
    #         query_ref_len = sum_table[sum_table["Query"] == name].iloc[0]["Sequence_length"]
    #         missing_left = sum_table[sum_table["Query"] == name].iloc[0]["Missing_left"]
    #         missing_right = sum_table[sum_table["Query"] == name].iloc[0]["Missing_right"]
    #         extended_hit_tables.append(extend_hits(group, max_len, query_ref_len, missing_left, missing_right))
    #     extended_hit_tables = pd.concat(extended_hit_tables)
    #     extended_hit_tables.to_csv(Path(wd) / Path("hits_selected.txt"), index=False, sep="\t")
    #     time1 = datetime.now()
    #     print("Results saved in hits_selected.txt")
    #     print("Running time: %s Seconds" % (time1 - time0))
    #     return extended_hit_tables
    # else:
    #     return hits_selected
    return hits_selected


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


# def align_new_query(wd, blast_round):
#     # print("Align new queries...")
#     # fr = open(Path(tmp_wd) / Path("new_queries.fasta"), "r")
#     # seq = fr.read()
#     # fr.close()
#     # fw = open(Path(wd) / Path("parameters")/Path("all_queries.fasta"), "a")
#     # fw.write("\n" + seq)
#     # fw.close()
#
#     # first check if the file exist and not zero size
#     # allowing human adjustment of the alignments?
#
#     mafft_cmd = "mafft --multipair --addfragments %s %s > %s" % \
#                 (Path(wd) / Path("parameters") / Path("ref_seq") / Path("queries_%d.fasta" % blast_round),
#                  Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta"),
#                  Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_%d.fasta" % blast_round))
#     os.system(mafft_cmd)
#     msa1 = AlignIO.read(Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_1.fasta"), "fasta")
#     msa2 = AlignIO.read(Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_%d.fasta" % blast_round),
#                         "fasta")
#     i = len(msa1)
#     fw = open(Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_%d.fasta" % blast_round), "w")
#     fw.close()
#     while i < len(msa2):
#         with open(Path(wd) / Path("parameters") / Path("ref_msa") / Path("msa_queries_%d.fasta" % blast_round),
#                   "a") as fw:
#             fw.write(">" + msa2[i].description + "\n")
#             fw.write(str(msa2[i].seq) + "\n")
#             i += 1


@func_set_timeout(600)
def my_efetch(accession, strand, seq_start, seq_stop):
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession,
                           strand=strand, seq_start=seq_start, seq_stop=seq_stop)
    print("Extended start: %d, " % seq_start, end="")
    print("Extended end: %d, " % seq_stop, end="")
    print("Strand: %d, " % strand, end="")
    return handle


def check_annotation(feature_list, key_annotations, exclude_sources):
    for exclude_source in exclude_sources:
        for feature in feature_list:
            if feature.find(exclude_source) >= 0:
                return False
    for key_annotation in key_annotations:
        for feature in feature_list:
            if feature.find(key_annotation) >= 0:
                return True
    return False


def write_seq_info(record, start, end, strand, wd, file):
    # todo: bould
    accession = record.id
    length = str(record.seq).lower().count("a")
    length += str(record.seq).lower().count("t")
    length += str(record.seq).lower().count("c")
    length += str(record.seq).lower().count("g")
    date = record.annotations["date"]
    description = record.description
    source = record.annotations["source"]
    organism = record.annotations["organism"]
    taxonomy = record.annotations["taxonomy"]
    taxonomy_str = ""
    for x in taxonomy:
        taxonomy_str += x + "|"
    taxonomy_str += organism
    if "references" in record.annotations.keys():
        reference = record.annotations["references"][0]
        title = reference.title
        authors = reference.authors
        journal = reference.journal
    else:
        title = "unknown"
        authors = "unknown"
        journal = "unknown"

    feature = record.features[0]

    if "organelle" in feature.qualifiers.keys():
        organelle = str(feature.qualifiers["organelle"][0])
    else:
        organelle = "unknown"

    if "mol_type" in feature.qualifiers.keys():
        mol_type = str(feature.qualifiers["mol_type"][0])
    else:
        mol_type = "unknown"

    # todo: the gb file may contain two db_xref, taxon and bold
    if "db_xref" in feature.qualifiers.keys():
        db_xref = "|".join(feature.qualifiers["db_xref"])
    else:
        db_xref = "unknown"

    if "specimen_voucher" in feature.qualifiers.keys():
        specimen_voucher = str(feature.qualifiers["specimen_voucher"][0])
    else:
        specimen_voucher = "unknown"

    if "country" in feature.qualifiers.keys():
        country = str(feature.qualifiers["country"][0])
    else:
        country = "unknown"

    if "lat_lon" in feature.qualifiers.keys():
        lat_lon = str(feature.qualifiers["lat_lon"][0])
    else:
        lat_lon = "unknown"

    if "collection_date" in feature.qualifiers.keys():
        collection_date = str(feature.qualifiers["collection_date"][0])
    else:
        collection_date = "unknown"

    if "collected_by" in feature.qualifiers.keys():
        collected_by = str(feature.qualifiers["collected_by"][0])
    else:
        collected_by = "unknown"

    if "identified_by" in feature.qualifiers.keys():
        identified_by = str(feature.qualifiers["identified_by"][0])
    else:
        identified_by = "unknown"

    with open(Path(wd) / Path(file), "a") as fw:
        fw.write(accession + "\t" + str(start) + "\t" + str(end) + "\t" + str(strand) + "\t" + str(
            length) + "\t" + date + "\t")
        fw.write(description + "\t" + source + "\t" + organism + "\t" + taxonomy_str + "\t")
        fw.write(title + "\t" + authors + "\t" + journal + "\t" + organelle + "\t" + mol_type + "\t" + db_xref + "\t")
        fw.write(specimen_voucher + "\t" + country + "\t" + lat_lon + "\t")
        fw.write(collection_date + "\t" + collected_by + "\t" + identified_by + "\n")


def write_fas_file(record, start, end, strand, wd, file):
    description = record.description
    organism = record.annotations["organism"].replace(" ", "_")
    if start <= 0:
        start = 1
    end = len(record.seq) + start - 1
    if strand == 1:
        # if start <= 0:
        #     start = 1
        # end = len(record.seq) + start - 1
        fas_description = ">%s:%d-%d|%s|%s" % (record.id, start, end, organism, description)
        fas_seq = str(record.seq)
    else:
        # if end <= 0:
        #     end = 1
        # start = len(record.seq) + end - 1
        # fas_description = ">%s:%d-%d_reverse_complement|%s|%s" % (record.id, end, start, organism, description)
        fas_description = ">%s:%d-%d_reverse_complement|%s|%s" % (record.id, start, end, organism, description)
        fas_seq = str(record.seq)  # 659
    with open(Path(wd) / Path(file), "a") as fw:
        fw.write(fas_description + "\n")
        fw.write(fas_seq + "\n")
    write_seq_info(record, start, end, strand, wd, file=Path(file).stem + "_seq_info.txt")


def seq_check_download(wd, acc_file, out_file, key_annotations, exclude_sources, entrez_email):
    """download fasta files from Genbank according to given accessions"""
    # todo: use user provided email
    Entrez.email = entrez_email
    acc_list = list(acc_file.index)
    while len(acc_list) != 0:
        accession = acc_list[0]
        try:
            t0 = datetime.now()
            seq_start = int(acc_file.loc[accession, "start"])
            seq_stop = int(acc_file.loc[accession, "end"])
            strand = int(acc_file.loc[accession, "strand"])
            if strand == 2:
                seq_start, seq_stop = (seq_stop, seq_start)
            # print("start： %d" % start)
            # print("end： %d" % end)
            handle = my_efetch(accession=accession, strand=strand, seq_start=seq_start, seq_stop=seq_stop)
            record = SeqIO.read(handle, "gb")
            print("Actual sequence length: %d" % len(record.seq))
            feature_list = []
            for feature in record.features:
                feature_list.extend(feature.qualifiers.values())
                # print(feature.qualifiers.values())
            feature_list = [x[0].replace(" ", "").lower() for x in feature_list]
            key_annotations = [x.replace(" ", "").lower() for x in key_annotations]
            exclude_sources = [x.replace(" ", "").lower() for x in exclude_sources]
            if check_annotation(feature_list, key_annotations, exclude_sources):
                write_fas_file(record, seq_start, seq_stop, strand, wd, file=out_file)
            else:
                if not os.path.exists(Path(wd) / Path("erroneous_" + Path(out_file).stem + "_seq_info.txt")):
                    with open(Path(wd) / Path("erroneous_" + Path(out_file).stem + "_seq_info.txt"), "w") as fw:
                        fw.write(
                            "accession\tstart\tend\tstrand\tlength\tdate\tdescription\tsource\torganism\ttaxonomy\t")
                        fw.write("title\tauthors\tjournal\torganelle\tmol_type\tdb_xref\t")
                        fw.write("specimen_voucher\tcountry\tlat_lon\t")
                        fw.write("collection_date\tcollected_by\tidentified_by\n")
                write_fas_file(record, seq_start, seq_stop, strand, wd, file="erroneous_" + out_file)
            t1 = datetime.now()
            print("%s downloaded in %s seconds" % (accession, t1 - t0))
            acc_list.pop(0)
        except ValueError:  # features location no correct
            acc_list.pop(0)
            with open(Path(wd) / Path("value_error_list.txt"), "a") as fw:
                fw.write(accession + "\n")
            print("ValueError: CompoundLocation should have at least 2 parts, skip %s" % accession)
        except urllib.error.HTTPError:  # HTTP Error 400, wrongly parsed accession
            acc_list.pop(0)
            with open(Path(wd) / Path("bad_request_list.txt"), "a") as fw:
                fw.write(accession + "\n")
            print("Bad request: %s. Move on to the next sequence." % accession)
        except func_timeout.exceptions.FunctionTimedOut:
            print("Time out, try downloading %s again..." % accession)
        except Exception as result:
            print("Error: %s, try downloading %s again..." % (result, accession))


def filter_duplicate_key(wd, file):
    if os.path.getsize(Path(wd) / Path(file)) > 0:
        os.system("copy %s %s" % (str(Path(wd) / Path(file)), str(Path(wd) / Path("tmp_" + file))))
    key_list = []
    with open(Path(wd) / Path(file), "w") as fw:
        for record in SeqIO.parse(Path(wd) / Path("tmp_" + file), "fasta"):
            if record.description not in key_list:
                fw.write(">" + record.description + "\n")
                fw.write(str(record.seq) + "\n")
                key_list.append(record.description)
            else:
                print("Filtered duplicate sequence: %s" % record.description)
    os.system("del %s" % str(Path(wd) / Path("tmp_" + file)))


def seq_check_download_main(wd, acc_file, out_file, key_annotations, exclude_sources, entrez_email, extend=False):
    # print_line()
    print("Downloading sequences...")
    df = pd.read_table(Path(wd) / Path(acc_file), sep="\t", engine="python")
    # drop duplicate
    df.index = df["subject_acc.ver"]
    file_list = os.listdir(wd)

    # def get_accession(record):
    #     accession = record.description.split(" ")[0].split("|")[0].split(":")[0]
    #     return accession

    index_list = list(df.index)
    if out_file in file_list:
        filter_duplicate_key(wd, out_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(Path(wd) / Path(out_file), "fasta"), key_function=get_query_accession)
        print("Correct sequences already downloaded: %d" % len(seq_dict.keys()))
        index_list = list(set(index_list) - set(seq_dict.keys()))
    if "erroneous_" + out_file in file_list:
        filter_duplicate_key(wd, "erroneous_" + out_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(Path(wd) / Path("erroneous_" + out_file), "fasta"),
                                 key_function=get_query_accession)
        print("Erroneous sequences already downloaded: %d" % len(seq_dict.keys()))
        index_list = list(set(index_list) - set(seq_dict.keys()))
    if Path(out_file).stem + "_seq_info.txt" not in file_list:
        with open(Path(wd) / Path(Path(out_file).stem + "_seq_info.txt"), "w") as fw:
            fw.write("accession\tstart\tend\tstrand\tlength\tdate\tdescription\tsource\torganism\ttaxonomy\t")
            fw.write("title\tauthors\tjournal\torganelle\tmol_type\tdb_xref\t")
            fw.write("specimen_voucher\tcountry\tlat_lon\t")
            fw.write("collection_date\tcollected_by\tidentified_by\n")
    else:
        seq_info = pd.read_table(Path(wd) / Path(Path(out_file).stem + "_seq_info.txt"), sep="\t", engine="python")
        seq_info.drop_duplicates(subset=['accession'], keep='first', inplace=True)
        seq_info.to_csv(Path(wd) / Path(Path(out_file).stem + "_seq_info.txt"), sep="\t", index=False)

    index_list.sort()
    print("Sequences to download: %d" % len(index_list))
    if not extend:
        df1 = df.loc[index_list][["s_start", "s_end", "s_strand"]].copy()
    else:
        df1 = df.loc[index_list][["s_extstart", "s_extend", "s_strand"]].copy()
    df1.columns = ['start', 'end', 'strand']
    df1["strand"] = df1["strand"].replace([True, False], [1, 2])
    # df1[["s_extstart", "s_extend", "s_strand"]] = df1[["s_extstart", "s_extend", "s_strand"]].apply(pd.to_numeric)

    fw = open(Path(wd) / Path("value_error_list.txt"), "w")
    fw.close()
    fw = open(Path(wd) / Path("bad_request_list.txt"), "w")
    fw.close()

    seq_check_download(wd=wd, acc_file=df1, out_file=out_file,
                       key_annotations=key_annotations, exclude_sources=exclude_sources, entrez_email=entrez_email)

    fw = open(Path(wd) / Path("value_error_list.txt"), "r")
    value_error_list = fw.read().splitlines()
    if len(value_error_list) > 0:
        print("%d value errors, the accession numbers were save in value_error_list.txt" % len(value_error_list))
    fw.close()
    fw = open(Path(wd) / Path("bad_request_list.txt"), "r")
    bad_request_list = fw.read().splitlines()
    if len(bad_request_list) > 0:
        print("%d bad requests, the accession numbers were save in bad_request_list.txt" % len(bad_request_list))
    fw.close()

    print("All sequences successfully downloaded!")
    print("Correct sequences were saved in %s." % out_file)
    print("Erroneous sequences were saved in %s" % "erroneous_" + out_file)


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
        if blast_count > count * 1.1:
            print("Total sequences number is greater than Entrez search results number. Stop BLAST iteration.")
            break
        if blast_count_new1 < 3 and last_new < 3:
            print("Cannot find more than 2 correct new sequences. Stop BLAST iteration.")
            break
        last_new = blast_count_new1

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
    calculate_missing_length(wd, ref_msa_file)

    blast_results_extend_main(wd, max_length)

    seq_check_download_main(wd=Path(wd) / Path("results"), acc_file=r"blast_results.txt",
                            out_file=r"blast_results_checked.fasta",
                            key_annotations=key_annotations, exclude_sources=exclude_sources, entrez_email=entrez_email,
                            extend=True)
    print("Stop thread.")


########################################################################################################################
#  under development

def my_concatenation(in_path, out_path):
    t0 = datetime.now()
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    file_list = os.listdir(in_path)
    file_list = [file for file in file_list if Path(file).suffix in [".fasta", ".fas", ".fa"]]
    # todo: the description of trimal output seqs contain blank space?
    for file in file_list:
        fr = open(Path(in_path) / Path(file), "r")
        records = fr.read().replace(" ", "")
        fr.close()
        fw = open(Path(in_path) / Path(file), "w")
        fw.write(records)
        fw.close()
    taxon_completion(in_path=in_path, out_path=out_path)
    concat(in_path=os.path.join(out_path, "completion.result"), out_path=out_path)
    fas2phy(in_path=os.path.join(out_path, "concat.result"), out_path=out_path)
    t1 = datetime.now()
    print("Running time: %s seconds" % (t1 - t0))


def check_inpath_validity(path):
    """check if an input path is valid (an existing path except shortcuts and .result files)
    ----------
    Parameters
    - path - the path whose validity is to be checked
    -------
    Returns
    True if path is valid
    False if path is NOT valid"""
    if not isinstance(path, str):  # invalid input: value type
        return False
    if path.endswith('.lnk'):
        return False  # do nothing with shortcut files
    # if path.endswith('.result'):  # any folder or file endswith .result will not be loaded
    #     return False
    if not os.path.exists(path):  # invalid input: no such path
        return False
    return True


def check_outpath_validity(path):
    """check if an output path is valid (an existing folder except shortcuts)
    ----------
    Parameters
    - path - the path whose validity is to be checked
    -------
    Returns
    True if path is valid
    False if path is NOT valid"""
    if not isinstance(path, str):
        return False
    if path.endswith('.lnk'):
        return False  # do nothing with shortcut files
    if not os.path.exists(path):  # invalid input: no such path
        return False
    if not os.path.isdir(path):
        return False  # must be an existing folder
    return True


def create_folder(path):
    """create a new folder if one with the same name hasn't been created.
    ----------
    Parameters
    - path - the path of the folder to be created
    -------
    Returns
    - path - if folder is created
    0 - if not
    -1 if illegal characters are in the path"""
    try:
        if not os.path.exists(path):
            os.mkdir(path)
            return path
        else:
            return 0
    except:
        return -1


def get_file_handles(in_path):
    """in_path could be a file/folder, a list of files/folders or a list of both
    Therefore, there is need to consider each condition and get all file handles
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    -------
    Returns
    - file_handles - all file handles in a list"""
    # STEP 1: if in_path is not a list, then make it a list
    if not isinstance(in_path, list):
        in_path = [in_path]

    # STEP 2: for all elements in the given list, decide whether it's a file or a folder:
    abs_handles = {}  # in case duplicates but in different forms like abs/relative
    file_handles = []  # result, returned value
    for file_path in in_path:
        if not check_inpath_validity(file_path):
            continue
        # substep 1 : if an element is a file, add the file to result
        if os.path.isfile(file_path):
            file_handles.append(file_path)
        # substep 2 : if an element is a folder, get files directly in the folder
        else:  # is a folder
            file_handles += [os.path.join(file_path, file) for file in os.listdir(file_path)
                             if (os.path.isfile(os.path.join(file_path, file))
                                 and not file.endswith('.lnk'))]  # exclude shortcut

    # STEP 3: replace reverse slash, remove duplicates and sort
    file_handles = list(map(lambda x: x.replace('\\', '/'), file_handles))
    for handle in file_handles:
        abs_handles.setdefault(os.path.abspath(handle), [])
        abs_handles[os.path.abspath(handle)].append(handle)
    file_handles = []
    for handle_list in abs_handles.values():
        min_len = min(list(map(lambda x: len(x), handle_list)))
        for handle in handle_list:
            if len(handle) == min_len:
                file_handles.append(handle)
                break
    file_handles.sort()
    file_handles = [x for x in file_handles if x.split(".")[-1] in ["fasta", "fas", "fa", "phylip", "phy"]]

    return file_handles


def filter_by_length(in_path, out_path='./'):
    """filter fasta file(s) to keep only the longest sequence for each taxon
    When multiple files are input, do filtering for each file.
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    -------
    Returns
    -written_lengths - a list of number of sequences written (for each file)"""
    # STEP 1: get ready: check paths and create folders, make list of inputs
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return []  # write nothing if output folder is invalid
    written_lengths = []  # to store the number of records written, returned value

    # STEP 2: for each file in the list, for each taxon,
    # keep information of records with the longest sequence using dictionary
    for in_file in file_handles:
        t0 = datetime.now()
        longest_records = {}
        record_iter = SeqIO.parse(in_file, 'fasta')
        for record in record_iter:
            if len(record.description.split('|')) == 1:
                taxon_name = record.description
            else:
                taxon_name = record.description.split('|')[1]
            longest_records.setdefault(taxon_name, SeqRecord.SeqRecord(''))
            if len(record.seq) > len(longest_records[taxon_name].seq):
                longest_records[taxon_name] = record

        # STEP 3: write a csv file to store the detailed information of kept records
        # temp_output_folder = os.path.join(out_path, 'filtered.result')
        temp_output_folder = out_path
        create_folder(temp_output_folder)
        cols = ['taxon name', 'description']
        logs = []  # kept records, used to create a dataframe and write to csv
        for key, value in longest_records.items():  # key:taxon_name, value:record itself
            logs.append([key, value.description])
        df = pd.DataFrame(logs, columns=cols)
        df.to_csv(os.path.join(temp_output_folder,
                               os.path.splitext(os.path.basename(in_file))[0]
                               + '_kept_records.csv'),
                  index=False)

        # STEP 4: write a fasta file to store only name and sequence of kept records
        out_file = os.path.join(temp_output_folder, os.path.basename(in_file))
        for taxon_name in longest_records.keys():  # avoid unwanted output
            longest_records[taxon_name].name = ''
            longest_records[taxon_name].id = taxon_name.strip()
            longest_records[taxon_name].description = ''
        written_lengths.append(SeqIO.write(longest_records.values(), out_file, 'fasta'))
        t1 = datetime.now()
        print("Filtered results: %s " % out_file)
        print("Running time: %s seconds" % (t1 - t0))

    return written_lengths  # if no file is written, return value is an empty list


def fas2phy(in_path, out_path='./'):
    """ convert fasta file to phylip file.
    When multiple files are input, do transformation for each file.
    (If a fasta is an unaligned collection of sequences, nothing will be done.
    If a fasta is an empty fasta, write an empty file accordingly.)
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    -------
    Returns
    - written_results - a list of pair [species_number, sequence_length]
        where species_number is the number of species in a fasta
        and sequence length is the length of sequences in this fasta"""
    # STEP 0: check validity of in_path and out_path
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return []  # write nothing if output folder is invalid
    out_path = os.path.join(out_path, 'phylip.result')
    create_folder(out_path)
    written_results = []  # to store the pair [species_number, sequence_length]

    # STEP 1: read and get length of the longest name, species number and sequence length.
    # and also checking whether all sequences are of the same length
    for in_file in file_handles:
        record_iter = SeqIO.parse(in_file, 'fasta')
        aligned = True
        len_longest_name = 0
        species_number = 0
        sequence_length = 0

        for record in record_iter:
            taxon_name = record.description
            taxon_length = len(record.seq)
            if species_number == 0:  # initiate sequence length
                sequence_length = taxon_length
            if taxon_length != sequence_length:  # phylip requires sequences in same length
                aligned = False
                break
                # raise Exception('Invalid sequence length: ', str(in_file))
            if len(taxon_name) > len_longest_name:
                len_longest_name = len(taxon_name)
            species_number += 1

        # STEP 2: re-read input file and write to phylip format.
        if not aligned:  # fasta is not aligned, cannot transform to phylip, so skip
            continue
        illegal_characters = ["\t", "\n", " ", ":", ",", ")", "(", ";", "]", "[", "'"]
        filename = os.path.splitext(os.path.basename(in_file))[0] + '.phy'
        fptr = open(os.path.join(out_path, filename), 'w')
        fptr.write(f' {species_number}  {sequence_length}' + '\n')  # header
        record_iter = SeqIO.parse(in_file, 'fasta')
        name_room = len_longest_name + 1  # room for each taxon's name and following spaces
        for record in record_iter:
            taxon_name = record.description
            taxon_seq = str(record.seq)

            # replace illegal characters in phylip by underlines if necessary.
            if any(substring in taxon_name for substring in illegal_characters):
                illegals = [substring for substring in taxon_name
                            if substring in illegal_characters]
                for char in illegals:
                    taxon_name = taxon_name.replace(char, '_')
                # to avoid gathering '_' in taxon name
                taxon_name = '_'.join([part for part in taxon_name.split('_') if part != '']) \
                    .strip()

            space = ' ' * (name_room - len(taxon_name))
            fptr.write(taxon_name + space + taxon_seq + '\n')
        fptr.close()
        written_results.append([species_number, sequence_length])
    print("Phylip format results saved in: %s " % out_path)
    return written_results  # pair [species_number, sequence_length]


def taxon_completion(in_path, out_path='./'):
    """Complete fasta files where there are missing records
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    -------
    Returns
    - missing_records - a dictionary of {'input file name': [missing taxa]}"""

    # STEP 0: check validity of in_path and out_path, get legal input handles
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return {}  # write nothing if output folder is invalid
    out_path = os.path.join(out_path, 'completion.result')
    create_folder(out_path)
    present_records = {}  # to store the missing taxon for each gene
    missing_records = {}  # the '-'-completed taxa
    sequence_lengths = {}  # to store the lengths of sequences, -1 for varied length

    # STEP 1: create a list of all taxa(descriptions) in input files
    total_taxa = []
    for in_file in file_handles:
        file_basename = os.path.splitext(os.path.basename(in_file))[0]
        record_iter = SeqIO.parse(in_file, 'fasta')
        recorded_length = list(set([len(record.seq) for record in record_iter]))
        if len(recorded_length) == 1:  # aligned or same-length sequences:
            sequence_lengths.setdefault(file_basename, recorded_length[0])
        elif len(recorded_length) == 0:  # 0 means no records
            sequence_lengths.setdefault(file_basename, 0)
        else:  # not 0/1 means differ in length,
            sequence_lengths.setdefault(file_basename, -1)
        record_iter = SeqIO.parse(in_file, 'fasta')
        recorded_taxa = [record.description for record in record_iter]
        present_records.setdefault(file_basename, recorded_taxa)  # recorded taxa
        total_taxa += recorded_taxa
    total_taxa = list(set(total_taxa))  # remove duplicates
    total_taxa.sort()

    # STEP 2: copy source files to out_path
    for in_file in file_handles:
        try:
            shutil.copyfile(in_file, os.path.join(out_path, os.path.basename(in_file)))
        except IOError as e:
            print("Unable to copy file. %s" % e)
        except:
            print("Unexpected error:", sys.exc_info())

    # STEP 3: append '-'s to the end of these copied files if a marker is not present
    for in_file in file_handles:
        file_basename = os.path.splitext(os.path.basename(in_file))[0]
        recorded_taxa = present_records[file_basename]
        missing_records.setdefault(file_basename,
                                   [taxon for taxon in total_taxa if taxon not in recorded_taxa])
        unrecorded_taxa = missing_records[file_basename]
        recorded_length = sequence_lengths[file_basename]
        # condition 1: if file is aligned, then add equivalent amount of - to keep them aligned
        if recorded_length > 0:
            f = open(os.path.join(out_path, os.path.basename(in_file)), 'a')
            for taxon in unrecorded_taxa:
                f.write(f"\n>{taxon}\n{'-' * recorded_length}")
            f.close()
        # condition 2: this indicates that this file contains nothing, thus ignore it
        elif recorded_length == 0:
            try:
                os.remove(os.path.join(out_path, os.path.basename(in_file)))
            except:
                pass
            finally:
                continue
        # condition 3: if file is not aligned, then for each gap, only three - will be added
        else:
            f = open(os.path.join(out_path, os.path.basename(in_file)), 'a')
            for taxon in unrecorded_taxa:
                f.write(f"\n>{taxon}\n{'-' * 3}")
            f.close()

    # STEP 4: write added records to log file
    f = open(os.path.join(out_path, 'completion.log'), 'w')
    f.write("Appended taxon for each fasta file are listed below:\n")
    for key, value in missing_records.items():
        value.sort()  # order the result
        f.write(f"{key}:  {', '.join([name for name in value])}\n")
    f.close()
    print("Taxon completion results saved in: %s " % out_path)

    return missing_records  # {marker: [missing taxa]}


def concat(in_path, out_path='./', filename='concat.fasta'):
    """
    Concatenate different fasta files according to taxon name,
        to link sequences from different files together.
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
        DESCRIPTION.
    - out_path - destination folder of output, where log files and result will be written.
    - filename - the name of concatenated file and the file name in provided cfg file
    -------
    Returns
    - marker_record - the start and end position of each marker"""
    # STEP 0: check validity of in_path and out_path, get legal input handles
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return {}  # write nothing if output folder is invalid
    out_path = os.path.join(out_path, 'concat.result')
    create_folder(out_path)

    # STEP 1: check whether there is a 'marker gap' for all taxa and whether same length
    # substep 1: collect all taxa from input files
    all_taxa = []
    for in_file in file_handles:
        record_iter = SeqIO.parse(in_file, 'fasta')
        all_taxa += [record.description for record in record_iter]
        record_iter = SeqIO.parse(in_file, 'fasta')
        gene_length = [len(record.seq) for record in record_iter]
        if len(list(set(gene_length))) > 1:
            return {}  # markers with same length required

    # substep2: test missing taxon
    all_taxa = list(set(all_taxa))
    for in_file in file_handles:
        record_iter = SeqIO.parse(in_file, 'fasta')
        num_of_taxa = len([record for record in record_iter])
        if len(all_taxa) != num_of_taxa:
            return {}  # cannot concat when there is missing taxon

    # STEP 2: concat all sequences using 'description' as symbol, write concat result
    concat_result = {}  # {record.description: record}
    marker_record = {}  # {marker: position} where position: [start, end]
    for in_file in file_handles:
        marker_record.setdefault('total_length', 0)
        marker_name = os.path.splitext(os.path.basename(in_file))[0]
        marker_record.setdefault(marker_name, [])
        record_iter = SeqIO.parse(in_file, 'fasta')
        for record in record_iter:
            concat_result.setdefault(record.description,
                                     SeqRecord.SeqRecord('', id=record.id,
                                                         description=''))
            concat_result[record.description].seq += record.seq
            if marker_record[marker_name] == []:
                length = len(record.seq)
                marker_record[marker_name] = [marker_record['total_length'] + 1,
                                              marker_record['total_length'] + length]
                marker_record['total_length'] += length
    del marker_record['total_length']

    SeqIO.write(concat_result.values(), os.path.join(out_path, filename), 'fasta')

    # STEP 3: prepare log and cfg file to record gene positions
    cfg_filename = os.path.splitext(filename)[0]
    header = f"""## ALIGNMENT FILE ##
alignment = {cfg_filename}.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = GTR, GTR+G, GTR+I+G;

# MODEL SELECTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
"""

    tail = """
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;"""

    body = ""
    for key, value in marker_record.items():
        if value == []:
            continue
        for i in range(3):
            body += fr"{key}_pos{i + 1} = {i + value[0]}-{value[1]}\3;" + "\n"
    cfg = header + body + tail
    f = open(os.path.join(out_path, f'{cfg_filename}.cfg'), 'w')
    f.write(cfg)
    f.close()
    print("Concatenation results saved in: %s " % out_path)
    return marker_record


def mafft_add(in_path, in_file, out_path, cmd_str):
    print("Aligning %s..." % in_file)
    len_list = []
    for record in SeqIO.parse(Path(in_path) / Path(in_file), "fasta"):
        len_list.append((record.description, len(record.seq)))

    def take_2(elem):
        return elem[1]

    len_list.sort(key=take_2, reverse=True)

    if len(len_list) > 100:
        # if len(len_list) > 100:
        a = len_list[99][1]
        b = a * 0.5

        file1 = Path(out_path) / Path("long_" + in_file)
        file2 = Path(out_path) / Path("add1_" + in_file)
        file3 = Path(out_path) / Path("add2_" + in_file)

        fw1 = open(file1, "w")
        fw2 = open(file2, "w")
        fw3 = open(file3, "w")

        # for record in SeqIO.parse(Path(in_path)/Path(in_file), "fasta"):
        #     if len(record.seq) >= a:
        #         SeqIO.write(record, fw1, "fasta")
        #     elif len(record.seq) > b:
        #         SeqIO.write(record, fw2, "fasta")
        #     else:
        #         SeqIO.write(record, fw3, "fasta")

        # addfragments first
        long_count = 0
        for record in SeqIO.parse(Path(in_path) / Path(in_file), "fasta"):
            if len(record.seq) >= a and long_count < 100:
                SeqIO.write(record, fw1, "fasta")
                long_count += 1
            elif len(record.seq) > b:
                SeqIO.write(record, fw3, "fasta")
            else:
                SeqIO.write(record, fw2, "fasta")

        fw1.close()
        fw2.close()
        fw3.close()

        msa1 = Path(out_path) / Path("msa1_" + in_file)
        msa2 = Path(out_path) / Path("msa2_" + in_file)
        msa3 = Path(out_path) / Path("msa_" + in_file)

        # os.system("mafft --localpair --maxiterate 1000 %s > %s" % (file1, msa1))
        # os.system("mafft --auto --add %s %s > %s" % (file2, msa1, msa2))  # FFT - NS - 2(Fast but rough)
        # os.system("mafft --auto --addfragments %s %s > %s" % (file3, msa2, msa3))  # Multi-INS-fragment
        os.system(" %s %s > %s" % (cmd_str[0], file1, msa1))

        # os.system("%s --auto --add %s %s > %s" % (cmd_str[1], file2,  msa1, msa2))  # FFT - NS - 2(Fast but rough)
        # os.system("%s --auto --addfragments %s %s > %s" % (cmd_str[2], file3, msa2, msa3))  # Multi-INS-fragment

        # add fragments first
        os.system(
            "%s --auto --addfragments %s %s > %s" % (cmd_str[1], file2, msa1, msa2))  # FFT - NS - 2(Fast but rough)
        if os.path.getsize(msa2) == 0:
            os.system("%s --auto --add %s %s > %s" % (cmd_str[1], file2, msa1, msa2))
        os.system("%s --auto --add %s %s > %s" % (cmd_str[2], file3, msa2, msa3))  # Multi-INS-fragment

        os.remove(file1)
        os.remove(file2)
        os.remove(file3)
        os.remove(msa1)
        os.remove(msa2)
        print("Aligned results: %s" % msa3)
    else:
        os.system(" %s %s > %s" % (cmd_str[0], Path(in_path) / Path(in_file), Path(out_path) / Path("msa_" + in_file)))


def mafft(in_path, out_path='', add_choice='', add_path='', algorithm='auto',
          thread=-1, reorder=True, additional_params='',
          pure_command_mode=False, pure_command=''):
    """
    call mafft to do multiple sequence alignment
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    - add_choice - corresponding to mafft param [add, addfragments, addfull]
    - add_path - corresponding to mafft param new_sequences and in_path becomes aligned sequences
    - algorithm - corresponding to mafft param Algorithm,
        see manual https://mafft.cbrc.jp/alignment/software/manual/manual.html
    - thread - the number of threads, -1 if unsure
    - reorder - reorder the aligned sequences, False to keep input oreder
    - additional_params - additional parameters in the form of command
    - pure_command_mode - if True, only run commands in the textbox, one command per line
    - pure_command - run only if pure_command_mode is True, replace the GUI operations
    -------
    Returns
    - file_handles - the valid input files if not in pure command mode
    - commands - the commands in pure command mode
    [] if path invalid"""
    # STEP 0: if pure command, then only execute input command
    if pure_command_mode == True:
        pure_command = pure_command.split('\n')
        for command in pure_command:
            os.system(command)
        return pure_command

    # STEP 1: check path validity and get file handles
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return []
    if add_choice and not check_inpath_validity(add_path):  # if add is True, then the following path should be valid
        return []

    # STEP 2: get parameters and call mafft
    # for in_file in file_handles:
    #     basename = os.path.basename(in_file)
    #     out_file = os.path.join(out_path, basename)
    #     if add_choice:
    #         command = f"mafft --{algorithm} --{add_choice} {add_path} --thread {thread} {'--reorder' * reorder} {additional_params} {in_file} > {out_file}"
    #     else:
    #         command1 = f"mafft --{algorithm} --thread {thread} {'--reorder' * reorder} {additional_params} {in_file} > {out_file}"
    #     os.system(command)

    # multiple files in a directory
    if os.path.isdir(in_path):
        file_list = os.listdir(in_path)
    # one file
    else:
        file_list = [os.path.basename(in_path)]
        in_path = os.path.dirname(in_path)
    file_list = [x for x in file_list if Path(x).suffix in [".fasta", ".fas", ".fa"]]
    if len(file_list) == 0:
        print("Could not find any fasta file in the input.")
        return

    if algorithm == "auto":
        for file in file_list:
            t0 = datetime.now()
            in_file = os.path.join(in_path, file)
            out_file = os.path.join(out_path, "msa_" + file)
            command = f"mafft --auto --thread {thread} {'--reorder' * reorder} {additional_params} {in_file} > {out_file}"
            # print(command)
            print("Aligning %s..." % file)
            os.system(command)
            t1 = datetime.now()
            print("Running time: %s seconds" % (t1 - t0))
    else:
        command1 = f"mafft --localpair --maxiterate 1000 --thread {thread} {'--reorder' * reorder} {additional_params}"
        command2 = f"mafft --thread {thread} {'--reorder' * reorder} {additional_params}"
        command3 = f"mafft --thread {thread} {'--reorder' * reorder} {additional_params}"
        cmd_str = [command1, command2, command3]

        for file in file_list:
            t0 = datetime.now()
            mafft_add(in_path, file, out_path, cmd_str)
            t1 = datetime.now()
            print("Running time: %s seconds" % (t1 - t0))

    return file_handles


def trimal(in_path, out_path='',
           htmlout=True, bp_length=False,
           implement_methods='automated1', gt='', st='', ct='', cons='',
           additional_params='', pure_command_mode=False, pure_command=''):
    """
    call trimal to trim a set of aligned sequences
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    - htmlout - if True, produce report(s) in html format.
    - bp_length - if True, the description of each record in output file will be followed by kept length after trimal
    - implement_methods - corresponding to trimal params [gappyout, strict, strictplus, automated1]
    - gt - corresponding to trimal params gt: gapthreshold, 1 - (fraction of sequences with a gap allowed).
    - st - corresponding to trimal params st: simthreshold, Minimum average similarity allowed.
    - ct - corresponding to trimal params ct: conthreshold, Minimum consistency value allowed.
    - cons - corresponding to trimal params cons: Minimum percentage of the positions in the original alignment to conserve.
    - additional_params - additional parameters in the form of command
    - pure_command_mode - if True, only run commands in the textbox, one command per line
    - pure_command - run only if pure_command_mode is True, replace the GUI operations
    -------
    Returns
    - file_handles - the valid input files if not in pure command mode
    - commands - the commands in pure command mode
    [] if path invalid"""
    # STEP 0: if pure command, then only execute input command
    if pure_command_mode:
        pure_command = pure_command.split('\n')
        for command in pure_command:
            os.system(command)
        return pure_command

    # STEP 1: check path validity and get file handles
    file_handles = get_file_handles(in_path)
    if not check_outpath_validity(out_path):
        return []

    # STEP 2: get parameters and call trimal

    for in_file in file_handles:
        t0 = datetime.now()
        basename = os.path.basename(in_file)
        out_file = os.path.join(out_path, basename)
        command = ["trimal -in %s -out %s" % (in_file, out_file)]
        # command = [f"{trimal_path} -in {in_file} -out {out_file}"]

        if htmlout:
            html_folder = os.path.join(out_path, 'htmlout')
            create_folder(html_folder)
            html_out_file = os.path.join(html_folder, os.path.splitext(basename)[0] + '.html')
            command.append(f"-htmlout {html_out_file}")

        if implement_methods:
            command.append(f"-{implement_methods}")
        if gt:
            command.append(f"-gt {gt}")
        if st:
            command.append(f"-st {st}")
        if ct:
            command.append(f"-ct {ct}")
        if cons:
            command.append(f"-cons {cons}")

        command.append(additional_params)
        command = " ".join(command)
        # print(command)
        print("Trimming %s..." % basename)
        os.system(command)

        # STEP 3: if bp_length is False:
        if not bp_length:
            records = []
            record_iter = SeqIO.parse(out_file, 'fasta')
            records = [SeqRecord.SeqRecord(record.seq,
                                           id=''.join(record.description.split()[0]),  # id=' '
                                           description='')
                       for record in record_iter]
            SeqIO.write(records, out_file, 'fasta')
        t1 = datetime.now()

        print("Running time: %s seconds" % (t1 - t0))
    return file_handles


def iqtree(in_path, out_path='./', bootstrap='UFBoot', bootstrap_number=1000, threads='AUTO',
           partition_file='', iteration_number=1000, constrain='', redo=False,
           additional_params='',
           pure_command_mode=False, pure_command=''):
    '''
    call iqtree to construct phylogenetic tree
    ----------
    Parameters
    - in_path - the input file(s) and folder(s). could be a list or a string.
    - out_path - destination folder of output, where log files and result will be written.
    - bootstrap - the algorithm of bootstrap [UFBoot, SNB] or empty string to cancle bootstrap.
        UFBoot for Ultra Fast Bootstrap, SNB for Standard Nonparamertric bootstrap
    - bootstrap_number - the number of bootstrap.
    - threads - the number of threads, default is 'AUTO'.
    - partition_file - the handle of a partition file in nexus. If not, automatically call model finder
    - iteration_number - the iteration number of bootstrap.
    - redo - if there is complete output result in the output path, set redo to rerun the analysis.
    - additional_params - additional parameters in the form of command.
    - pure_command_mode - if True, only run commands in the textbox, one command per line.
    - pure_command - run only if pure_command_mode is True, replace the GUI operations.
    -------
    Returns
    - file_handles - the valid input files if not in pure command mode
    - commands - the commands in pure command mode
    [] if path invalid"""
    '''

    # STEP 1: if pure command, then only execute input command
    # print("# STEP 1: if pure command, then only execute input command")
    if pure_command_mode == True:
        pure_command = pure_command.split('\n')
        for command in pure_command:
            os.system(command)
        return pure_command

    # STEP 2: check path validity and get file handles
    # print("# STEP 2: check path validity and get file handles")
    file_handles = get_file_handles(in_path)
    # print(file_handles)
    if not check_outpath_validity(out_path):
        return []

    # STEP 3: get parameters and call iqtree
    # print("# STEP 3: get parameters and call iqtree")
    for in_file in file_handles:
        basename = os.path.basename(in_file)
        out_file = os.path.join(out_path, basename)
        command = [f"iqtree -s {in_file}"]

        if bootstrap:
            if bootstrap == 'UFBoot':  # short for ultra fast boostrap (2)
                command.append(f"-bb {bootstrap_number}")
            elif bootstrap == 'SNB':  # short for standard nonparametric bootstrap
                command.append(f"-b {bootstrap_number}")

        if threads:
            command.append(f"-nt {threads}")

        if partition_file:
            command.append(f"-spp {partition_file}")

        if iteration_number:
            command.append(f"-nm {iteration_number}")

        if constrain:
            constrain = get_file_handles(constrain)
            if len(constrain) > 1:
                return "more than one constrain files found."
            command.append(f"-g {constrain[0]} -pre {constrain[0]}")

        if redo:
            command.append("-redo")

        if additional_params:
            command.append(str(additional_params))

        command = ' '.join(command)
        print(command)
        os.system(command)
    return file_handles
