# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:57
# @Author:Ruijing Cheng
# @File:blast_put_get.py
# @Software:PyCharm

import os
import re
import time
import pandas as pd
import numpy as np
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request
from pathlib import Path
import func_timeout.exceptions
from func_timeout import func_set_timeout
from Bio import SeqIO
from tools import get_query_accession


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
