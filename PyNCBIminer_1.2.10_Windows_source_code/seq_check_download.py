# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:50
# @Author:Ruijing Cheng
# @File:seq_check_download.py
# @Software:PyCharm

import urllib.error
import func_timeout.exceptions
from func_timeout import func_set_timeout
from Bio import Entrez
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
import os
import pandas as pd
from tools import get_query_accession


@func_set_timeout(600)
def my_efetch(accession, strand, seq_start, seq_stop):
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession,
                           strand=strand, seq_start=seq_start, seq_stop=seq_stop)
    print("Extended start: %d, " % seq_start, end="")
    print("Extended end: %d, " % seq_stop, end="")
    print("Strand: %d, " % strand, end="")
    return handle


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