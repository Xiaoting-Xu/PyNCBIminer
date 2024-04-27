# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:55
# @Author:Ruijing Cheng
# @File:tools.py
# @Software:PyCharm


def print_line(character="#"):
    print(character * 50)


def get_query_accession(record):
    parts = record.description.split(" ")[0].split("|")[0].split(":")
    # assert len(parts) == 2
    return parts[0]