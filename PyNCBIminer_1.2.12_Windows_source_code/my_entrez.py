# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:46
# @Author:Ruijing Cheng
# @File:my_entrez.py
# @Software:PyCharm

from func_timeout import func_set_timeout
from Bio import Entrez


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
def entrez_count(entrez_email, organisms, entrez_qualifier="", date_from="", date_to=""):
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
    else:
        print("Entrez email: %s" % entrez_email)
    Entrez.email = entrez_email
    entrez_query=format_entrez_query(organisms, entrez_qualifier, date_from, date_to)
    print("Entrez query: %s" % entrez_query)
    handle = Entrez.esearch(db="nucleotide", term=entrez_query)
    record = Entrez.read(handle)
    print("Entrez search results count: %s" % record["Count"])
    return int(record["Count"])


def entrez_summary(entrez_email, organisms, target_region_dict, date_from="", date_to=""):
    """
    gives the user a summary of the widely used markers for their focal taxa
    :return:
    """
    if len(entrez_email) == 0:
        print("Warning: email address is not specified.")
        print("To make use of NCBI's E-utilities, NCBI requires you to specify your email address with each request. ")
        print("In case of excessive usage of the E-utilities, "
              "NCBI will attempt to contact a user at the email address provided "
              "before blocking access to the E-utilities.")
    else:
        print("Entrez email: %s" % entrez_email)
    Entrez.email = entrez_email
    print("PyNCBIminer will search for the approximate number of the following markers. \n"
          "It will take a few minutes. \n"
          "Please wait until it notifies you that the search is complete.")
    for target_region in target_region_dict.keys():
        print(target_region, end="\t")
        entrez_query = format_entrez_query(organisms, target_region_dict[target_region], date_from, date_to)
        # print(entrez_query)
        handle = Entrez.esearch(db="nucleotide", term=entrez_query)
        record = Entrez.read(handle)
        print(record["Count"])
    print("The search is complete.")
    # print("The most widely used marker is: ")
