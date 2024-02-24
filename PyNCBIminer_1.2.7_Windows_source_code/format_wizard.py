import os
import shutil
import sys
import pandas as pd
from Bio import SeqIO, SeqRecord


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
    file_handles = [x for x in file_handles if x.split(".")[-1] in ["fasta", "fas", "fa"]]

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
        temp_output_folder = os.path.join(out_path, 'filtered.result')
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
    return marker_record
