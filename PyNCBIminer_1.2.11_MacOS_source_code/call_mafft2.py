import os
from format_wizard import check_inpath_validity, check_outpath_validity, get_file_handles, create_folder
from pathlib import Path
from Bio import SeqIO
from datetime import datetime


def mafft_add(in_path, in_file, out_path, cmd_str):
    print("Aligning %s..." % Path(in_path) / Path(in_file))
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
        # print("Aligned results: %s" % msa3)
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
            print("Aligning %s..." % in_file)
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
