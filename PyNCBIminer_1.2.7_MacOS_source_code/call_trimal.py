import os
from Bio import SeqIO, SeqRecord
from format_wizard import check_inpath_validity, check_outpath_validity, get_file_handles, create_folder
from datetime import datetime


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
    if len(file_handles) == 0:
        print("Could not find any fasta file in the input.")
        return
    if not check_outpath_validity(out_path):
        return []
    
    # STEP 2: get parameters and call trimal
    for in_file in file_handles:
        print("Trimming %s..." % in_file)
        t0 = datetime.now()
        basename = os.path.basename(in_file)
        out_file = os.path.join(out_path, basename)
        command = [f"trimal -in {in_file} -out {out_file}"]
        
        if htmlout:
            html_folder = os.path.join(out_path, 'htmlout')
            create_folder(html_folder)
            html_out_file = os.path.join(html_folder, os.path.splitext(basename)[0]+'.html')
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
        os.system(command)
        
        # STEP 3: if bp_length is False:
        if not bp_length:
            records = []
            record_iter = SeqIO.parse(out_file, 'fasta')
            records = [SeqRecord.SeqRecord(record.seq, 
                                           # id=' '.join(record.description.split()[0]),
                                           id=''.join(record.description.split()[0]),
                                           description='')
                       for record in record_iter]
            SeqIO.write(records, out_file, 'fasta')
        t1 = datetime.now()
        print("Running time: %s seconds" % (t1 - t0))

    return file_handles
