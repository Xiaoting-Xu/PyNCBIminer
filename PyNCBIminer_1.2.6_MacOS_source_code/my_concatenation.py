# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:36
# @Author:Ruijing Cheng
# @File:my_concatenation.py
# @Software:PyCharm

from datetime import datetime
import os
from pathlib import Path
from format_wizard import taxon_completion, concat, fas2phy


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