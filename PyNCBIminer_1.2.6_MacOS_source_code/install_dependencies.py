# *-* coding:utf-8 *-*
# @Time:2024/2/2 14:15
# @Author:Ruijing Cheng
# @File:install_dependencies.py
# @Software:PyCharm

import os
from pathlib import Path
import zipfile
from urllib.request import urlretrieve


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
