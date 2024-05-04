import os
import pandas as pd
from pathlib import Path

path0 = r"./PyNCBIminer_1.2.5_Windows_source_code/TNRS_dep"
os.environ["R_HOME"] = path0
# 创建临时环境变量，path0应指向解压后的TNRS_dep目录，可以把TNRS_dep目录放在pyNCBIminer程序所在目录下
print(os.environ["R_HOME"])  # 输出验证环境变量目录是否正确

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def TNRS(path1):
    pandas2ri.activate()
    tnrs = importr("TNRS")
    # path1 = '你的路径/你的marker/results' # 用来获取当前marker的results目录，可以和主程序对接
    os.chdir(Path(path1)/Path("results"))  # 这里方便起见使用相对路径，需要后续对接主程序修改
    # 核心就是读取当前marker的results目录下的blast_results_checked_seq_info.txt
    df = pd.read_table('blast_results_checked_seq_info.txt')
    df2 = df[['accession', 'organism']]
    result = tnrs.TNRS(df2, sources="wcvp", classification="wfo", mode="resolve", matches="best",
                       skip_internet_check=True)
    with localconverter(ro.default_converter + pandas2ri.converter):
        pdf = ro.conversion.rpy2py(result)
    dfc = pd.merge(df, pdf, left_on='organism', right_on='Name_submitted', how='left')
    columns_to_keep = list(range(0, 22)) + [26, 33, 35, 54]
    dfc = dfc.iloc[:, columns_to_keep]
    dfc.rename(columns={'organism': 'organism_ori', 'Name_matched': 'organism'}, inplace=True)
    dfc.fillna('NA', inplace=True)
    dfc.to_csv('blast_results_checked_seq_info_pro.txt', index=False, sep='\t')
    print("操作完成", "处理完成,结果已保存在blast_results_checked_seq_info_pro.txt")
    # 将TNRS更正完成后形成的blast_results_checked_seq_info_pro.txt保存在当前marker的results目录下
    """关于blast_results_checked_seq_info_pro.txt：在原来blast_results_checked_seq_info.txt的基础上追加了Name_matched,
    Name_matched_accepted_family,Genus_matched,Taxonomic_status四列,分别对应更正后名称,对应科名,对应属名,名称状态。
    处理过程中已将Name_matched更名为organism,取代了原本的organism,原本的organism更名为organism_ori,方便后续filter操作
    建议在filter过程中加入对Taxonomic_status的考量
    后续filter操作就换成读取blast_results_checked_seq_info_pro.txt，其他不变"""
