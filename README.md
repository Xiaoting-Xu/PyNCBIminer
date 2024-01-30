# PyNCBIminer

PyNCBIminer is a user-friendly graphical interface software designed for efficient and precise retrieval of GenBank data. Its simple operation does not require a background in bioinformatics. PyNCBIminer automatically performs BLAST iterations and optimizations based on user-provided nucleotide sequence names, target taxa, and initial reference sequences. This enhances the integrity of both sequence quantity and length, facilitating the identification and retrieval of specific gene sequence data for multiple species or taxa as specified by the user.

# 1. Download

**For windows**: Download the corresponding version of the packaged from(), users only need to unzip the package and run the.exe file directly.
**For macOS**:

# 2. Running PyNCBIminer

## 2.1 Sequence Retrieving Module

The Sequence Retrieving Module is instrumental in accomplishing the primary function of PyNCBIminer, enabling the identification and download of specific gene sequence data for multiple species or taxa as specified by the user.

### **Example**

**Taxa**: Saxifragales, **Gene**: ITS

### Step 1: Set target region

**Working diretory**: Enter the absolute path of the working directory in the blank box under **'Working directory'**. You can also browse and select by clicking the `View` button. 

**Basic setting**: Under the **'Basic Settings'** section, choose **'ITS'** from the dropdown options and click the `Set Target Region` button.


**Advanced setting**: After clicking the `Set Target Region` button, the **'Advanced Settings'** section will automatically configure default reference sequences and BLAST parameters, requiring no further modifications. If you wish to modify parameters, you can edit the text boxes; all text boxes are editable. For more details on parameters, please refer to the [manual]().

### Step 2: Submit BLAST

Click the `Submit New BLAST` button in the 'Working directory' section to initiate the BLAST process. The process continues until no new sequences can be found, at which point the BLAST stops, and sequence downloading begins.

### Step 3: Load previous job

If you wish to resume an incomplete task, enter the working directory and click the `Load Previous Job` button for loading. The program will automatically assess the progress and resume execution from the point of interruption.

**View results**

Three folders, namely `parameters`,`results`, and `tmp_files` , will be generated in the working directory. They primarily store initial parameters, BLAST results, and intermediate files, respectively.

`parameters <folder>`: The BLAST parameters (`blast_parameters.txt <file>`), the fasta file of initial reference sequences (`initial_queries.fasta <file>`), information about the newly selected reference sequences in each round (`all_new_queries_info.txt <file>`), and the sequences in fasta format (`ref_seq <folder>`) along with the sequence alignment results (`ref_msa <folder>`) are stored in the respective folders in the working directory.

`results`: The final BLAST results (`blast_results.txt <file>`), the annotated fasta file and sequence information for correctly annotated sequences (`blast_results_checked.fasta <file>` and `blast_results_checked_seq_info.txt <file>`), as well as the fasta file and sequence information for sequences with annotation errors (`erroneous_blast_results_checked.fasta <file>` and `erroneous_blast_results_checked_seq_info.txt <file>`), are all available in the output folders.

`tmp_files`: For each round of BLAST, intermediate results are stored in separate subfolders with a prefix 'BLAST_' followed by numerical identifiers, created within the working directory. These subfolders individually contain the results of each round of BLAST. For more details, please refer to the [manual]()



