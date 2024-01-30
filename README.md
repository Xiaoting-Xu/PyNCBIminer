<img width="273" alt="249798617c0d2757f5d202764f05730" src="https://github.com/Xiaoting-Xu/PyNCBIminer/assets/158121092/a3a5b2b5-bf03-43fc-8078-624add350cd0"># PyNCBIminer

PyNCBIminer is a user-friendly graphical interface software designed for efficient and precise retrieval of GenBank data. Its simple operation does not require a background in bioinformatics. PyNCBIminer automatically performs BLAST iterations and optimizations based on user-provided nucleotide sequence names, target taxa, and initial reference sequences. This enhances the integrity of both sequence quantity and length, facilitating the identification and retrieval of specific gene sequence data for multiple species or taxa as specified by the user.

# 1. Download

**For windows**: Download the corresponding version of the packaged from(), users only need to unzip the package and run the.exe file directly.
**For macOS**:

# 2. Running PyNCBIminer

## 2.1 Sequence Retrieving Module

The Sequence Retrieving Module is instrumental in accomplishing the primary function of PyNCBIminer, enabling the identification and download of specific gene sequence data for multiple species or taxa as specified by the user.

### **Example**

**Taxa**: Saxifragales, **Gene**: ITS

### Step 1: set working directory

Enter the absolute path of the working directory in the blank box under 'Working directory.' You can also browse and select by clicking the 'View' button. Ensure that the correct input is the directory path, not the file path. The path should not contain any spaces. If the input is a non-existent directory path, a new folder will be created.

### Step 2: set target region

Under the 'Basic Settings' section, choose 'ITS' from the dropdown options and click the 'Set Target Region' button.

After clicking the 'Set Target Region' button, the 'Advanced Settings' section will automatically configure default reference sequences and BLAST parameters, requiring no further modifications. If you wish to modify parameters, you can edit the text boxes; all text boxes are editable. For more details on parameters, please refer to the [manual]().

