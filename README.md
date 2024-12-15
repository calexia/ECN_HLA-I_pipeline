# ECN90 HLA-peptidomics data processing pipeline


## Description
This script processes HLA-peptidomics datasets. It filters out peptides based on their length (including only 8-14 mers) and assign HLA-restriction according to NetMHCpan 4.1a. It also filters out peptides that are not enriched in endocrine pancreas based on the comparison of the peptide source protein/gene to existing databases and datasets and related tissue enrichment information. 

## Prerequisites

* Linux environment
* Python v3
* NetMHCpan 4.1
To install netMHCpan please download and follow instructions on:  https://services.healthtech.dtu.dk/services/NetMHCpan-4.1a/


## Installation
Creating and setting up a python virtual environment (venv) should be done in the folder that contains the requirements.txt

Create a venv
````bash
python -m venv venv
````
Assign the venv to a terminal
````bash
source venv/bin/activate #Linux
````
Install dependencies
````bash
pip install -r requirements.txt
````

## Instructions for use
Several pipelines (MAIN files) are available. 
* ECN90 basal vs IFNa  
This pipeline takes as an input the output files from Peaks (DB search psm) of basal and IFN-alpha exposed ECN90 cells. It concatenates the information so that the output contains only one row per identified peptide. The final output also informs about the number of biological replicates (n=4 per condition) and the number of technical replicates (n=2 per biological replicate) among which the peptide was found. The HLA-binding is predicted for HLA-A02:01, -A03:01, -B40:01, -B49:01, -C03:04, -C07:01, -E01:01.  
 
```bash
python3 ./src/IFNa/HLA-I_IFNa_MAIN_NatureComm.py <basal_file_path> <ifna_file_path> <output_folder_path>
```

* human islet  
This pipeline takes as an input the output files from Peaks (DB search psm) of each human islet dataset. The islet sample (1 or 2) must be specified.

```bash
python3 ./src/Human_islets/HLA-I_islets_MAIN_NatureComm.py <islet_file_path> <islet_list_selection> <output_folder_path>
```

## Demo

From the folder of this README, run the following commands to test proper installation of the scripts.
* ECN90 basal vs IFNa
```bash
python3 ./src/IFNa/HLA-I_IFNa_MAIN.py ./demo/ECN90_basal_vs_IFNa/Input/DB_search_psm_Basal.csv ./demo/ECN90_basal_vs_IFNa/Input/DB_search_psm_IFNa.csv ./demo/ECN90_basal_vs_IFNa/Output/;
diff ./demo/ECN90_basal_vs_IFNa/Expected_output/CleanResultTable.txt ./demo/ECN90_basal_vs_IFNa/Output/CleanResultTable.txt
```

* Human islet

Islet 1:
```bash
python3 ./src/Human_islets/HLA-I_islets_MAIN.py ./demo/Human_islets/Input/DB_search_psm_islet1.csv 1 ./demo/Human_islets/Output/Islet1/;
diff ./demo/Human_islets/Expected_output/islet1/CleanResultTable.txt ./demo/Human_islets/Output/Islet1/CleanResultTable.txt
```

Islet 2:
```bash
python3 ./src/Human_islets/HLA-I_islets_MAIN.py ./demo/Human_islets/Input/DB_search_psm_islet2.csv 2 ./demo/Human_islets/Output/Islet2/;
diff ./demo/Human_islets/Expected_output/islet2/CleanResultTable.txt ./demo/Human_islets/Output/Islet2/CleanResultTable.txt
```