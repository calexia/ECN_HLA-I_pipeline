import concatenation as concat
import pandas as pd

FILE_PATH = "D:\\DataAlexia\Projets\Scripts\HLA-I\Peptide_list_cleaning\Output\\"

file = pd.read_csv(FILE_PATH + "20211007_CleanPeptideList.txt", sep='\t', dtype="str")


file = concat.rawPeptide(file)

file.to_csv(FILE_PATH + "20211007_CleanPeptideList_rawPep.txt", sep='\t', index=False)
