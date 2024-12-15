import re
import pandas
import numpy as np


#clean the result table
def cleanResults(dfconcat):
    print("\nCONCAT : CLEANING RESULT TABLE...")


    cols = [
        "Gene name",
        "Gene synonyms",
        "Gene description",
        "Accession",
        "Length",
        "Peptide",
        "PTM",
        "HLA restriction",
        "Rank",
        "Total occurences",
        "Sample occurences",


    ]
    dfconcat = dfconcat[cols]
    dfconcat = dfconcat.apply(pandas.to_numeric, downcast='integer', errors="ignore")  #applies to_numeric to the df = changes type to integer
    dfconcat = dfconcat.sort_values(by=['Gene name', 'Length'])

    return dfconcat



def rawPeptide(finalTable):
    print("CONCAT : Adding Raw peptide column")
    poolValue = re.compile('[a-zA-Z]{1,}')

    for index, row in finalTable.iterrows():
        peptide = str(row["Peptide"])
        peptide = poolValue.findall(peptide)
        rawPeptide = ''.join(peptide)
        rawPeptide = str(rawPeptide)
        finalTable.at[index, "Raw_Peptide"] = str(rawPeptide)
    return finalTable
