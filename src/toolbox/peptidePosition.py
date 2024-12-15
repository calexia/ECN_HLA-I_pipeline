import pandas as pd
import src.toolbox.concatenation as concat
import re
import numpy as np
from Bio import SeqIO


def addPeptidePosition(data, DB_NAME):
    print("Adding peptide position")

    # data = pd.read_csv(FILE_NAME,
    #                    sep="\t",
    #                    dtype={
    #                        "Peptide": "str",
    #                        "Length": np.int64,
    #                        "Accession": "str"
    #                    }
    #                    )

    """ ADD RAW PEPTIDE INFORMATION TO THE DATA TABLE """
    data_rawPep = concat.rawPeptide(data)

    """ Iteration over raw peptide and adds start and end"""
    regex = "[A-Z]{1}\w{5}"
    database = SeqIO.parse(DB_NAME, "fasta")
    databaseDict = SeqIO.to_dict(database, key_function=lambda x: x.id.split('|')[1])

    for index, row in data_rawPep.iterrows():
        rawPep = row["Raw_Peptide"]
        accessionID = row["Accession"]
        pepLength = row["Length"]
        firstAccession = re.findall(regex, accessionID)[0]

        print("Checking peptide : " + rawPep)

        if firstAccession in databaseDict:
            proteinSeq = databaseDict[firstAccession].seq
            indexStart = int(proteinSeq.find(rawPep))

            if indexStart != -1:
                startPos = str(indexStart + 1)
                endPos = str(indexStart + pepLength)
                data_rawPep.at[index, "Position"] = startPos + '-' + endPos

            else:
                print("Peptide not in this protein")
        else:
            print("Accession not in database")

    cols = [
        "Gene description",
        "Gene synonyms",
        "Gene name",
        "Accession",
        "Position",
        "Peptide",
        "Length",
        "PTM",
        "Presence",
        "HLA restriction",
        "Rank",
        "FractionBasal",
        "FractionIFNa",
        "Total occurencesBasal",
        "Sample occurencesBasal",
        "Total occurencesIFNa",
        "Sample occurencesIFNa",
        "Comments"
    ]

    data_rawPep = data_rawPep[cols]
    return data_rawPep

# FILE_PATH = "D://DataAlexia/Projets/Scripts/HLA-I/IFNa_dataset/"
# FILE_NAME = FILE_PATH + "Output/CleanResultTable_NoTERT_specificPTMadded.txt"
# DB_NAME = FILE_PATH + "20190122_SPhs_AC_mRNAvariants_IFNalpha_ERV_SV40_validated.fasta"
# data=pd.read_csv(FILE_NAME, sep='\t')
#
# data_rawPep = addPeptidePosition(data, DB_NAME)
# data_rawPep.to_csv(path_or_buf=FILE_PATH + 'Output/CleanPeptideList.txt', index=False, sep='\t')


