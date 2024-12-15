import pandas as pd
import os
import re


def multipleTable (finalTable):
    print("MHC : Separating data into multiple tables")
    # Create a list with 7 dataframes for peptides lengths from 8 to 14
    dfLengths = [pd.DataFrame() for i in range(7)]

    for index, row in finalTable.iterrows():
        pepLength = int(row["Length"])
        dfLengths[pepLength -8] = dfLengths[pepLength - 8].append(row)

    return dfLengths

DEFAULT_DEBUG_NETMHC_OUTPUT_PATH = '/mnt/109E1E2A9E1E08BC/DataAlexia/Manuscripts/Research/IFNalpha/Revision_round1/Code_availability/data/example/Output/'
def computeBindingPredictions(finalTable, filePath, type, allelesList, outputFilePath= DEFAULT_DEBUG_NETMHC_OUTPUT_PATH):
    print("MHC : Initializing predictions")
    dfLengths = multipleTable(finalTable)
    allelesStr= ','.join(allelesList)

    for idx, dfLength in enumerate(dfLengths):
        dfLengthNoDupli = dfLength.drop_duplicates(subset='Raw_Peptide', keep='first')
        nmersFile = filePath + '{}mers.txt'.format(8 + idx)
        dfLengthNoDupli.to_csv(path_or_buf=nmersFile, columns=['Raw_Peptide'], index=False, header=None, sep='\t')
        #Command line for NetMHCpan for IFNa dataset
        outputFile = outputFilePath + "{}mersPredictions{}.txt".format(8 + idx, type)
        netMHCCmd = "netMHCpan4.1a -a {} -f {} -p -xls -xlsfile {}".format(allelesStr, nmersFile, outputFile)
        print("Outputing netMHCpan in {}".format(outputFile))
        print("NetMHCCmd is {}".format(netMHCCmd))
        os.system(netMHCCmd)


def bindingPredictionsTreatment(predictiondf, allelesList):
    """

    :param predictiondf:
    :return:
    """
    print ("MHC : Processing NetMHC data...")
    restrictiondf = pd.DataFrame()

    for index, row in predictiondf.iterrows():
        NBValue = row['NB']

        restrictiondf = restrictiondf.append(row)

        colNameList = ["Rank_" + s for s in allelesList]

        if NBValue != 0:

            scores = list(map(lambda colName: row[colName], colNameList))
            # is equivalent to
            # scoreAllele1 = row[colNameList[0]]
            # scoreAllele2 = row[colNameList[1]]
            # scoreAllele3 = row[colNameList[2]]
            # scoreAllele4 = row[colNameList[3]]
            # scoreAllele5 = row[colNameList[4]]
            # scoreAllele6 = row[colNameList[5]]
            # scoreAllele7 = row[colNameList[6]]
            # ...
            #
            # scores = [scoreAllele1, scoreAllele2, scoreAllele3, scoreAllele4, scoreAllele5, scoreAllele6, scoreAllele7, ...]

            #creates a dict for indexing the alleles such as index:allele
            allelesDict = { i : allelesList[i].split('-')[1] for i in range (0, len(allelesList)) }

        if NBValue == 1:
            minScore = min(scores)
            alleleIndex = scores.index(minScore)
            restrictiondf.at[index, "HLA restriction"] = allelesDict[alleleIndex]

        elif NBValue >= 2: #considers all cases
            lowest = secondLowest = 3.0
            for float in scores:
                if float < lowest:
                    secondLowest = lowest
                    lowest = float

                elif secondLowest != lowest and secondLowest>float:
                    secondLowest = float

            if lowest in scores and secondLowest in scores:
                if lowest == secondLowest:
                    allelesIndexes = list(idx for idx,value in enumerate(scores) if value == lowest)
                    alleleIndexLowest = allelesIndexes[0]
                    alleleIndexSecondLowest = allelesIndexes[1]
                else :
                    alleleIndexLowest = scores.index(lowest)
                    alleleIndexSecondLowest = scores.index(secondLowest)

                if round((3.0*lowest),4) < secondLowest:
                    restrictiondf.at[index, "HLA restriction"] = allelesDict[alleleIndexLowest]
                else:
                    if allelesDict[alleleIndexLowest] == "E01:01":
                        restrictiondf.at[index, "HLA restriction"] = allelesDict[alleleIndexSecondLowest]
                    elif allelesDict[alleleIndexSecondLowest] == "E01:01":
                        restrictiondf.at[index, "HLA restriction"] = allelesDict[alleleIndexLowest]
                    else:
                        restrictiondf.at[index, "HLA restriction"] = allelesDict[alleleIndexLowest] + '/' + allelesDict[alleleIndexSecondLowest]

    return restrictiondf


def finalTableHLA_Assignment(finalTable, concatRestrictiondf, filePath):
    finalTableAssigned = finalTable.set_index("Raw_Peptide").join(concatRestrictiondf.set_index("Peptide"), how='left')
    finalTableAssigned = finalTableAssigned.reset_index()
    #finalTableAssigned.to_csv(path_or_buf=filePath +"assigned.txt", sep='\t' )

    for index, row in finalTableAssigned.iterrows():
        HLARestrictions = row["HLA restriction"]
        HLARestrictions = str(HLARestrictions).split('/')

        #Thalfs = []
        Ranks = []

        for restriction in HLARestrictions:

            if str(restriction) != "nan":

                #Thalf = str(row["Thalf(h)_HLA-" + str(restriction)])
                Rank = str(row["Rank_HLA-" + str(restriction)])

                #Thalfs.append(Thalf)
                Ranks.append(Rank)

            else:
                restriction = ''

        #finalTableAssigned.at[index, "Thalf(h)"] = "/".join(Thalfs)
        finalTableAssigned.at[index, "Rank"] = "/".join(Ranks)


    colsToKeep =[
        #"Previously validated ?",
        "Peptide",
        "Length",
        "Accession",
        "Fraction",
        "PTM",
        "HLA restriction",
        "Rank"
    ]

    #finalTableAssigned = finalTableAssigned[colsToKeep]

    return finalTableAssigned
