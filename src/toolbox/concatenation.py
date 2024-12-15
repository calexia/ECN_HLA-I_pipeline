import re
import pandas
import numpy as np

#Concatenate the result tables
def concatResults(dfbasal, dfifna):
    print ("\nCONCAT : CONCATENATING THE RESULT TABLE...")
    dfconcat = dfbasal.set_index("Peptide").join(dfifna.set_index("Peptide"), how='outer', lsuffix="Basal",
                                                 rsuffix="IFNa")
    # dfconcat = pandas.merge(dfbasal, dfifna, on='Peptide', how='outer', suffixes=("Basal", "IFNa"))
    dfconcat = dfconcat.reset_index()

    # Remove NaN added by merging and recomputing dtypes
    dfconcat = dfconcat.fillna("0")
    dfconcat = dfconcat.apply(pandas.to_numeric, downcast='integer', errors="ignore") #applies to_numeric to the df = changes type to integer

    for index, row in dfconcat.iterrows():
        occBasal = int(row["Total occurencesBasal"])
        occIFNa = int(row["Total occurencesIFNa"])
        if (occBasal !=0) & (occIFNa !=0):
            dfconcat.at[index,"Presence"] = "Both"
        elif (occBasal!=0) & (occIFNa==0):
            dfconcat.at[index, "Presence"] = "Basal"
        elif (occBasal==0) & (occIFNa!=0):
            dfconcat.at[index, "Presence"] = "IFNa"

    return dfconcat

#clean the result table
def cleanResults(dfconcat):
    print("\nCONCAT : CLEANING RESULT TABLE...")
    poolValueAll = re.compile('[\w*\s\(\)]{2,}')
    poolValuePTM = re.compile('[\w*\s\(\)\-]{1,}')

    # removes the 0 added in the previous function
    dfconcat.replace("0", np.nan, inplace=True)

    #initiates empty col...
    dfconcat["PTM"] = ""
    dfconcat["Accession"] = ""
    dfconcat["Gene name"] = ""
    dfconcat["Gene description"] = ""
    dfconcat["Gene synonyms"] = ""
    dfconcat["HLA restriction"] = ""
    dfconcat["Rank"] = ""

    #... with specific dtype
    dfconcat.astype(
        {
            "PTM" : 'str',
            'Accession': 'str',
            'Gene name': 'str',
            'Gene description': 'str',
            'Gene synonyms': 'str',
            'HLA restriction': 'str',
            'Rank': 'str',
        }
    )

    for index, row in dfconcat.iterrows():
        lengthBasal = row["LengthBasal"]
        lengthIFNa = int(row["LengthIFNa"])
        if lengthBasal !=0:
            dfconcat.at[index, "Length"] = lengthBasal
        else:
            dfconcat.at[index, "Length"] = lengthIFNa

        PTMBasal = row["PTMBasal"]
        PTMIFNa = row["PTMIFNa"]

        accBasal = row["AccessionBasal"]
        accIFNa = row["AccessionIFNa"]

        geneBasal = row["Gene nameBasal"]
        geneIFNa = row["Gene nameIFNa"]

        genedescriBasal = row["Gene descriptionBasal"]
        genedescriIFNa = row["Gene descriptionIFNa"]

        geneSynonymsBasal = row["Gene synonymsBasal"]
        geneSynonymsIFNa = row["Gene synonymsIFNa"]

        HLA_restrictionBasal = row["HLA restrictionBasal"]
        HLA_restrictionIFNa = row["HLA restrictionIFNa"]

        RankBasal = row["RankBasal"]
        RankIFNa = row["RankIFNa"]

        presence = row["Presence"]


        if presence == "IFNa":
            # Refactor this as function1
            dfconcat.at[index, "PTM"] = PTMIFNa
            dfconcat.at[index, "Accession"] = accIFNa
            dfconcat.at[index, "Gene name"] = geneIFNa
            dfconcat.at[index, "Gene description"] = genedescriIFNa
            dfconcat.at[index, "Gene synonyms"] = geneSynonymsIFNa
            dfconcat.at[index, "HLA restriction"] = HLA_restrictionIFNa
            dfconcat.at[index, "Rank"] = str(RankIFNa)
            ## Refactor this as function1
        elif presence == "Basal":
            # Refactor this as function1
            dfconcat.at[index, "PTM"] = PTMBasal
            dfconcat.at[index, "Accession"] = accBasal
            dfconcat.at[index, "Gene name"] = geneBasal
            dfconcat.at[index, "Gene description"] = genedescriBasal
            dfconcat.at[index, "Gene synonyms"] = geneSynonymsBasal
            dfconcat.at[index, "HLA restriction"] = HLA_restrictionBasal
            dfconcat.at[index, "Rank"] = str(RankBasal)
            ## Refactor this as function1
        elif presence == "Both":
            # Refactor this as function2(regex)
            PTMstr = str(PTMBasal) + ";" + str(PTMIFNa)
            PTMList = list(PTMstr.split(';'))
            PTMListUnique = str(set(PTMList))
            PTMListUnique = poolValuePTM.findall(PTMListUnique)
            truePTMListUnique = ';'.join(PTMListUnique)
            ## Refactor this as function2

            # Refactor this as function2(regex)
            accessionstr = accBasal + ';' + accBasal
            accessionList = list(accessionstr.split(';'))
            accessionListUnique = str(set(accessionList))
            accessionListUnique = poolValueAll.findall(accessionListUnique)
            trueAccessionListUnique = ';'.join(accessionListUnique)
            ## Refactor this as function2

            # Refactor this as function2(regex)
            genestr = geneBasal + ';' + geneBasal
            geneList = list(genestr.split(';'))
            geneListUnique = str(set(geneList))
            geneListUnique = poolValueAll.findall(geneListUnique)
            trueGeneListUnique = ';'.join(geneListUnique)
            ## Refactor this as function2

            # Refactor this as function2(regex)
            genedesristr = genedescriBasal + ';' + genedescriIFNa
            genedescriList = list(genedesristr.split(';'))
            genedescriListUnique = str(set(genedescriList))
            genedescriListUnique = poolValueAll.findall(genedescriListUnique)
            trueGeneDescriListUnique = ';'.join(genedescriListUnique)
            ## Refactor this as function2

            # Refactor this as function2(regex)
            geneSynostr = geneSynonymsBasal + "," + geneSynonymsIFNa
            geneSynoList = list(geneSynostr.split(','))
            geneSynoListUnique = str(sorted(set(geneSynoList)))
            geneSynoListUnique = poolValueAll.findall(geneSynoListUnique)
            trueGeneSynoListUnique = ','.join(geneSynoListUnique)
            ## Refactor this as function2


            # Refactor this as function1
            dfconcat.at[index, "PTM"] = truePTMListUnique
            dfconcat.at[index, "Accession"] = trueAccessionListUnique
            dfconcat.at[index, "Gene name"] = trueGeneListUnique
            dfconcat.at[index, "Gene description"] = trueGeneDescriListUnique
            dfconcat.at[index, "Gene synonyms"] = trueGeneSynoListUnique
            dfconcat.at[index, "HLA restriction"] = HLA_restrictionBasal
            dfconcat.at[index, "Rank"] = RankBasal
            ## Refactor this as function1

    cols = [
        "Gene name",
        "Gene synonyms",
        "Gene description",
        "Accession",
        "Length",
        "Peptide",
        "PTM",
        "Presence",
        "HLA restriction",
        "Rank",
        "FractionBasal",
        "FractionIFNa",
        "Total occurencesBasal",
        "Sample occurencesBasal",
        "Total occurencesIFNa",
        "Sample occurencesIFNa"
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
