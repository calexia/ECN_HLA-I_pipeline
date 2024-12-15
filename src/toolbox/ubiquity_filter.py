import pandas as pd

"""STEP 5 : separate non-ubiquitous peptides from ubiquitous peptides, based on protein atlas"""
"""STEP 5a : use the 3 dictionnaries to match the RNA tissue cat / RNA TS TPM / gene synonym to the gene """
def addGeneSyno(data, geneSynonymDict):
    print('UBI : Adding gene synonyms...')
    for index, row in data.iterrows():  # iterates over rows in df
        geneList = row["Gene name"].split(";")  # get every gene from gene list
        geneSynonym = []  # intializes an empty list
        for gene in geneList:  # consider each gene
            sepsyn = ','
            geneSynonymList = ''
            if gene in geneSynonymDict:  # check if gene is a key in gene synonym dict
                geneSynonym.append(str(geneSynonymDict[gene]))  # append the corresponding value to a list
                geneSynonymList = sepsyn.join(geneSynonym)  # join all the gene synonyms
            data.at[index, "Gene synonyms"] = geneSynonymList

    return data

def nonUbiExpression(data, RNATSTPMDict):
    print('UBI : Adding tissue and RNA TS TPM expression...')
    sep = ';'

    for index, row in data.iterrows(): #iterates over rows in df
        geneList = row["Gene name"].split(";") #get every gene from gene list
        geneSynoList = row["Gene synonyms"].split(',')
        TSTPM = [] #intializes an empty list

        for gene in geneList:
            if gene in RNATSTPMDict:
                TSTPM.append(RNATSTPMDict[gene]) #appends the corresponding value to a list
            elif gene not in RNATSTPMDict:
                for geneSyno in geneSynoList:
                    if geneSyno in RNATSTPMDict:
                        TSTPM.append(RNATSTPMDict[geneSyno])  # appends the corresponding value to a list

            TSTPMList = sep.join(TSTPM)  # join all the values into a list
            data.at[index, "RNA TS TPM"] = TSTPMList  # replace the list into the df, in a new col

    return data

