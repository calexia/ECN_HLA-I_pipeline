import re

"add the normal tissue information"
def addNormalTissue(data, normalTissueDict):
    print ('ENRICH : add normal tissue enrichment information')
    sep = ','

    data["Normal Tissue Expression"] = ""
    for index, row in data.iterrows():
        geneList = row["Gene name"].split(";")  # get every gene from gene list
        geneSynoList = row["Gene synonyms"].split(',')
        normalTissueExp = []

        for gene in geneList:
            if gene in normalTissueDict:
                normalTissueExp.extend(normalTissueDict[gene])

        if len(normalTissueExp) == 0:
            for geneSyno in geneSynoList:
                if geneSyno in normalTissueDict:
                    normalTissueExp.extend(normalTissueDict[geneSyno])

        normalTissueExpList = sep.join(normalTissueExp)
        data.at[index, "Normal Tissue Expression"] = normalTissueExpList

    return data


"checks the RNA evidence"
def RNAevidence(data):

    for index, row in data.iterrows():
        RNATSTPM = row["RNA TS TPM"]

        if "pancreas" not in RNATSTPM:
            data.at[index, "RNA evidence"] = "N"

        else:
            data.at[index, "RNA evidence"] = "Y"

    return data

"checks protein evidence"
def protEvidence(data):
    for index, row in data.iterrows():

        NormalTissueExp = row["Normal Tissue Expression"]

        if "pancreas" not in NormalTissueExp:
            data.at[index, "Prot evidence"] = "N"

        else:
            data.at[index, "Prot evidence"] = "Y"

    return data

"checks literature evidence"
def litEvidence(data):
    for index, row in data.iterrows():

        expTerm = row["Expression_term"]

        if "Pancreas" not in expTerm \
                and "Islets of Langerhans" not in expTerm \
                and 'Beta cell' not in expTerm:

            data.at[index, "Lit evidence"] = "N"

        else:
            data.at[index, "Lit evidence"] = "Y"

    return data

"filers RNA, prot and literature evidence"

def RNAprotLitEvidenceFilter(data, ubidf, type):
    indexesToDrop=[]
    for index, row in data.iterrows():
        RNAevidence = row["RNA evidence"]
        protEvidence = row["Prot evidence"]
        litEvidence = row["Lit evidence"]

        if RNAevidence == protEvidence == litEvidence == "N":
            ubidf = ubidf.append(row)
            indexesToDrop.append(index)
        else:
            continue

    data = data.drop(index = indexesToDrop)


    uniqueValuesEnriched = data["Peptide"].nunique()
    uniqueValuesNonEnriched = ubidf["Peptide"].nunique()

    print("The {} table contains {} unique peptides that are enriched based on protein atlas and {} unique ubiquitous peptides."
          .format(type, uniqueValuesEnriched, uniqueValuesNonEnriched))


    return data, ubidf


"""STEP 6 : identify beta-cell enriched expressed peptides, based on protein atlas and HPRD"""
"""STEP 6a : use the dictionary to match the expression_term to the gene in the result table"""
def expTermExpression(data, expTermDict):
    print ('ENRICH : Adding expression term...')

    sep = ";"
    expTermDictWithLists = { key: value.split(';') for key, value in expTermDict.items() }
    for index, row in data.iterrows():
        geneList = row["Gene name"].split(';')
        geneSynonymsList = row["Gene synonyms"].split(',')
        expTerm = []
        expTermSynonym = []
        for gene in geneList:
            if gene in expTermDictWithLists:
                expTermStr = sep.join(expTermDictWithLists[gene])
                expTerm.append(str(expTermStr))
                expTermList = sep.join(expTerm)
                data.at[index, "Expression_term"] = expTermList
            elif gene not in expTermDictWithLists:
                for synonym in geneSynonymsList:
                    if synonym in expTermDictWithLists:
                        expTermStr = sep.join(expTermDictWithLists[synonym])
                        expTermSynonym.append(str(expTermStr))
                        expTermSynonymList=sep.join(expTermSynonym)
                        data.at[index, "Expression_term"] = expTermSynonymList


    return data

"""STEP 6b : reformating RNA TS TPM : changes the "_" to ";", remove all the "." and "," / remove all the numbers
                AND remove duplicate values"""
def reformatRNATSTPM_expTerm(data):
    print ('ENRICH : Reformatting RNA TS TPM and expression term columns...')
    poolValues = re.compile('[a-zA-Z\s]{2,}') #regex for any word, space included at least 2 characters

    for index, row in data.iterrows():
        RNATSTPM = str(row["RNA TS TPM"])
        RNATSTPM = poolValues.findall(RNATSTPM)
        RNATSTPMlist = ';'.join(RNATSTPM)
        RNATSTPMlist = list(RNATSTPMlist.split(';'))
        RNATSTPMstr = str(list(sorted(set(RNATSTPMlist))))
        newRNATSTPM = poolValues.findall(RNATSTPMstr)
        truenewRNATSTPM = ';'.join(newRNATSTPM)
        data.at[index, "RNA TS TPM"] = truenewRNATSTPM

        expression_term = str(row["Expression_term"])
        expression_term = poolValues.findall(expression_term)
        exp_termList = ';'.join(expression_term)
        exp_termList = list(exp_termList.split(';'))
        exp_termstr = str(list(sorted(set(exp_termList))))
        newExp_term = poolValues.findall(exp_termstr)
        truenewExp_term = ';'.join(newExp_term)
        data.at[index, "Expression_term"] = truenewExp_term

    return data


""" STEP 6e : filter based on Sandberg's dataset, that has been reprocessed"""
def filterSandberg(data, ubidf, sandbergdf, type):
    print("ENRICH : Filtering genes on Sandberg's data...")

    genestoDrop = [] #initiates an empty list

    for index, row in sandbergdf.iterrows():
        enrichment = row["Cell enrichment"]
        gene = row["Genes"]

        if "not enriched" in enrichment or \
            "acinar cell" in enrichment or \
            "ductal cell" in enrichment or \
            "epithelial cells" in enrichment:
            genestoDrop.append(index)
        if "SST" in gene:
            genestoDrop.append(index)
    sandbergdf = sandbergdf.drop(index = genestoDrop)

    genesValidated = list(sandbergdf["Genes"])
    genesValidated.append('HSPA5') #even if HSPA5 is ubiquitous, it is of interest so manually included in valid genes

    indexesToDrop = []
    for idx, line in data.iterrows():
        geneList = line["Gene name"].split(';')
        geneSynonymsList = line["Gene synonyms"].split(',')


        if all (gene not in genesValidated for gene in geneList) and all (synonym not in genesValidated for synonym in geneSynonymsList):
            ubidf = ubidf.append(line)
            indexesToDrop.append(idx)

    data = data.drop(index=indexesToDrop)

    uniqueValuesEnriched = data["Peptide"].nunique()
    uniqueValuesNonEnriched = ubidf["Peptide"].nunique()

    print("The {} table contains {} unique peptides that are beta-cell enriched and {} unique ubiquitous peptides."
          .format(type, uniqueValuesEnriched, uniqueValuesNonEnriched))

    return data, ubidf


"""STEP 6f : add gene description using the corresponding dictionary"""
def addGeneDescri(data, geneDescriptionDict):
    print("ENRICH : Adding gene description...")
    for index, row in data.iterrows():
        geneList = row["Gene name"].split(';')
        descri = []
        descriList = []
        for gene in geneList:
            sep = ';'
            if gene in geneDescriptionDict:
                descri.append(str(geneDescriptionDict[gene]))
            descriList = sep.join(descri)
        data.at[index, "Gene description"] = descriList
    return data


def removeATHM(data, ATHMdf):
    ATHMindexes= []
    for index, row in data.iterrows():
        geneDescri = row["Gene description"]
        if "Myosin" in geneDescri \
                or "myosin" in geneDescri \
                or "tubulin" in geneDescri \
                or "Tubulin" in geneDescri \
                or "Actin" in geneDescri \
                or "Histone"in geneDescri \
                or "histone" in geneDescri \
                or "Telomerase reverse transcriptase" in geneDescri:
            # cannot consider actin without capital letter because all prot inter-actin-g with others will be discarded
            ATHMdf = ATHMdf.append(row)
            ATHMindexes.append(index)

    data = data.drop(index=ATHMindexes)
    return data, ATHMdf