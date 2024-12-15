import re

"""STEP 1a : create the dictionaries"""
def createAnnotDict (annotationdf):
    print('INIT : Generating gene annotation dictionary...')
    annotationDict = {} #initialize empty dict
    for index, row in annotationdf.iterrows(): #iterates over the rows of annotation df
        gene = row["Gene name"] #get gene name for each row
        accessionList = row["UniProt accession"].split(';') #get uniprot accessions list for each row
        for accession in accessionList: #get each accession from accession list
            annotationDict[accession] = gene #attributes each accession to the gene name, in a dict

    return annotationDict

def createAtlasBasedDict(atlas):
    print('INIT : Generating RNA Tissue category AND RNA TS TPM dictionaries...')
    atlas = atlas.fillna("")  # fills empty cells of the df // format issues
    #RNATissueCatDict = {}  # intialize empty dict
    RNATSTPMDict = {}  # intialize empty dict
    geneSynonymDict = {}  # intialize empty dict
    geneDescriptionDict = {}  # intialize empty dict
    for index, row in atlas.iterrows():  # iterates cover each row in atlas df
        gene = row["Gene"]  # get gene name for each row
        #RNATissueCatList = row["RNA tissue category"]  # get RNA tissue cat for each row
        RNATSTPMList = row["RNA TS TPM"]  # get RNA TS TPM for each row
        geneSynonymList = row["Gene synonym"]  # get Gene for each row
        geneDescriptionList = row["Gene description"]  # get gene descri for each row

        #RNATissueCatDict[gene] = RNATissueCatList  # attributes RNA tissue cat to each gene, in a dict
        RNATSTPMDict[gene] = RNATSTPMList  # attributes RNA TS TPM value to each gene, in a dict
        geneSynonymDict[gene] = geneSynonymList  # attributes gene synonym to each gene, in a dict
        geneDescriptionDict[gene] = geneDescriptionList  # attributes gene description to each gene, in a dict

    return RNATSTPMDict, geneSynonymDict, geneDescriptionDict

def createNormalTissueDict(normalTissuedf):
    print("INIT : Generating normal tissue dictionary...")

    normalTissueDict = {}
    for index, row in normalTissuedf.iterrows():
        gene = row["Gene name"]
        level = row["Level"]
        tissue = row["Tissue"]
        reliability = row["Reliability"]

        if reliability == "Uncertain":
            continue
        else:

            if level == "Medium"\
                or level == "High":

                if gene not in normalTissueDict:
                    normalTissueDict[gene] = [tissue]
                else:
                    normalTissueDict[gene].append(tissue)

    return normalTissueDict

def createExpTermDict(HPRDdf):

    print('INIT : Generating expression_term tissue category dictionary...')
    expTermDict = {} #initialize empty dict
    for index, row in HPRDdf.iterrows(): #iterates over rows of HPRDFile
        geneSymbol = row["GeneSymbol"] #get gene of each row
        expTerm = row["Expression_term"] #get expTerm of each row

        if geneSymbol not in expTermDict: #check if key not already in dict
            expTermDict[geneSymbol] = [expTerm] #if not, add key and value
        else:
            expTermDict[geneSymbol].append(expTerm) #if already in dict, add value to existing key

    return expTermDict

def createDictionaries(annotationdf, atlas, HPRDdf, normalTissuedf):
#creates 3 dictionaries
    annotationDict = createAnnotDict(annotationdf)
    RNATSTPMDict, geneSynonymDict, geneDescriptionDict = createAtlasBasedDict(atlas)
    expTermDict = createExpTermDict(HPRDdf)
    normalTissueDict = createNormalTissueDict(normalTissuedf)

    return annotationDict, RNATSTPMDict, geneSynonymDict, geneDescriptionDict, expTermDict, normalTissueDict


"""STEP 1b : Reformat value of each key in expTermDict"""
def reformatDictValue(expTermDict):
    print ("INIT : Reformating expTerm dictionary...")

    poolValues = re.compile('[a-zA-Z\s]{2,}') #regex for any word, space included, 1 character min

    for key, value in expTermDict.items(): #iterates over each key
        value = poolValues.findall(str(value)) #find all the values corresponding to the regex
        newValue = ';'.join(value)#put all the values in a str, joined by ;
        expTermDict[key] = newValue #change the value for each key by the new value

    return expTermDict
