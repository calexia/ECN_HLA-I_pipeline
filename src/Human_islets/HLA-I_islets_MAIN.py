import os

import sys
# Set path as the root project folder to allow relative imports
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

import pandas as pd
import numpy as np

import src.toolbox.initialization as init
import src.toolbox.reformat_entering as reformat
import src.toolbox.length_filter as length
import src.toolbox.pepOcc_geneAnnot_single_file as occAnnot
import src.toolbox.ubiquity_filter as ubi
import src.toolbox.enrichment_filter as enrich
import src.toolbox.concatenation_singleFile as concat
import src.toolbox.netMHC_I_pipeline as MHC
import src.toolbox.AAoccurence as AA

#islet 22:
LIST_ISLETS_1 = ["HLA-A02:01", "HLA-A03:01", "HLA-B07:02", "HLA-B13:01", "HLA-C06:02", "HLA-C15:02", "HLA-E01:01"]
#islet 27:
LIST_ISLETS_2 = ["HLA-A02:01", "HLA-A03:01", "HLA-B50:01", "HLA-B51:01", "HLA-C02:02", "HLA-C06:02", "HLA-E01:01"]

command_arguments = sys.argv[1:]
if (len(command_arguments) == 3):
    ISLET_FILE_PATH = command_arguments[0]
    ALLELES_LIST_SELECTION = command_arguments[1]
    if (int(ALLELES_LIST_SELECTION) == 1):
        ALLELES_LIST = LIST_ISLETS_1
    elif (int(ALLELES_LIST_SELECTION) == 2):
        ALLELES_LIST = LIST_ISLETS_2
    else:
        print('Wrong islet_list_selection parameter value, set 1 or 2')
        exit(0)
    OUTPUT_FOLDER_PATH = command_arguments[2]
else:
    print('Missing arguments : python3 src/Human_islets/HLA-I_islets_MAIN_NatureComm.py <islet_file_path> <islet_list_selection> <output_folder_path>')
    exit(0)

DATA_TYPES = [("islet1", ISLET_FILE_PATH)]

SANDBERG_FILE_PATH = './resources/Sandberg_cellType_enrichment.txt'
ANNOTATION_FILE_PATH = './resources/mainAnnot_homo_sapiens.txt'
ATLAS_FILE_PATH = './resources/proteinatlas.tsv'
HPRD_FILE_PATH = './resources/TISSUE_EXPRESSIONS.txt'
NORMAL_TISSUE_FILE_PATH = './resources/normal_tissue.tsv'

"""STEP 0: manually remove and isolate peptides explained only by a mRNA variant"""

"""INITIALIZATION : creation of the needed dictionaries"""
print("Generating the needed dictionaries...")
# Reading of files for dictionaries generation
annotationFile = pd.read_csv(ANNOTATION_FILE_PATH,
                             sep='\t',
                             dtype = {
                                 "Gene name":"str",
                                 "UniProt accession" : "str"
                             }
                             )
atlasFile = pd.read_csv(ATLAS_FILE_PATH,
                        sep='\t',
                        dtype = {
                            "Gene" : 'str',
                            "Gene description" : "str",
                            "Gene synonym" : "str",
                            "RNA TS TPM" : "str",
                            "RNA tissue category" : "str"
                        }
                        )
HPRDFile = pd.read_csv(HPRD_FILE_PATH,
                       sep='\t',
                       dtype = {
                           "GeneSymbol" : 'str',
                           "Expression_term" : "str"
                       }
                       )

normalTissueFile = pd.read_csv(NORMAL_TISSUE_FILE_PATH,
                               sep = '\t',
                               dtype={
                                   "Gene name" : 'str',
                                   "Tissue" : 'str',
                                   'Level' : 'str'
                               }
                               )

def processDict(annotationFile, atlasFile, HPRDFile, normalTissueFile):
#Import files as dfs
    annotations = pd.DataFrame(annotationFile)
    atlas = pd.DataFrame(atlasFile)
    HPRDdf = pd.DataFrame(HPRDFile)
    normalTissuedf = pd.DataFrame(normalTissueFile)

#Apply processing functions = step 1 of the pipeline
    annotationDict, RNATSTPMDict, geneSynonymDict, geneDescriptionDict, expTermDict, normalTissueDict = init.createDictionaries(annotations, atlas, HPRDdf, normalTissuedf)
    expTermDict = init.reformatDictValue(expTermDict)

    return annotationDict, RNATSTPMDict, geneSynonymDict, geneDescriptionDict, expTermDict, normalTissueDict


"""REFORMATTING ENTERING FILE"""
def reformat_entering_data(data, type):
    data = reformat.reformatAccession(data, type)
    data = data[data.Accession != 0]
    data = reformat.reduceTable(data)

    return data

"""FILTERING FOR BETA CELL ENRICHMENT"""
def betaCell_enrichment(data, ubidf, ATHMdf, sandbergdf, geneDescriptionDict, expTermDict, normalTissueDict, type):
    data = enrich.expTermExpression(data, expTermDict)
    data = enrich.addNormalTissue(data, normalTissueDict)
    data = enrich.reformatRNATSTPM_expTerm(data) #reformats only RNATSTPM col
    data = enrich.litEvidence(data)
    data = enrich.RNAevidence(data)
    data = enrich.protEvidence(data)
    data, ubidf = enrich.RNAprotLitEvidenceFilter(data, ubidf, type)
    data, ubidf = enrich.filterSandberg(data, ubidf, sandbergdf, type)
    data = enrich.addGeneDescri(data, geneDescriptionDict)
    data, ATHMdf = enrich.removeATHM(data, ATHMdf)

    return data, ubidf, ATHMdf

"""CONCATENATION OF RESULTS TABLES"""
def clean(allFormattedData):
    finalTable = allFormattedData[0]

    finalTable = concat.cleanResults(finalTable)

    return finalTable

"""GENERTING MHC PREDICTIONS"""

def MHCPredictionGeneration(data, filePath, type, allelesList):
    MHC.computeBindingPredictions(data, filePath, CURRENT_TYPE, allelesList, OUTPUT_FOLDER_PATH)
    restrictiondfList = []  # initializes empty list

    first3colToName = ["Pos", "Peptide", "ID"]
    first2colToDrop = ["Pos", "ID"]
    last2col = ["Ave", "NB"]
    HLAcolNameList = ["Rank_" + s for s in allelesList]
    repeatedCol = ["core", "icore", "EL-score"]

    preColNames = []
    preColtoDrop = []

    for i in range(0, len(HLAcolNameList)):
        preColNames.extend(repeatedCol)
        preColNames.append(HLAcolNameList[i])

        preColtoDrop.extend(repeatedCol)

    first3colToName.extend(preColNames)
    colNames = first3colToName
    colNames.extend(last2col)

    first2colToDrop.extend(preColtoDrop)
    colsToDrop = first2colToDrop
    colsToDrop.append(last2col[0])

    for j in range(7):
        predictionFiles = pd.read_csv(OUTPUT_FOLDER_PATH + '{}mersPredictions{}.txt'.format(8 + j, type), sep='\t', header=[1])

        predictiondf = pd.DataFrame(predictionFiles)
        predictiondf.columns = colNames
        predictiondf = predictiondf.drop(columns=colsToDrop)
        for col, content in predictiondf.iteritems():
            if "Rank_" in str(col):
                predictiondf = predictiondf.astype(
                    { col: np.float64 }
                )
        restrictiondf = MHC.bindingPredictionsTreatment(predictiondf, allelesList)
        restrictiondfList.append(restrictiondf)

    concatRestrictiondf = pd.concat(restrictiondfList)

    dataAssigned = MHC.finalTableHLA_Assignment(data, concatRestrictiondf, OUTPUT_FOLDER_PATH)
    return dataAssigned


#if __name__ == "main":
"""PIPELINE : FILTERING STEPS"""
#Call of the function to generate the dictionaries
annotationDict, RNATSTPMDict, geneSynonymDict, geneDescriptionDict, expTermDict, normalTissueDict = processDict(annotationFile, atlasFile, HPRDFile, normalTissueFile)

"""Consideration of datatypes """
#Takes into account both datatypes
allFormattedData = [] #initializes empty list to keep all results table
for type in DATA_TYPES:
    global CURRENT_TYPE  # defines a global variable for type
    CURRENT_TYPE = type[0]
    CURRENT_TYPE_FILE = type[1]
    dataFile = pd.read_csv(CURRENT_TYPE_FILE.format(CURRENT_TYPE),
                           sep=',',
                           dtype = {
                               "Peptide" : "str",
                               "Length" : np.int64,
                               "Fraction" : np.int64,
                               "Accession" : "str",
                               "PTM" : "str"
                           }
                           )
    sandbergFile = pd.read_csv(SANDBERG_FILE_PATH,
                               sep='\t',
                               dtype = {
                                   "Genes": 'str',
                                   "sampleID": 'str',
                                   "gamma cell": np.float64,
                                   "alpha cell": np.float64,
                                   "beta cell": np.float64,
                                   "acinar cell": np.float64,
                                   "epsilon cell": np.float64,
                                   "ductal cell": np.float64,
                                   "delta cell": np.float64,
                                   "Cell enrichment" : "str"
                               }
                               )

    print("\nRUNNING THE ANALYSIS OF THE {} FILE".format(CURRENT_TYPE))

    # Import files as dfs
    data = pd.DataFrame(dataFile)
    sandbergdf = pd.DataFrame(sandbergFile)

    data = data.fillna(0)  # pandas cannot deal with empty cells, fill all with 0

    #FUNCTIONS CALLING
        # call of function to reformat data to enter pipeline
    data = reformat_entering_data(data, CURRENT_TYPE)
    data.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + 'DB search psm_{}_Entering.txt'.format(CURRENT_TYPE),
                         index=False, sep='\t')
    # counts the AA occurences of all peptides and saves it into a specific file
    countAll = AA.countAAoccurence(data)
    countPercentageAll = AA.computeAAocc(countAll)
    countPercentageAll.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + '{}_AAcounts_all.txt'.format(CURRENT_TYPE),
                              index=False, sep='\t')
    # call of function to filter on length of peptide
    data = length.lengthFilter(data, 8, 14, CURRENT_TYPE)
    # counts the AA occurences of all peptides and saves it into a specific file
    countLengthFiltered = AA.countAAoccurence(data)
    countPercLengthFiltered = AA.computeAAocc(countLengthFiltered)
    countPercLengthFiltered.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + '{}_AAcounts_LengthFiltered.txt'.format(CURRENT_TYPE),
                              index=False, sep='\t')
    #calls the function to add "raw peptide" column
    data = concat.rawPeptide(data)
    #call of function to generate MHC binding predictions
    data = MHCPredictionGeneration(data, OUTPUT_FOLDER_PATH, CURRENT_TYPE, ALLELES_LIST)
    #call of function to add peptide occurence AND gene annotation
    data = occAnnot.countOccurences(data)
    data = occAnnot.annotGene(data, annotationDict)
    data.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + '{}_Length_Filtered.txt'.format(CURRENT_TYPE),
                        index=False, sep='\t')
    #call of function to filter ubi from non ubi peptides

    data = ubi.addGeneSyno(data, geneSynonymDict)
    data = ubi.nonUbiExpression(data, RNATSTPMDict)

    #call of function to filter peptides enriched in beta cells
    ubidf = pd.DataFrame()
    ATHMdf = pd.DataFrame()
    formattedData, ubidf, ATHMdf = betaCell_enrichment(data, ubidf, ATHMdf, sandbergdf, geneDescriptionDict, expTermDict, normalTissueDict, type)

    allFormattedData.append(formattedData)

    #saves intermediate datframes as csv
    formattedData.to_csv(path_or_buf= OUTPUT_FOLDER_PATH + 'DB search psm_{}_Formatted.txt'.format(CURRENT_TYPE), index=False, sep='\t')
    ubidf.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + '{}_ubiquitous_peptides.txt'.format(CURRENT_TYPE), index=False, sep='\t')
    ATHMdf.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + '{}_ATHM_peptides.txt'.format(CURRENT_TYPE), index=False, sep='\t')

    # dataFilteredOut.to_csv(path_or_buf=OUTPUT_FILE_PATH + '{}curated_out_peptides.txt'.format(CURRENT_TYPE), index=False, sep='\t')


"""PIPELINE : CONCATENATION OF RESULTS"""
#calls the function to concat and clean the df
finalTable = clean(allFormattedData)
# counts the AA occurences of all peptides and saves it into a specific file
countFinal = AA.countAAoccurence(finalTable)
countPercFinal = AA.computeAAocc(countFinal)
countPercFinal.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + 'AAcounts_betaEnriched.txt',
                                   index=False, sep='\t')
#saves if to a csv
finalTable.to_csv(path_or_buf=OUTPUT_FOLDER_PATH + 'CleanResultTable.txt', index=False, sep='\t')

