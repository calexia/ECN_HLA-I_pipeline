import pandas as pd

FILE_PATH = "D://DataAlexia/Projets/Scripts/HLA-I/IFNa_dataset/"
INPUT_FILE_PATH = FILE_PATH + "Input/"


atlasFile = pd.read_csv(INPUT_FILE_PATH + 'proteinatlas.tsv',
                        sep='\t',
                        dtype = {
                            "Gene" : 'str',
                            "Gene description" : "str",
                            "Gene synonym" : "str",
                            "RNA TS TPM" : "str",
                            "RNA tissue category" : "str"
                        }
                        )
HPRDFile = pd.read_csv(INPUT_FILE_PATH + 'TISSUE_EXPRESSIONS.txt',
                       sep='\t',
                       dtype = {
                           "GeneSymbol" : 'str',
                           "Expression_term" : "str"
                       }
                       )

normalTissueFile = pd.read_csv(INPUT_FILE_PATH + "normal_tissue.tsv",
                               sep = '\t',
                               dtype={
                                   "Gene name" : 'str',
                                   "Tissue" : 'str',
                                   'Level' : 'str'
                               }
                               )


RNAevidenceOfInterest = atlasFile["RNA TS TPM"] == "pancreas"

LitEvidenceOfPancreas = HPRDFile["Expression_term"] == "Pancreas"
LitEvidenceOfIslets = HPRDFile["Expression_term"] == "Islets of Langerhans"
LitEvidenceOfBetaCell = HPRDFile["Expression_term"] == "Beta cell"
LitEvidenceOfInterest = HPRDFile[LitEvidenceOfBetaCell | LitEvidenceOfIslets | LitEvidenceOfPancreas]

NormalTissueLevelHigh = normalTissueFile["Level"] == "High"
NormalTissueLevelMed = normalTissueFile["Level"] == "Medium"
NormalTissueReliability = normalTissueFile["Reliability"] != "Uncertain"
NormalTissuePancreas = normalTissueFile["Tissue"] == "pancreas"
NormalTissueEvidenceOfInterest = normalTissueFile[(NormalTissueLevelHigh | NormalTissueLevelMed) & NormalTissueReliability & NormalTissuePancreas]


LitEvidenceOfInterest.to_csv(FILE_PATH + "Output/LitEvidence.txt", sep = "\t")
RNAevidenceOfInterest.to_csv(FILE_PATH + "Output/RNAEvidence.txt", sep = "\t")
NormalTissueEvidenceOfInterest.to_csv(FILE_PATH + "Output/ProtEvidence.txt", sep = "\t")