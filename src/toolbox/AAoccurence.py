import pandas as pd
from collections import Counter
#
# dataFile = pd.read_csv("/home/alexia/Documents/Codes/HLA-I/Output/DB search psm_Basal_Entering.txt", sep="\t")
# data=pd.DataFrame(dataFile)

def countAAoccurence(data):

    print("AA : counting the occurences of each aa")
    AA= [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V"
    ]

    data = data.drop_duplicates(subset = "Peptide", keep= 'first') #drop duplicates peptides so that AA count is not biased by most abundant peptides
    data = data.reset_index()
    aaOccCount = data[["Peptide"]].join(pd.DataFrame([Counter(peptide) for peptide in data["Peptide"]]).reindex(list(AA), axis=1).fillna(0).astype(int))
    # aaOccCount = data[["Peptide"]]\
    #     .join( #considers the col peptide from the data AND will join
    #     pd.DataFrame( #a dataframe
    #         [Counter(peptide) for peptide in data["Peptide"]]) #which contains the count of peptide in the data
    #         .reindex(list(AA), axis=1) #indexed with the AA list, as a col
    #         .fillna(0) #fills empty value with 0
    #         .astype(int)) #make sure the counts are all int
    return aaOccCount

def computeAAocc(aaOccCount):
    print("AA : Computing percentages")

    aaOcc = aaOccCount.sum(axis=0) #count the total number of AA
    aaOcc = aaOcc[1:]
    aaOccdf = pd.DataFrame(data=aaOcc, columns = ["Total counts"])

    somme = float(aaOccdf.sum(axis=0))

    for index, row in aaOccdf.iterrows():
        totCount = float(row["Total counts"])
        percentage = float((totCount / somme) * 100)
        aaOccdf.at[index, "Percentage (%)"] = percentage



    return aaOccdf


# countAll = countAAoccurence(data)
# countPercentageAll = computeAAocc(countAll)