import re
import pandas as pd

""" STEP 2b : keep all the accessions ID """
def reformatAccession (data, type):
    uniqueValues = data["Peptide"].nunique() #count the number of unique peptides
    print ("The {} table contains {} unique peptides, before any filtering".format(type, uniqueValues))

    print('REFORMAT : Reformatting accession...')
    poolIDs = re.compile('(\w*)\|') #compile the word of any length written in front of the |

    for index, row in data.iterrows():
        accessionID = row['Accession'] #localize the accession column
        if accessionID != 0:
            accessionIDs = str(accessionID) #make it a string to apply findall
            # if "mRNA" in accessionIDs:
                # print("mRNA contained in accession {} at row {}".format(accessionIDs, index))
            accessionIDs = poolIDs.findall(accessionIDs)
            newAccessionIDs = ';'.join(accessionIDs)
            data.at[index, 'Accession'] = newAccessionIDs
        # else:
        #     print("Indexed {} peptide {} has no accession".format(index, row['Peptide']))

    return data

""" STEP 2c: reduce the table to have a single line per peptide"""
def reduceTable(data):
    print("REFORMAT : Reducing table...")
    dataReduced = pd.DataFrame(
        columns=[
            "Peptide", "Length", "Fraction", "Accession", "PTM"
        ]
    ) #creates a new dataframe with set col labels
    dataReduced.set_index("Peptide") #set the col peptide as the df index

    for index, row in data.iterrows(): #iterates the df over rows
        peptide = row["Peptide"] #get the peptide of each row
        accession = row["Accession"] #get the accession of each row
        fraction = row["Fraction"] #get the fraction of each row
        if peptide not in dataReduced.index: #check peptide not already in the index
            dataReduced.loc[peptide] = [
                peptide,
                row["Length"],
                str(fraction),
                str(accession),
                row["PTM"],
            ] #fills up all the cells for the considered peptide, based on current row
        else:
            alreadyExistingRow = dataReduced.loc[peptide] #if peptide already in the df, creates a variable for the corresponding row

            #concat Accessions, gene names and gene descriptions (as accession goes with both gene name and gene description)
            alreadyExistingAccession =  alreadyExistingRow["Accession"] #get the accession of exisiting row
            newAlreadyExistingAccession = alreadyExistingAccession

            accession_str = (str(accession)).split(';') #split current accession in str
            alreadyExistingAccession_split = alreadyExistingAccession.split(';') #split exisiting accession in str
            for i in accession_str: #for ele in current accession
                if i not in alreadyExistingAccession_split:
            # if accession not in alreadyExistingAccession:
                    newAlreadyExistingAccession = newAlreadyExistingAccession + ';' + i

            dataReduced.at[peptide, "Accession"]= newAlreadyExistingAccession #add the current accession, that is not already in existing, in existing accession

            #concat Fractions
            alreadyExistingFraction = alreadyExistingRow["Fraction"] #defines the existing fraction
            if str(fraction) not in alreadyExistingFraction.split(";"):	#if current row fraction not in existing fraction
                dataReduced.at[peptide, "Fraction"] = (alreadyExistingFraction + ";" + str(fraction)) #append current fraction

    return dataReduced
