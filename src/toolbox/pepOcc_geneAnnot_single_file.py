"""STEP 3 : count occurences"""
def countOccurences(data):
    print("OCCANNOT : Counting the occurences of each peptides...")
    data["Total occurences"] = 0 #initialize value of occurences
    data["Sample occurences"] = 0 #initialize value of occurences


    for index, row in data.iterrows(): #iterates over row in df
        fractions = row["Fraction"] #get fraction for each row
        fractions_split = fractions.split(';') #put fraction into a list
        totalOccurences = len(fractions_split) #count the number of total occurences
        data.at[index, "Total occurences"] = totalOccurences #put the value in the corresponding cell in the df
        data.at[index, "Sample occurences"] = totalOccurences #put the value in the sample occurence cell in the df

    return data



def countOccurences_Sergio(data, type):
    print("OCCANNOT : Counting the occurences of each peptides...")
    data["Total occurences"] = 0 #initialize value of occurences
    data["Sample occurences"] = 0 #initialize value of occurences

    if "B" in str(type): #check data type
        pairedFractionsList = [
            ['1', '2'],
            ['3', '4'],
            ['5', '6'],
            ['7', '8']
        ] #defines the fraction numbers corresponding to the data type
    elif "IFNa" in str(type): #check data type
        pairedFractionsList = [
        ['9' , '10'],
        ['11' , '12'],
        ['13' , '14'],
        ['15' , '16']
    ] #defines the fraction numbers corresponding to the data type
    else:
        print('ERROR : no type found')
        quit() #anticipates potential error

    for index, row in data.iterrows(): #iterates over row in df
        fractions = row["Fraction"] #get fraction for each row
        fractions_split = fractions.split(';') #put fraction into a list
        totalOccurences = len(fractions_split) #count the number of total occurences
        data.at[index, "Total occurences"] = totalOccurences #put the value in the corresponding cell in the df
        data.at[index, "Sample occurences"] = totalOccurences #put the value in the sample occurence cell in the df

        for pairedFraction in pairedFractionsList: #for ele in list of fraction pairs
            if all(elem in fractions_split for elem in pairedFraction): #check if ALL the elements of the fraction list are in the paire fraction list
               data.at[index, "Sample occurences"] -= 1 #recount the sample occurence and put it back in the corresponding cell

    return data

""" STEP 4 : use the gene annotation dict to match the gene name to accession"""
def annotGene(data, annotationDict):
    print('OCCANNOT : Adding gene annotation...')
    for index, row in data.iterrows(): #iterates over row
        accessionList = row["Accession"].split(';') #get all the accession for each row as a list
        genes = [] #initializes empty list
        geneList = [] #initializes empty list
        for accession in accessionList: #call each ele in accession list
            sep = ';'
            if accession in annotationDict: #check if accession is in previously defined dict (as a key)
                genes.append(str(annotationDict[accession])) #add str of dict key in list // added str function otherwise
            geneList=sep.join(genes) #join the list in a bigger list
        data.at[index, "Gene name"] = geneList #put the bigger list back in df
    return data