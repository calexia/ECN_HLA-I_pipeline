
"""STEP 2a: remove several columns, keeps only peptide of specific length"""
def lengthFilter (data, min, max, type):
    print('FILTERING : filtering length...')
    # uniqueValues = data["Peptide"].nunique() #count the number of unique peptides
    # print ("The {} table contains {} unique peptides, before any filtering".format(type, uniqueValues))

    #colToDrop = ['-10lgP', 'Mass', 'ppm', 'm/z', 'Z', 'RT', 'Area', 'Scan', 'Id', 'from Chimera','Source File','AScore', 'Found By'] #col to remove
    #data=data.drop(columns=colToDrop) #drop cols
    data=data[(data["Length"] >= min) & (data["Length"] <= max)] #peptide length to keep

    uniqueValuesLengthRestrict = data["Peptide"].nunique() #count number of unique peptides, after length restriction
    print ("The {} table contains {} unique peptides, of length between {} and {} amino acids".format(type, uniqueValuesLengthRestrict, min, max))

    return data
