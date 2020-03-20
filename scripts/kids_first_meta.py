#This function takes 4 arguments: manifest file, biospecimen file, clinical file, and the name/location of output file.

# Import necessary packages
import time
from sys import argv
import urllib.request as url
import os
import re

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Input/output files
manifestFile = argv[1]
biospecimenFile = argv[2]
clinicalFile = argv[3]
outputFile = argv[4]

# Dictionary to store needed information in
outputDict = {}

# Obtain information from manifest file and add to initial dictionary
with open(manifestFile) as manifest:
    manifestColumnNames = manifest.readline()
    manifestColumnNames = manifestColumnNames.rstrip().split("\t")
    # Information needed from this file are the "File Name", "Family Id", "Sample External ID", and "Proband" (Yes, or No).
    fileNameIndex = manifestColumnNames.index("File Name")
    familyIdIndex = manifestColumnNames.index("Family Id")
    manifestSampleIdIndex = manifestColumnNames.index("Sample External ID")
    probandIndex = manifestColumnNames.index("Proband")
    
    # Add to the initial dictionary where the key is the "Aliquot External ID" as this is common across all input files.
    # The value is a list where value[0] the "File Name", value[1] is "Family Id", and value[2] is "Proband"
    # as "Yes" or "No"
    for sample in manifest:
        sample = sample.rstrip().split("\t")
        if sample[probandIndex] == "Yes" or sample[probandIndex] == "No":
            externalID = sample[manifestSampleIdIndex]
            if "Schiffman" not in externalID:
                outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[familyIdIndex], sample[probandIndex]]
            else:
                externalID = re.findall(r"Schiffman\-\w+", externalID)[0]
                outputDict[externalID] = [sample[fileNameIndex], sample[familyIdIndex], sample[probandIndex]]
        # Sometimes "Yes" or "No" is not listed under "Proband". Therefore, this else statement will check NCBI for affected
        # status based on "Aliquot External ID"
        else:
            getAffectedStatus = str(url.urlopen(f"https://www.ncbi.nlm.nih.gov/biosample/?term={sample[manifestSampleIdIndex]}").read())
            if "subject is affected</th><td>Yes" in getAffectedStatus:
                outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[familyIdIndex], "Yes"]
            elif "subject is affected</th><td>No" in getAffectedStatus:
                outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[familyIdIndex], "No"]

# Obtain information from the biospecimen file and add to initial dictionary
with open(biospecimenFile) as biospecimen:
    biospecimenColumnNames = biospecimen.readline()
    biospecimenColumnNames = biospecimenColumnNames.rstrip().split("\t")
    # Information needed from this file are the "External Aliquot Id", and "Biospecimens Id"
    biospecimenExternalIdIndex = biospecimenColumnNames.index("External Sample Id")
    sampleIdIndex = biospecimenColumnNames.index("Biospecimens Id")

    # Use the patient "External Aliquot Id" as the key to append the "Biospecimens Id" to the value list.
    for line in biospecimen:
        line = line.rstrip().split("\t")
        externalID = line[biospecimenExternalIdIndex]
        if "Schiffman" not in externalID:
            outputDict[line[biospecimenExternalIdIndex]].append(line[sampleIdIndex])
        else:
            try:
                externalID = re.findall(r"Schiffman\-\w+", externalID)[0]
                outputDict[externalID].append(line[sampleIdIndex])
            except:
                print(f"WARNING: {externalID} not found in biospecimen file!")
                continue

# Obtain gender information from the clinical file and add to initial dictionary
with open(clinicalFile) as clinical:
    clinicalColumnNames = clinical.readline()
    clinicalColumnNames = clinicalColumnNames.rstrip().split("\t")
    #Information needed from this file are the "External Id", and "Gender"
    clinicalExternalIdIndex = clinicalColumnNames.index("External Id")
    genderIndex = clinicalColumnNames.index("Gender")

    # Use the patient "External Id" as the key to append the "Biospecimens Id" to the value list.
    for line in clinical:
        line = line.rstrip().split("\t")
        externalID = line[clinicalExternalIdIndex]
        if "Schiffman" not in externalID and externalID in outputDict:
            outputDict[line[clinicalExternalIdIndex]].append(line[genderIndex])
        else:
            try:
                externalID = re.findall(r"Schiffman\-\w+", externalID)[0]
                outputDict[line[clinicalExternalIdIndex]].append(line[genderIndex])
            except:
                print(f"WARNING: {externalID} not found in clinical file!")
                continue

# Create list of samples where only trios are included
familyDict = {}
for key, value in outputDict.items():
    if value[1] not in familyDict:
        familyDict[value[1]] = [value]
    else:
        familyDict[value[1]].append(value)

trioList = []
for key, value in familyDict.items():
    if len(value) == 3:
        trioList.append(value[0])
        trioList.append(value[1])
        trioList.append(value[2])

# Output information to outputFile
with open(outputFile, 'w') as output:
    output.write("file_name\tfamily_id\tsample_id\tproband\tsex\n")
    for item in trioList:
        if item[-1] == "Female":
            output.write(f"{item[0]}\t{item[1]}\t{item[3]}\t{item[2]}\t2\n")
        else:
            output.write(f"{item[0]}\t{item[1]}\t{item[3]}\t{item[2]}\t1\n")

# Output message and time complete
print(f"{char}There are {len(trioList) / 3} trios.{char}")
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')