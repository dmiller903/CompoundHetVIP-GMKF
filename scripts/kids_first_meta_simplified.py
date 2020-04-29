#This function takes 4 arguments: tsv file, biospecimen file, clinical file, and the name/location of output file.

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
tsvTable = argv[1]
clinicalFile = argv[2]
outputFile = argv[3]

# Dictionary to store needed information in
genderDict = {}
outputDict = {}

# Obtain gender information from the clinical file and add to initial dictionary
with open(clinicalFile) as clinical:
    clinicalColumnNames = clinical.readline()
    clinicalColumnNames = clinicalColumnNames.rstrip().split("\t")
    #Information needed from this file are the "Participant ID" or "Kf Id", and "Gender"
    try:
        clinicalParticipantIdIndex = clinicalColumnNames.index("Participant ID")
    except:
        clinicalParticipantIdIndex = clinicalColumnNames.index("Kf Id")
    genderIndex = clinicalColumnNames.index("Gender")
    for line in clinical:
        lineList = line.rstrip().split("\t")
        genderDict[lineList[clinicalParticipantIdIndex]] = lineList[genderIndex]

# Obtain information from tsv file and add to initial dictionary
with open(tsvTable) as table:
    tableColumnNames = table.readline()
    tableColumnNames = tableColumnNames.strip().split("\t")
    # Information needed from this file are the "File Name", "Family Id", "Participants ID", and "Proband" (Yes, or No)
    # and "Biospecimen ID".
    fileNameIndex = tableColumnNames.index("File Name")
    familyIdIndex = tableColumnNames.index("Family Id")
    tableSampleIdIndex = tableColumnNames.index("Participants ID")
    probandIndex = tableColumnNames.index("Proband")
    biospecimenIndex = tableColumnNames.index("Biospecimen ID")
    # Add to the initial dictionary where the key is the "Participants ID" as this is common across all input files.
    # The value is a list where value[0] the "File Name", value[1] is "Participants ID, value[2] is "Family Id", 
    # and value[3] is "Proband" as "Yes" or "No"
    for sample in table:
        sample = sample.rstrip().split("\t")
        gender = genderDict[sample[tableSampleIdIndex]]
        outputDict[sample[tableSampleIdIndex]] = [sample[fileNameIndex], sample[tableSampleIdIndex], sample[familyIdIndex], sample[probandIndex], sample[biospecimenIndex], gender]


# Create list of samples where only trios are included
familyDict = {}
for key, value in outputDict.items():
    if value[2] not in familyDict:
        familyDict[value[2]] = [value]
    else:
        familyDict[value[2]].append(value)

trioList = []
for key, value in familyDict.items():
    if len(value) == 3:
        trioList.append(value[0])
        trioList.append(value[1])
        trioList.append(value[2])

# Output information to outputFile
with open(outputFile, 'w') as output:
    output.write("file_name\tfamily_id\tsample_id\tproband\tsex\n")
    total_excluded = 0
    for item in trioList:
        if item[-1] == "Female":
            output.write(f"{item[0]}\t{item[2]}\t{item[4]}\t{item[3]}\t2\n")
        else:
            output.write(f"{item[0]}\t{item[2]}\t{item[4]}\t{item[3]}\t1\n")

# Output message and time complete
print(f"{char}There are {(len(trioList) / 3) - total_excluded} trios.{char}")
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')