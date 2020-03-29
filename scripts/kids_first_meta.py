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
clinicalFile = argv[2]
biospecimenFile = argv[3]
outputFile = argv[4]

# Dictionary to store needed information in
outputDict = {}

# Obtain information from manifest file and add to initial dictionary
with open(manifestFile) as manifest:
    manifestColumnNames = manifest.readline()
    manifestColumnNames = manifestColumnNames.rstrip().split("\t")
    # Information needed from this file are the "File Name", "Family Id", "Participants ID", and "Proband" (Yes, or No)
    # and "Sample External ID".
    fileNameIndex = manifestColumnNames.index("File Name")
    familyIdIndex = manifestColumnNames.index("Family Id")
    manifestSampleIdIndex = manifestColumnNames.index("Participants ID")
    probandIndex = manifestColumnNames.index("Proband")
    externalIdIndex = manifestColumnNames.index("Sample External ID")
    # Add to the initial dictionary where the key is the "Participants ID" as this is common across all input files.
    # The value is a list where value[0] the "File Name", value[1] is "Participants ID, value[2] is "Family Id", 
    # and value[3] is "Proband" as "Yes" or "No"
    participantID_externalID = {}
    for sample in manifest:
        sample = sample.rstrip().split("\t")
        participantID_externalID[sample[manifestSampleIdIndex]] = sample[externalIdIndex]
        if sample[probandIndex] == "Yes" or sample[probandIndex] == "No":
            outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[manifestSampleIdIndex], sample[familyIdIndex], sample[probandIndex]]
        # Sometimes "Yes" or "No" is not listed under "Proband". Therefore, this else statement will check NCBI for affected
        # status based on "Sample External ID"
        else:
            externalID = sample[externalIdIndex]
            getAffectedStatus = str(url.urlopen(f"https://www.ncbi.nlm.nih.gov/biosample/?term={externalID}").read())
            if "subject is affected</th><td>Yes" in getAffectedStatus:
                outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[manifestSampleIdIndex], sample[familyIdIndex], "Yes"]
            elif "subject is affected</th><td>No" in getAffectedStatus:
                outputDict[sample[manifestSampleIdIndex]] = [sample[fileNameIndex], sample[manifestSampleIdIndex], sample[familyIdIndex], "No"]

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

    # Use the Participant ID as the key to append the "Gender" to the value list.
    for line in clinical:
        line = line.rstrip().split("\t")
        outputDict[line[clinicalParticipantIdIndex]].append(line[genderIndex])

# Obtain biospecimen ID
biospecimenDict = {}
with open(biospecimenFile) as biospecimen:
    headerList = biospecimen.readline().rstrip("\n").split("\t")
    print(headerList)
    biospecimenIDIndex = headerList.index("Biospecimens Id")
    externalSampleIDIndex = headerList.index("External Sample Id")
    for line in biospecimen:
        lineList = line.rstrip("\n").split("\t")
        if lineList[externalSampleIDIndex] in participantID_externalID.values():
            biospecimenDict[lineList[externalSampleIDIndex]] = lineList[biospecimenIDIndex]
print(biospecimenDict)


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
    for item in trioList:
        externalID = participantID_externalID[item[1]]
        biospecimenID = biospecimenDict[externalID]
        if item[-1] == "Female":
            output.write(f"{item[0]}\t{item[2]}\t{biospecimenID}\t{item[3]}\t2\n")
        else:
            output.write(f"{item[0]}\t{item[2]}\t{biospecimenID}\t{item[3]}\t1\n")

# Output message and time complete
print(f"{char}There are {len(trioList) / 3} trios.{char}")
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')