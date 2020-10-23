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
outputFile = argv[2]

disease = outputFile.split("/")[0]

# Dictionary to store needed information in
outputDict = {}

# Obtain information from tsv file and add to initial dictionary
with open(tsvTable) as table:
    tableColumnNames = table.readline()
    tableColumnNames = tableColumnNames.strip().split("\t")
    # Information needed from this file are the "File Name", "Family Id", "Participants ID", and "Proband" (Yes, or No)
    # and "Biospecimen ID".
    familyIdIndex = tableColumnNames.index("Family Id")
    tableSampleIdIndex = tableColumnNames.index("Participants ID")
    probandIndex = tableColumnNames.index("Proband")
    diagnosisIndex = tableColumnNames.index("Diagnosis (Source Text)")
    # Add to the initial dictionary where the key is the "Participants ID" as this is common across all input files.
    # The value is a list where value[0] the "File Name", value[1] is "Participants ID, value[2] is "Family Id", 
    # and value[3] is "Proband" as "Yes" or "No"
    for sample in table:
        sample = sample.rstrip().split("\t")
        proband = sample[probandIndex]
        diagnosis = sample[diagnosisIndex]
        if proband == "--" and diagnosis != "--":
            proband = "Yes"
        elif proband == "--" and diagnosis == "--":
            proband = "No"
        if proband == "Yes":
            outputDict[sample[tableSampleIdIndex]] = [sample[familyIdIndex], proband, sample[diagnosisIndex]]

# Output information to outputFile
with open(outputFile, 'w') as output:
    numberProbands = 0
    output.write("file_name\tdiagnosis\n")
    total_excluded = 0
    for key, value in outputDict.items():
        familyName = value[0]
        diagnosis = value[-1]
        proband = value[1]
        if os.path.exists(f"{disease}/gVCF/{familyName}") and proband == "Yes":
            output.write(f"{key}\t{diagnosis}\n")
            numberProbands += 1

print(str(numberProbands))
# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')