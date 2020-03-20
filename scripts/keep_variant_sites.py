import gzip
import re
import os
import time
from sys import argv
import concurrent.futures
import math

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argv information 
inputFile = argv[1]
pathToFiles = argv[2]
numCores = int(argv[3])
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Create a list of proband files that need to have non-variant sites removed. Create a list of parent files that need sites removed
probandDict = {}
familyDict = {}
with open(inputFile) as tsvFile:
    header = tsvFile.readline()
    header = header.rstrip().split("\t")
    fileNameIndex = header.index("file_name")
    familyIdIndex = header.index("family_id")
    probandIndex = header.index("proband")
    sampleIdIndex = header.index("sample_id")
    for sample in tsvFile:
        sample = sample.rstrip().split("\t")
        if os.path.exists(f"{pathToFiles}/{sample[fileNameIndex]}"):
            if sample[probandIndex] == "Yes":
                probandDict[sample[familyIdIndex]] = [sample[sampleIdIndex], sample[fileNameIndex]]
                if sample[familyIdIndex] not in familyDict:
                    familyDict[sample[familyIdIndex]] = [[sample[sampleIdIndex], sample[fileNameIndex]]]
                else:
                    familyDict[sample[familyIdIndex]].append([sample[sampleIdIndex], sample[fileNameIndex]])
            else:
                if sample[familyIdIndex] not in familyDict:
                    familyDict[sample[familyIdIndex]] = [[sample[sampleIdIndex], sample[fileNameIndex]]]
                else:
                    familyDict[sample[familyIdIndex]].append([sample[sampleIdIndex], sample[fileNameIndex]])

# Create new familyDict, called trioDict, to include only the files with trios
trioDict = {}
familyList = []
for key, value in familyDict.items():
    if len(value) == 3:
        trioDict[key] = value
        familyList.append(key)

# Filter each proband file, remove  variants-only sites, create a dictionary of variant-only sites
def filterVariantOnly(familyID, probandID, fileName):
    os.system(f"mkdir {pathToFiles}/{familyID}")
    os.system(f"mkdir {pathToFiles}/{familyID}/{probandID}")
    outputName = f"{pathToFiles}/{familyID}/{probandID}/{probandID}_parsed.vcf"
    positionDict = {}
    with gzip.open(f"{pathToFiles}/{fileName}", 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith('#'):
                parsed.write(line.encode())
            elif "END" not in line:
                parsed.write(line.encode())
                line = line.split("\t")
                chrom = line[0]
                pos = line[1]
                if chrom not in positionDict:
                    positionDict[chrom] = {pos}
                else:
                    positionDict[chrom].add(pos)
    finalName = f"{outputName}.gz"
    os.system(f"zcat {outputName} | /root/miniconda2/bin/bgzip > {finalName}")
    os.system(f"rm {outputName}")
    return(positionDict)

#Filter each parent file for sites that occur in proband of that family
def filterParents(familyID, parentID, fileName, positionDict):
    os.system(f"mkdir {pathToFiles}/{familyID}/{parentID}")
    outputName = f"{pathToFiles}/{familyID}/{parentID}/{parentID}_parsed.vcf"
    with gzip.open(f"{pathToFiles}/{fileName}", 'rt') as gVCF, gzip.open(outputName, 'wb') as parsed:
        for line in gVCF:
            if line.startswith("#"):
                parsed.write(line.encode())
            else:
                lineList = line.split("\t")
                chrom = lineList[0]
                pos = lineList[1]
                if pos in positionDict[chrom]:
                    parsed.write(line.encode())
                else:
                    if "END" in line:
                        for i in range(int(pos), int(lineList[7].lstrip("END=")) + 1):
                            if str(i) in positionDict[chrom]:
                                parsed.write(line.encode())
    finalName = f"{outputName}.gz"
    os.system(f"zcat {outputName} | /root/miniconda2/bin/bgzip > {finalName}")
    os.system(f"rm {outputName}")

# Iterate through familyList and remove variant sites from proband first, while creating a position dictionary. 
# Then use the position dictionary to iterate through each parent file and keep positions that are in common with proband.
def filterFiles(familyID):
    for sample in trioDict[familyID]:
        if sample == probandDict[familyID]:
            probandID = sample[0]
            fileName = sample[1]
            positionDict = filterVariantOnly(familyID, probandID, fileName)
        
    for sample in trioDict[familyID]:
        if sample != probandDict[familyID]:
            parentID = sample[0]
            fileName = sample[1]
            filterParents(familyID, parentID, fileName, positionDict)
    
# Use concurrent.futures to filter through multiple trios at a time using the filterFiles function
for i in range(0, len(familyList), numCores):
    familyListSlice = familyList[i:(i+numCores)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
        executor.map(filterFiles, familyListSlice)

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours) {char}')