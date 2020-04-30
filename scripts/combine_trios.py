import gzip
import re
import os
import time
from sys import argv
import concurrent.futures
import glob

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
numCores = int(argv[3])
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Download reference files if needed
if not os.path.exists("/references/Homo_sapiens_assembly38.fasta"):
    os.system("wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9f54d2a7bc73000ee99fd6/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references \
    && rm /tmp/references.zip \
    && gzip -d /tmp/references/*.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# Create a dictionary of files that need to be combined into one vcf file
fileDict = {}
familyList = []
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        sampleId = sampleData[sampleIdIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        actualFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleId}/{sampleId}_parsed.vcf.gz"
        outputName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio.vcf.gz"
        if sampleFamilyId not in fileDict and os.path.exists(f"{pathToFiles}/{sampleFamilyId}") and not os.path.exists(f"{outputName}"):
            fileDict[sampleFamilyId] = [actualFileName]
            familyList.append(sampleFamilyId)
        elif sampleFamilyId in fileDict and os.path.exists(f"{pathToFiles}/{sampleFamilyId}") and not os.path.exists(f"{outputName}"):
            fileDict[sampleFamilyId].append(actualFileName)

probandDict = {}
parentDict = {}
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    probandIndex = headerList.index("proband")
    genderIndex = headerList.index("sex")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        probandStatus = sampleData[probandIndex]
        gender = sampleData[genderIndex]
        if probandStatus == "Yes":
            probandDict[sampleId] = sampleFamilyId
        else:
            if sampleFamilyId not in parentDict:
                parentDict[sampleFamilyId] = {sampleId: gender}
            else:
                parentDict[sampleFamilyId][sampleId] = gender

# Create fam files
def createFamFiles(proband):
    familyId = probandDict[proband]
    familyDict = parentDict[familyId]
    paternal = ""
    maternal = ""
    outputString = ""
    sampleDict = {}
    for key, value in familyDict.items():
        if value == "1":
            paternal = key
        else:
            maternal = key
    with open(inputFile) as sampleFile:
        header = sampleFile.readline()
        headerList = header.rstrip().split("\t")
        fileNameIndex = headerList.index("file_name")
        familyIdIndex = headerList.index("family_id")
        sampleIdIndex = headerList.index("sample_id")
        probandIndex = headerList.index("proband")
        genderIndex = headerList.index("sex")
        for sample in sampleFile:
            sampleData = sample.rstrip("\n").split("\t")
            fileName = sampleData[fileNameIndex]
            sampleFamilyId = sampleData[familyIdIndex]
            sampleId = sampleData[sampleIdIndex]
            probandStatus = sampleData[probandIndex]
            gender = sampleData[genderIndex]
            if probandStatus == "Yes" and familyId == sampleFamilyId:
                sampleDict[sampleId] = f"{sampleFamilyId}\t{sampleId}\t{paternal}\t{maternal}\t{gender}\t2\n"
            elif probandStatus == "No" and familyId == sampleFamilyId:
                sampleDict[sampleId] = f"{sampleFamilyId}\t{sampleId}\t0\t0\t{gender}\t1\n"
    with open(f"{pathToFiles}/{familyId}/{familyId}_trio.fam", "w") as outputFile:
        for key, value in sorted(sampleDict.items()):
            outputFile.write(value)

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    executor.map(createFamFiles, probandDict)

# Use GATK to combine all trios into one vcf
def combineTrios(trio):
    files = fileDict[trio]
    fileString = ""
    os.system(f"mkdir {pathToFiles}/{trio}/{trio}_trio")
    outputName = f"{pathToFiles}/{trio}/{trio}_trio/{trio}_trio.vcf.gz"
    for file in files:
        fileString += f"-V {file} "
        os.system(f"/root/miniconda2/bin/gatk IndexFeatureFile -F {file}")
    os.system(f"/root/miniconda2/bin/gatk CombineGVCFs -R /references/Homo_sapiens_assembly38.fasta {fileString} -O {outputName}")

for i in range(0, len(familyList), numCores):
    familyListSlice = familyList[i:(i+numCores)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
        executor.map(combineTrios, familyListSlice)

# Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')
