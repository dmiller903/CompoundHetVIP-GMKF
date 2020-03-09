import re
import os
import concurrent.futures
from sys import argv
import gzip
import time

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

#Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]
diseaseName = re.findall(r"[\w_\-]+", pathToFiles)[0]

#Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileDict = dict()
concatFiles = list()
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        if sampleFamilyId not in fileDict:
            fileDict[sampleFamilyId] = list()
            # concatFileName is what each trio will be named after each trio is combined into one trio
            concatFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio_phased_combined.vcf.gz"
            concatFiles.append(concatFileName)
            for chromosome in chromosomes:
                trioFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio_{chromosome}_phased_reverted.vcf"
                fileDict[sampleFamilyId].append(trioFileName)

#Concatenate individual chromosomes into one file for each trio
def concatMerge(trio):
    files = fileDict[trio]

    for index, file in enumerate(files):
        os.system(f"bgzip -f {file} && tabix -fp vcf {file}.gz")
        files[index] = f"{file}.gz"


    fileName = re.findall(r"([\w\-\/_]+\/[\w\-_]+)_chr[A-Z0-9][A-Z0-9]?_phased_reverted\.vcf", files[0])[0]
    outputName = f"{fileName}_phased_combined.vcf"
    files = " ".join(files)
    os.system(f"bcftools concat {files} -o {outputName}")
    os.system(f"bgzip -f {outputName} && tabix -fp vcf {outputName}.gz")

with concurrent.futures.ProcessPoolExecutor(max_workers=35) as executor:
    executor.map(concatMerge, fileDict)

# Merge all phased, concatenated, trio files into one    
concatFilesString = " ".join(concatFiles)
outputName = f"{pathToFiles}/{diseaseName}_phased_samples.vcf"
os.system(f"bcftools merge -m both {concatFilesString} -o {outputName}")
os.system(f"bgzip -f {outputName} && tabix -fp vcf {outputName}.gz")

# Create a merged family file
# create a proband dictionary where the key is the sampleId and the value is the familyId
# also create a parent dictionary where the key is familyId and the value is a dictionary that has a key of the sampleId and value of gender
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

# Create a dictionary where each sample has the rest of the family information needed for the family file
sampleDict = dict()
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
        paternal = ""
        maternal = ""
        if probandStatus == "Yes":
            familyDict = parentDict[sampleFamilyId]
            for key, value in familyDict.items():
                if value == "1":
                    paternal = key
                else:
                    maternal = key
            sampleDict[sampleId] = f"{sampleFamilyId}\t{sampleId}\t{paternal}\t{paternal}\t{gender}\t2\n"
        else:
            sampleDict[sampleId] = f"{sampleFamilyId}\t{sampleId}\t0\t0\t{gender}\t1\n"
            
# create a sample list in the order of the vcf file
with gzip.open(f"{pathToFiles}/{diseaseName}_phased_samples.vcf.gz", "rt") as vcfFile:
    for line in vcfFile:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            sampleList = line.rstrip().split("\t")[9:]
        else:
            break

# use the sample order in the list to output each sample in order as found in the vcf file
with open(f"{pathToFiles}/{diseaseName}.fam", "w") as outputFile:
    for sample in sampleList:
        outputFile.write(sampleDict[sample])

#Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')