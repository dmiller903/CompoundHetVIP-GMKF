import gzip
import re
import os
import time
from sys import argv
import concurrent.futures

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# argv information
inputFile = argv[1]
pathToFiles = argv[2]
numCores = int(argv[3])
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileSet = set()
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        trioFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio_liftover_parsed.vcf.gz"
        fileSet.add(trioFileName)

plinkFileSet = set()
# Separate combined trio files and individual participant files by chromosome
def separateByChr(file):
    with gzip.open(file, "rt") as vcf:
        fileName = re.findall(r"([\w\-/_]+)_liftover_parsed\.vcf\.gz", file)[0]
        outputName = f"{fileName}_"
        header = ""
        tmpFileSet = set()
        chromosomeSet = set()
        chromosomeNumber = ""
        
        for line in vcf:
            if line.startswith("#"):
                header = header + line
            elif not line.startswith("#") and line.split("\t")[0] not in chromosomeSet:
                chromosomeNumber = line.split("\t")[0]
                os.system(f"rm {outputName}{chromosomeNumber}.*")
                with open(f"{outputName}{chromosomeNumber}.vcf", "w") as chromosome:
                    chromosome.write(header)
                    chromosome.write(line)
                    chromosomeSet.add(chromosomeNumber)
                    tmpFileSet.add(f"{outputName}{chromosomeNumber}.vcf")
            else:
                with open(f"{outputName}{chromosomeNumber}.vcf", "a") as chromosome:
                    chromosome.write(line)
        
        return(tmpFileSet)

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    tmpFileList = executor.map(separateByChr, fileSet)
    for tmpList in tmpFileList:
        for tmpSet in tmpList:
            plinkFileSet.add(tmpSet)

plinkFileList = list(plinkFileSet)
plinkFileList.sort()

# Create bed, bim files for each chromosome of each trio
def createPlink(file):
    filePath, famFolder, sampleFolder, sampleFile, chrNumber = re.findall(r"([\w\-/_]+)\/([\w\-]+)\/([\w\-]+)\/([\w\-]+)_(chr[\w]+)\.?.*\.?.*\.vcf", file)[0]
    outputName = f"{filePath}/{famFolder}/{sampleFolder}/{sampleFile}_{chrNumber}"

    os.system(f"/plink2 --vcf {file} --fam {filePath}/{famFolder}/{famFolder}_trio.fam --make-bed --out {outputName}")

with concurrent.futures.ProcessPoolExecutor(max_workers=numCores) as executor:
    executor.map(createPlink, plinkFileList)

# Output Time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')