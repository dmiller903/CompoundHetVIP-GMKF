import re
import os
import time
from sys import argv
import concurrent.futures

startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

#Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

#Create a list of file(s) to be liftedover
fileSet = set()
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        trioFileName = "{}/{}/{}_trio/{}_trio.vcf.gz".format(pathToFiles, sampleFamilyId, sampleFamilyId, sampleFamilyId)
        fileSet.add(trioFileName)

#Liftover file(s)
def liftoverFiles(file):
    fileFolder, fileName = re.findall("([\w\-\/_]+)\/([\w\-_]+)\.?.*\.?.*\.gz", file)[0]
    os.system("gatk IndexFeatureFile -F {}".format(file))
    os.system('gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {} -O {}/{}_genotyped.vcf.gz'.format(file, fileFolder, fileName))
    os.system("java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I={}/{}_genotyped.vcf.gz O={}/{}_liftover.vcf.gz \
    CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={}/{}_rejected_variants.vcf".format(fileFolder, fileName, fileFolder, fileName, fileFolder, fileName))

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    executor.map(liftoverFiles, fileSet)

#Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print('{}Liftover Complete. Time elapsed: {} minutes ({} hours){}'.format(char, timeElapsedMinutes, timeElapsedHours, char))