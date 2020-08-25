import glob
import re
import os
import concurrent.futures
from sys import argv
import time

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# argv information
inputFile = argv[1]
pathToFiles = argv[2]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Download reference files if necessary
if not os.path.exists("/references/1000GP_Phase3/1000GP_Phase3.sample"):
    os.system("wget --no-check-certificate https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz -P /tmp/references \
    && tar zxf /tmp/references/1000GP_Phase3.tgz -C /tmp/references/ \
    && rm /tmp/references/1000GP_Phase3.tgz \
    && wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5db8af02f3bb87000b85b76e/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references/1000GP_Phase3 \
    && rm /tmp/references.zip")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
    os.system("chmod -R 777 /references/1000GP_Phase3")

# Modify genetic map files to format required for beagle if not already done
if not os.path.exists("/references/1000GP_Phase3/genetic_map_chr1_combined_b37_beagle.txt"):
    def updateFiles(file):
        fileName = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][0]
        chrom = re.findall(r"([\w/_]+genetic_map_chr(\w+)_combined_b37)\.txt", file)[0][1]
        beagleOutput = f"{fileName}_beagle.txt"
        with open(file) as inputFile, open(beagleOutput, 'w') as beagleOut:
            header = inputFile.readline()
            header = "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
            for line in inputFile:
                beagleLineList = line.rstrip().split(" ")
                beagleLine = f"{chrom} {beagleLineList[1]} {beagleLineList[2]} {beagleLineList[0]}\n"
                beagleOut.write(beagleLine)
        os.system(f"chmod 777 {beagleOutput}")
    for file in glob.glob("/references/1000GP_Phase3/genetic_map_chr*_combined_b37.txt"):
        updateFiles(file)

#Create a list of file(s) that need to have unplaced and multiallelic sites removed
fileDict = dict()
with open(inputFile) as sampleFile:
    header = sampleFile.readline()
    headerList = header.rstrip().split("\t")
    fileNameIndex = headerList.index("file_name")
    familyIdIndex = headerList.index("family_id")
    sampleIdIndex = headerList.index("sample_id")
    chromosomes = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",\
"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"}
    for sample in sampleFile:
        sampleData = sample.rstrip("\n").split("\t")
        fileName = sampleData[fileNameIndex]
        sampleFamilyId = sampleData[familyIdIndex]
        sampleId = sampleData[sampleIdIndex]
        if sampleFamilyId not in fileDict:
            fileDict[sampleFamilyId] = set()
            for chromosome in chromosomes:
                trioFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio_{chromosome}.vcf"
                fileDict[sampleFamilyId].add(trioFileName)

for trio in fileDict:
    trioChr = fileDict[trio]
    def runEagle(file):
        chromosome = re.findall(r"[\w\-\/_]+chr([0-9][0-9]?)", file)[0]
        outputName = "/tmp/" + re.findall(r"\/([\w\-_]+chr[0-9][0-9]?)", file)[0] + "_update.vcf"
        filePrefix = re.findall(r"([\w\-\/_]+chr[0-9][0-9]?)", file)[0] + "_beagle_phased"
        
        # VCF files must first have chr# changed to # only
        with open(file) as inputFile, open(outputName, 'w') as outputFile:
            for line in inputFile:
                line = line.replace(f"chr{chromosome}", f"{chromosome}")
                outputFile.write(line)

        # Updated VCF needs to be bgzipped and tabixed
        os.system(f"/root/miniconda2/bin/bgzip {outputName}")
        os.system(f"/root/miniconda2/bin/tabix {outputName}.gz")

        # Phase with Eagle
        os.system(f"java -Xmx40g -jar /beagle.25Nov19.28d.jar gt={outputName}.gz \
        out={filePrefix} \
        chrom={chromosome} \
        map=/references/1000GP_Phase3/genetic_map_chr{chromosome}_combined_b37_beagle.txt \
        ref=/references/1000GP_Phase3/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        impute=false")
    with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
        executor.map(runEagle, trioChr)

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')