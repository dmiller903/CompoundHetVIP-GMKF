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

# Create a dictionary of trio files that need to be phased
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
                trioFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio_{chromosome}"
                fileDict[sampleFamilyId].add(trioFileName)

for trio in fileDict:
    trioChr = fileDict[trio]
    def runShapeit(file):
        chromosome = re.findall(r"[\w\-\/_]+(chr[0-9][0-9]?)", file)[0]
        # Check for alignment issues between sample and reference panel
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check -B {file} --output-log {file}_check -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
        --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
        /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3")
        # If alignment issues are found, remove problematic positions while phasing using the exclude file from check step
        if os.path.exists(f"{file}_check.snp.strand.exclude"):
            os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {file} --output-log {file}_phased.log -O {file}_phased -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
            --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
            /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc --exclude-snp {file}_check.snp.strand.exclude --force --seed 123456789")
        # If no alignment issues are found, do not remove positions while phasing
        else:
            os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B {file} --output-log {file}_phased.log -O {file}_phased -M /references/1000GP_Phase3/genetic_map_{chromosome}_combined_b37.txt \
            --input-ref /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.hap.gz /references/1000GP_Phase3/1000GP_Phase3_{chromosome}.legend.gz \
            /references/1000GP_Phase3/1000GP_Phase3.sample --thread 3 --no-mcmc --force --seed 123456789")
        # Convert phased files to vcf files and bgzip output vcf
        os.system(f"/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps {file}_phased --output-log {file}_phased_vcf.log --output-vcf {file}_phased.vcf")

    with concurrent.futures.ProcessPoolExecutor(max_workers=23) as executor:
        executor.map(runShapeit, trioChr)

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')