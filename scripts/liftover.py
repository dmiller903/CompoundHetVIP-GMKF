import re
import os
import time
from sys import argv
import concurrent.futures

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Input file or list of files
inputFile = argv[1]
pathToFiles = argv[2]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Download necessary reference files if needed
if not os.path.exists("/references/hg38ToHg19.over.chain"):
    os.system("wget --no-check-certificate \
        https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9ddd0ba7bc73000ce87e38/?zip= -O /tmp/references.zip \
        && unzip /tmp/references.zip -d /tmp/references \
        && rm /tmp/references.zip \
        && /root/miniconda2/bin/bgzip -d /tmp/references/human_g1k_v37_modified.fasta.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/hg38ToHg19.over.chain.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.fasta.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.dict.gz \
        && /root/miniconda2/bin/bgzip -d /tmp/references/Homo_sapiens_assembly38.fasta.fai.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}"):
            os.system(f"mv {file} /references/")
        elif fileName == "readme":
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# Create a list of file(s) to be liftedover
fileList = []
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
        trioFileName = f"{pathToFiles}/{sampleFamilyId}/{sampleFamilyId}_trio/{sampleFamilyId}_trio.vcf.gz"
        if os.path.exists(trioFileName):
            fileList.append(trioFileName)

# Liftover file(s)
def liftoverFiles(file):
    fileFolder, fileName = re.findall("([\w\-\/_]+)\/([\w\-_]+)\.?.*\.?.*\.gz", file)[0]
    os.system(f"gatk IndexFeatureFile -F {file}")
    os.system(f'gatk --java-options "-Xmx4g" GenotypeGVCFs -R /references/Homo_sapiens_assembly38.fasta -V {file} -O {fileFolder}/{fileName}_genotyped.vcf.gz')
    os.system(f"java -jar /root/miniconda2/share/picard-2.21.1-0/picard.jar LiftoverVcf I={fileFolder}/{fileName}_genotyped.vcf.gz O={fileFolder}/{fileName}_liftover.vcf.gz \
    CHAIN=/references/hg38ToHg19.over.chain R=/references/human_g1k_v37_modified.fasta REJECT={fileFolder}/{fileName}_rejected_variants.vcf")

for i in range(0, len(fileList), 24):
    fileListSlice = fileList[i:(i+24)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
        executor.map(liftoverFiles, fileListSlice)

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')