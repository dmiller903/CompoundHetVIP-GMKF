import re
import os
from sys import argv
import time
import glob

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Download reference files if needed
if not os.path.exists("/references/human_g1k_v37.fasta"):
    os.system("wget --no-check-certificate \
    https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5dc57b1f7f37e3000ecaed96/?zip= -O /tmp/references.zip \
    && unzip /tmp/references.zip -d /tmp/references \
    && rm /tmp/references.zip \
    && gzip -d /tmp/references/human_g1k_v37.fasta.gz")
    for file in glob.glob("/tmp/references/*"):
        fileName = file.split("/")[-1]
        if not os.path.exists(f"/references/{fileName}") and fileName != "readme":
            os.system(f"mv {file} /references/")
    os.system("chmod 777 /references/*")

# argv information
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

#Get disease name based on path
diseaseName = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Use VT to split, trim and left align the phased samples.
os.system(f"/root/miniconda2/bin/vt decompose -s {pathToFiles}/{diseaseName}_phased_mcmc_samples.vcf.gz \
| /root/miniconda2/bin/vt normalize -n -r /references/human_g1k_v37.fasta - > \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_vt.vcf")

# Output message and time complete
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')