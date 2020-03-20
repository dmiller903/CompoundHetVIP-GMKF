import re
import os
from sys import argv
import time

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Download annotation files and CADD files if they haven't been downloaded already
if not os.path.exists("/usr/local/share/gemini/gemini_data/hg19.vista.enhancers.20131108.bed.gz.tbi"):
    os.system("gemini update --dataonly --extra cadd_score")

# argv information
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]
numCores = int(argv[2])

# Get disease name based on path
diseaseName = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Load annotated file into a GEMINI database
os.system(f"gemini load -v {pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated.vcf \
-p {pathToFiles}/{diseaseName}.fam -t snpEff --cores {numCores} {pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db")

# Print output information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')