import re
import os
from sys import argv
import time

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Download annotation files if needed
if not os.path.exists("/snpEff/./data/GRCh37.75/sequence.HSCHR6_MHC_SSTO.bin"):
    os.system("java -jar /snpEff/snpEff.jar download -v GRCh37.75")

# argv information
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Get disease name based on path
diseaseName = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]

# Annotate the vt trimmed file
os.system(f"java -Xmx20g -jar /snpEff/snpEff.jar GRCh37.75 -v \
{pathToFiles}/{diseaseName}_phased_samples_vt.vcf > \
{pathToFiles}/{diseaseName}_phased_samples_annotated.vcf")

# Print output information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')