import gzip
import re
import statistics
import argparse
import os
import time

# Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# Argparse Information
parser = argparse.ArgumentParser(description='Adds GDI scores and gene lengths to a GEMINI query file. If working with \
controlled data, this script can also anonymize family ID and patient ID information.')

parser.add_argument('gemini_query_file', help='File generate by GEMINI query')
parser.add_argument('output_file', help='Name of output file')
parser.add_argument('--anonymize', help="If controlled data is used, and you need to transfer the output file off of a \
secure server, this script can anonymize identifying information", default='y')

args = parser.parse_args()

#Create variables of each argument from argparse
geminiQuery = args.gemini_query_file
outputFile = args.output_file
anonymizePatients = args.anonymize

# Function to get total gene lengths and gene exon lengths
def getLengths(tsvFile):
    """
    Create a gene dictionary where the key is the gene name, and the value is the length of the gene. Also, create  a cds 
    dictionary where the key is the gene name, and the value is a dictionary where the key is the exon number and the value
    is a list of exon lengths
    """
    cdsDict = {}
    geneDict = {}
    with open(tsvFile) as gencodeParsed:
        for line in gencodeParsed:
            if "#" not in line and (("\tCDS\t" in line and 'transcript_type "protein_coding"' in line) or "\tgene\t" in line) and 'gene_type "protein_coding"' in line:
                lineList = line.rstrip("\n").split("\t")
                end = lineList[4]
                start = lineList[3]
                geneOrCDS = lineList[2]
                geneName = re.findall(r' gene_name "([\w\-\.]+)";', line)[0]
                if re.search(r'\.', geneName):
                    continue
                if re.search(r' exon_number (\d+)', line):
                    exonNumber = re.findall(r' exon_number (\d+)', line)[0]
                if geneOrCDS == "gene" and geneName not in geneDict:
                    geneDict[geneName] = int(end) - int(start)
                elif geneOrCDS == "CDS" and geneName not in cdsDict:
                    cdsDict[geneName] = {exonNumber: [int(end) - int(start)]}
                elif geneOrCDS == "CDS" and geneName in cdsDict and exonNumber in cdsDict[geneName]:
                    cdsDict[geneName][exonNumber].append(int(end) - int(start))
                elif geneOrCDS == "CDS" and geneName in cdsDict and exonNumber not in cdsDict[geneName]:
                    cdsDict[geneName][exonNumber] = [int(end) - int(start)]

    """
    Iterate through cdsDict and average each exon length. When finished, the cdsDict will have the keys as gene names and
    the values as dictionaries where the key is the exon number and the value is the average length for that exon
    """
    for key, value in cdsDict.items():
        for key2, value2 in value.items():
            cdsDict[key][key2] = sum(value2) / len(value2)

    """
    Iterate through the cdsDict and create a new dict 'cdsAvgDict' where the key is the gene name and the value is a list of
    all the average exon lengths
    """
    cdsAvgDict = {}
    for key, value in cdsDict.items():
        for key2, value2 in value.items():
            if key not in cdsAvgDict:
                cdsAvgDict[key] = [value2]
            else:
                cdsAvgDict[key].append(value2)

    # Create a new dictionary "cdsSumDict" where the key is the gene name and the value is the sum of all the exon lengths
    cdsSumDict = {}
    for key, value in cdsAvgDict.items():
        cdsSumDict[key] = sum(value)
    # Return dictionaries
    return(geneDict, cdsSumDict)

# Update add_GDI_raw.py to python3 syntax
with open("/add_GDI_raw.py") as raw, open("/add_GDI.py", "w") as newFile:
    for line in raw:
        if "print" not in line:
            newFile.write(line)

# Create a list of genes that are in the GEMINI query and a file of genes
geneList = []
with open(geminiQuery) as queryFile, open("/gene_list.txt", 'w') as geneFile:
    header = queryFile.readline()
    headerList = header.rstrip("\n").split("\t")
    geneIndex = headerList.index("gene")
    for line in queryFile:
        lineList = line.rstrip("\n").split("\t")
        gene = lineList[geneIndex]
        geneList.append(gene)
        geneFile.write(gene + "\n")

# Create a new gencode file that only includes the genes that are in the CH review
with gzip.open("/gencode.v33.annotation.gtf.gz", 'rt') as gencode, open("/tmp/gencode_parsed.tsv", 'w') as output:
    for line in gencode:
        for gene in geneList:
            if "#" not in line and (("\tCDS\t" in line and 'transcript_type "protein_coding"' in line) or "\tgene\t" in line) and gene in line:
                output.write(line)

# Create a gene length dictionary and a cds gene length dictionary using the parsed GENCODE file
geneDict, cdsSumDict = getLengths("/tmp/gencode_parsed.tsv")

"""
Create a dictionary where the key is the gene name, the value is a list where the 0th item is the cds length and the 
1st item is the total gene length
"""
geneLengths = {}
for key, value in sorted(geneDict.items()):
    if key in cdsSumDict:
        geneLengths[key] = [cdsSumDict[key], value]
    else:
        geneLengths[key] = ["-", value]

gdiDict = {}
# Generate GDI scores for each gene
os.system("cd / && python3 add_GDI.py && cd /proj")
with open("/GDI_output.txt") as gdiScores:
    for line in gdiScores:
        lineList = line.rstrip("\n").split("\t")
        gdiDict[lineList[0]] = [lineList[1], lineList[2]]

# Create new sample names based on family id
familyDict = {}
with open(geminiQuery) as queryFile:
    header = queryFile.readline()
    headerList = header.rstrip("\n").split("\t")
    familyIndex = headerList.index("family_id")
    familyCount = 1
    for line in queryFile:
        lineList = line.rstrip("\n").split("\t")
        familyId = lineList[familyIndex]
        if familyId not in familyDict:
            familyDict[familyId] = "patient_{}".format(familyCount)
            familyCount += 1

if anonymizePatients == "y":
    # Anonymize patients, Add gene lengths, and Add GDI values to output file
    with open(geminiQuery) as queryFile, open(outputFile, 'w') as outputFile:
        header = queryFile.readline()
        headerList = header.rstrip("\n").split("\t")
        geneIndex = headerList.index("gene")
        familyIndex = headerList.index("family_id")
        familyMembersIndex = headerList.index("family_members")
        samplesIndex = headerList.index("samples")
        newHeaderList = []
        for i, item in enumerate(headerList):
            if i not in [familyIndex, familyMembersIndex, samplesIndex]:
                newHeaderList.append(item)
        header = "\t".join(newHeaderList)
        header = f"patient\tgender\t{header}\tcds_gene_length\ttotal_gene_length\tGDI-raw\tGDI-Phred\n"
        outputFile.write(header)
        for line in queryFile:
            lineList = line.rstrip("\n").split("\t")
            gene = lineList[geneIndex]
            patient = lineList[familyIndex]
            patient = familyDict[patient]
            gender = re.findall(r"affected;(\w+)\)", line)[0]
            newLineList = []
            for i, item in enumerate(lineList):
                if i not in [familyIndex, familyMembersIndex, samplesIndex]:
                    newLineList.append(item)
            line = "\t".join(newLineList)
            if gene in geneLengths and gene in gdiDict:
                outputFile.write(f"{patient}\t{gender}\t{line}\t{geneLengths[gene][0]}\t{geneLengths[gene][1]}\t{gdiDict[gene][0]}\t{gdiDict[gene][1]}\n")
            elif gene in geneLengths and gene not in gdiDict:
                outputFile.write(f"{patient}\t{gender}\t{line}\t{geneLengths[gene][0]}\t{geneLengths[gene][1]}\tNA\tNA\n")
            elif gene in gdiDict and gene not in geneLengths:
                outputFile.write(f"{patient}\t{gender}\t{line}\tNA\tNA\t{gdiDict[gene][0]}\t{gdiDict[gene][1]}\n")
            else:
                outputFile.write(f"{patient}\t{gender}\t{line}\tNA\tNA\tNA\tNA\n")
elif anonymizePatients == "n":
    #Add gene lengths, and Add GDI values to output file
    with open(geminiQuery) as queryFile, open(outputFile, 'w') as outputFile:
        header = queryFile.readline().rstrip("\n")
        headerList = header.rstrip("\n").split("\t")
        geneIndex = headerList.index("gene")
        header = f"{header}\tcds_gene_length\ttotal_gene_length\tGDI-raw\tGDI-Phred\n"
        outputFile.write(header)
        for line in queryFile:
            lineList = line.rstrip("\n").split("\t")
            gene = lineList[geneIndex]
            if gene in geneLengths and gene in gdiDict:
                outputFile.write(f"{line}\t{geneLengths[gene][0]}\t{geneLengths[gene][1]}\t{gdiDict[gene][0]}\t{gdiDict[gene][1]}\n")
            elif gene in geneLengths and gene not in gdiDict:
                outputFile.write(f"{line}\t{geneLengths[gene][0]}\t{geneLengths[gene][1]}\tNA\tNA\n")
            elif gene in gdiDict and gene not in geneLengths:
                outputFile.write(f"{line}\tNA\tNA\t{gdiDict[gene][0]}\t{gdiDict[gene][1]}\n")
            else:
                outputFile.write(f"{line}\tNA\tNA\tNA\tNA\n")

# Print message and how long the previous steps took
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')