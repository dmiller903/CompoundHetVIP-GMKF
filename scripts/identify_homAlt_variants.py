import argparse
import re
import os
from sys import argv
import time
import concurrent.futures

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# argv information
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]
# Get disease name based on path and set MED HIGH variables
diseaseName = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]
familyFile = f"{pathToFiles}/{diseaseName}.fam"
inputCadd = float(argv[2])
inputCaddStr = argv[2]
inputFile = f"{pathToFiles}/{diseaseName}_gemini.tsv"
inputMaf = argv[3]
if inputMaf != "None":
    inputMaf = float(argv[3])
    inputMafStr = argv[3].replace("0.", "")
    homAltFile = f"{pathToFiles}/{diseaseName}_homAlt_cadd{inputCaddStr}_maf{inputMafStr}.tsv"
else:
    homAltFile = f"{pathToFiles}/{diseaseName}_homAlt_cadd{inputCaddStr}.tsv"

#Function to get convert sample genotype from alpha to numeric
def getNumericGenotype(genotype, ref, alt):
    if "|" in genotype and "." not in genotype:
        genotypeList = genotype.split("|")
        firstAllele = ""
        secondAllele = ""
        if genotypeList[0] == ref:
            firstAllele = "0"
        elif genotypeList[0] == alt:
            firstAllele = "1"
        else:
            firstAllele = "."
        if genotypeList[1] == ref:
            secondAllele = "0"
        elif genotypeList[1] == alt:
            secondAllele = "1"
        else:
            secondAllele = "."
        newGenotype = f"{firstAllele}|{secondAllele}"
        return(newGenotype)
    else:
        return(".|.")

#Function to grab information from header of input file
def getHeaderInfo(headerList):
    startIndex = headerList.index("start")
    geneIndex = headerList.index("gene")
    refIndex = headerList.index("ref")
    altIndex = headerList.index("alt")
    impactIndex = headerList.index("impact_severity")
    caddIndex = headerList.index("cadd_scaled")
    mafIndex = headerList.index("aaf_1kg_all")
    lofIndex = headerList.index("is_lof")
    exonicIndex = headerList.index("is_exonic")
    rsIndex = headerList.index("rs_ids")
    clinVarIndex = headerList.index("clinvar_sig")
    samples = headerList[15:]
    return(startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, rsIndex, clinVarIndex, samples)

#Function to grab information from line of input file
def getLineInfo(lineList):
    start = lineList[startIndex]
    gene = lineList[geneIndex]
    ref = lineList[refIndex]
    alt = lineList[altIndex]
    impact = lineList[impactIndex]
    cadd = lineList[caddIndex]
    maf = lineList[mafIndex]
    lof = lineList[lofIndex]
    exonic = lineList[exonicIndex]
    rs = lineList[rsIndex]
    clinVar = lineList[clinVarIndex]
    return(start, gene, ref, alt, impact, cadd, maf, lof, exonic, rs, clinVar)

def iterateThroughSamples():
    for sampleIndex in sampleIndexes:
        sample = headerList[sampleIndex]
        genotype = lineList[sampleIndex]
        newGenotype = getNumericGenotype(genotype, ref, alt)
        if gene not in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene] = [newGenotype]
            samplePositions[sample][gene] = [start]
        elif gene in sampleGenotype[sample] and "." not in newGenotype:
            sampleGenotype[sample][gene].append(newGenotype)
            samplePositions[sample][gene].append(start)

# Create a .tsv that has all pertinent information for compound heterozygous identification
impactSeverity = "'LOW'"
if not os.path.exists(f"{pathToFiles}/{diseaseName}_gemini.tsv"):
    os.system(f'gemini query --header -q "select chrom, start, vcf_id, ref, alt, gene, is_exonic, impact_severity, \
        is_lof, aaf_1kg_all, cadd_scaled, impact, biotype, rs_ids, clinvar_sig, (gts).(*) from variants where impact_severity != {impactSeverity}" \
        {pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
        > {pathToFiles}/{diseaseName}_gemini.tsv')

# Use fam file to create a list of samples, list of parents, and a parent dictionary where each key is a parent ID and value is sample ID
parentDict = {}
parentList = []
patientList = []
familyDict = {}
with open(familyFile) as familyFile:
    for line in familyFile:
        lineList = line.rstrip("\n").split("\t")
        if lineList[-1] is "2":
            parentDict[f"gts.{lineList[2]}"] = f"gts.{lineList[1]}"
            parentDict[f"gts.{lineList[3]}"] = f"gts.{lineList[1]}"
            patientList.append(f"gts.{lineList[1]}")
            familyDict[f"gts.{lineList[1]}"] = [f"gts.{lineList[2]}", f"gts.{lineList[3]}"]
        else:
            parentList.append(f"gts.{lineList[1]}")
"""
Iterate through inputFile in order to create two dictionaries: sampleGenotype and samplePositions. The key of 
sampleGenotype is the sample ID and value is a dictionary where the key is a gene and the value is a list of all 
genotypes ("0|1", "1|0", or "0|0") for that gene that meet specific CADD score, minor allele frequency, and impact 
severity criteria. The samplePositions has the same information, except the list for each gene is genotype positions.
"""
sampleGenotype = {}
samplePositions = {}
sampleIndexes = []
with open(inputFile) as geminiFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, rsIndex, clinVarIndex, samples = getHeaderInfo(headerList)
    for sample in samples:
        sampleIndexes.append(headerList.index(sample))
        sampleGenotype[sample] = {}
        samplePositions[sample] = {}    
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, maf, lof, exonic, rs, clinVar = getLineInfo(lineList)
        if cadd != "None" and maf != "None":
            if ((impact == "HIGH" or lof == "1") or (impact == "MED" and float(cadd) >= inputCadd)) and float(maf) <= inputMaf:
                iterateThroughSamples()
        elif cadd == "None" and maf == "None":
            if impact == "HIGH" or lof == "1":
                iterateThroughSamples()
        elif cadd != "None" and maf == "None":
            if (impact == "HIGH" or lof == "1") or (impact == "MED" and float(cadd) >= inputCadd):
                iterateThroughSamples()
        elif cadd == "None" and maf != "None":
            if (impact == "HIGH" or lof == "1") and float(maf) <= inputMaf:
                iterateThroughSamples()
print("Sample Dictionaries Created.")
print("Sample Dictionaries Created.")

"""
Use sampleGenotype and samplePositions to generate a new dictionaries where the key is the sample ID and the value is 
a dictionary where the key is a gene and the value is a list of genotypes (or positions) where homozygous Alt variant(s) are found.
"""
homAltPositionDict = {}
homAltGenotypeDict = {}
for patient in patientList:
    homAltPositionDict[patient] = {}
    homAltGenotypeDict[patient] = {}
    parent1 = familyDict[patient][0]
    parent2 = familyDict[patient][1]
    for gene, genotypes in sampleGenotype[patient].items():
        for i, genotype in enumerate(genotypes):
            positionList = samplePositions[patient][gene]
            position = positionList[i]
            #This part helps eliminate genotypes being added to the homAlt list where either parent is homozygous recessive
            parentGenotype1 = ""
            parentGenotype2 = ""
            if gene in samplePositions[parent1] and position in samplePositions[parent1][gene]:
                parentPosIndex = samplePositions[parent1][gene].index(position)
                parentGenotype1 = sampleGenotype[parent1][gene][parentPosIndex]
            if gene in samplePositions[parent2] and position in samplePositions[parent2][gene]:
                parentPosIndex = samplePositions[parent2][gene].index(position)
                parentGenotype2 = sampleGenotype[parent2][gene][parentPosIndex]
            if genotype == "1|1" and parentGenotype1 != "1|1" and parentGenotype2 != "1|1":
                if gene not in homAltPositionDict[patient]:
                    homAltPositionDict[patient][gene] = [position]
                    homAltGenotypeDict[patient][gene] = [genotype]
                elif gene in homAltPositionDict[patient]:
                    homAltPositionDict[patient][gene].append(position)
                    homAltGenotypeDict[patient][gene].append(genotype)

print("homozygous Alt variant dictionaries created.")

#Iterate through the input file and use the homAltPositionDict in order to output homozygous Alt variant data for each sample
with open(inputFile) as geminiFile, open(homAltFile, "w") as outputFile:
    header = geminiFile.readline()
    headerList = header.rstrip("\n").split("\t")
    startIndex, geneIndex, refIndex, altIndex, impactIndex, caddIndex, mafIndex, lofIndex, exonicIndex, rsIndex, clinVarIndex, samples = getHeaderInfo(headerList)
    columnInfo = headerList[0:15]
    newHeader = "\t".join(columnInfo) + "\tgenotype\tsample\n"
    outputFile.write(newHeader)
    
    sampleIndexes = []
    for patient in patientList:
        patientIndex = headerList.index(patient)
        sampleIndexes.append(patientIndex)
    for line in geminiFile:
        lineList = line.rstrip("\n").split("\t")
        start, gene, ref, alt, impact, cadd, maf, lof, exonic, rs, clinVar = getLineInfo(lineList)
        for sampleIndex in sampleIndexes:
            sample = headerList[sampleIndex]
            parent1 = familyDict[sample][0]
            parent2 = familyDict[sample][1]
            if gene in homAltPositionDict[sample] and start in homAltPositionDict[sample][gene]:
                genotype = lineList[sampleIndex]
                numericGenotype = getNumericGenotype(genotype, ref, alt)
                if numericGenotype == "1|1":
                    columnInfo = lineList[0:15]
                    columnStr = "\t".join(columnInfo)
                    newLine = f"{columnStr}\t{numericGenotype}\t{sample.lstrip('gts.')}\n"
                    outputFile.write(newLine)

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')