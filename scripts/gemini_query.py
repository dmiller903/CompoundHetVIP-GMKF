import re
import os
from sys import argv
import time

#Keep track of when the script began
startTime = time.time()
char = '\n' + ('*' * 70) + '\n'

# argv information
pathToFiles = argv[1]
if pathToFiles.endswith("/"):
    pathToFiles = pathToFiles[0:-1]

# Get disease name based on path and set MED HIGH variables
diseaseName = re.findall(r"([\w\-_]+)\/", pathToFiles)[0]
med = "'MED'"
high = "'HIGH'"

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf005_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.005 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf005_cadd15.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf01_cadd20.tsv')

# Filter CH Variants with minor allele frequency <= 0.01 and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_aaf01_cadd15.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_cadd20.tsv')

# Filter CH Variants with no regard for minor allele frequency and cadd >= 15
os.system(f'gemini comp_hets --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_ch_impactHM_cadd15.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.005 anc cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf005_cadd20.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.005 anc cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.005 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf005_cadd15.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.01 anc cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf01_cadd20.tsv')

# Filter for de novo variants with minor allele frequncey <= 0.01 anc cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and aaf_1kg_all <= 0.01 and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_aaf01_cadd15.tsv')

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 20
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=20)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_cadd20.tsv')

# Filter for de novo variants with no regard for minor allele frequency and cadd >= 15
os.system(f'gemini de_novo --columns "chrom, start, end, vcf_id, ref, alt, gene, impact_severity, aaf_1kg_all, cadd_scaled, impact, biotype" \
--filter "(impact_severity = {high} or is_lof = 1) or (impact_severity = {med} and cadd_scaled >=15)" \
{pathToFiles}/{diseaseName}_phased_mcmc_samples_annotated_cadd.db \
> {pathToFiles}/{diseaseName}_de_novo_impactHM_cadd15.tsv')

#Output time information
timeElapsedMinutes = round((time.time()-startTime) / 60, 2)
timeElapsedHours = round(timeElapsedMinutes / 60, 2)
print(f'{char}Done. Time elapsed: {timeElapsedMinutes} minutes ({timeElapsedHours} hours){char}')