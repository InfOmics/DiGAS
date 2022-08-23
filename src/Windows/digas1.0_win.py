#!/usr/bin/env python
# coding: utf-8


import re
import pandas as pd
import numpy as np
import os
import argparse
import csv
import subprocess
import shutil
import math
import random

parser = argparse.ArgumentParser(description='Description')


# BFILE

def bfile(string):
    if not os.path.exists(string + ".bim") or not os.path.exists(string + ".bed") or not os.path.exists(
            string + ".fam"):
        raise argparse.ArgumentTypeError('Not valid input plink binary files.')
    return string


parser.add_argument('--bfile', action="store", dest="bfile", type=bfile, required=True,
                    help='path of the input file in plink format without extension. Bim bed and fam file are required.')


#REGION FILE

# Region file
def regionFile(string):
    if not os.path.exists(string):
        raise argparse.ArgumentTypeError('Region file not exist.')

    filename, file_extension = os.path.splitext(string)

    if file_extension == '.txt':
        p = re.compile(r"\t")

        with open(string) as f:
            reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
            num_cols = len(next(reader))

    elif file_extension == '.csv':
        p = re.compile(r",")

        with open(string) as f:
            reader = csv.reader(f, delimiter=',', skipinitialspace=True)
            num_cols = len(next(reader))

    else:
        raise argparse.ArgumentTypeError('Required txt or csv file.')

    if num_cols == 2 or num_cols == 3:
        res = []
        with open(string) as f:
            for line in f.read().splitlines():
                if not line.startswith('#'):
                    tmp = p.split(line, maxsplit=2)
                    if num_cols == 2:
                        tmp.insert(0, tmp[0] + "-" + tmp[1])
                    res.append(tmp)

        df = pd.DataFrame(columns=['region', 'start', 'end'], data=res)

    else:
        print(num_cols)
        raise argparse.ArgumentTypeError(
            'txt or csv file must have two or three columns. Columns found:' + str(num_cols))

    print('Region file used: ' + string)
    return df.dropna()


parser.add_argument('--region', action="store", type=regionFile, dest="region",required=True, help='Region file.')

# ITERATION
parser.add_argument("--iterations", nargs='?', dest="iteration", type=int, const=1, default=1000,
                    help='number of iterations to calculate the p-values. Default 1000.')


# PHENO
def phenoFile(string):
    if not os.path.exists(string):
        raise argparse.ArgumentTypeError('Pheno file not exist.')

    filename, file_extension = os.path.splitext(string)

    if file_extension == '.txt':
        p = re.compile(r"\t")
    elif file_extension == '.csv':
        p = re.compile(r",")
    else:
        raise argparse.ArgumentTypeError('Required txt or csv file.')

    res = []
    with open(string) as f:
        for line in f.read().splitlines():
            if not line.startswith('#'):
                res.append(p.split(line, maxsplit=1))

    df = pd.DataFrame(columns=['subject', 'status'], data=res)

    print("Pheno file used: " + string)
    return df.dropna()


parser.add_argument('--pheno', action="store", type=phenoFile, dest="pheno",
                    help='optional tab separeted input file. Phenotype of the subjects of the analysis (overwrite data from original fam). It accepts either numeric or categorical labels. Only the subjects in common between fam and pheno are considered in the analysis. WARNING: prints the subjects excluded from the analysis.')

# TEMPORARY DIRECTORY
parser.add_argument("--saveIntermediate", action="store_false",
                    help = 'save all intermediate file in a directory called as --out option: list of subjects and raw data divided by phenotype, complete results (For each significative region generalized allele spectrum with respect to each phenotype category, fold change  with respect to two phenotype categies, p-value), etc  ')



#OUTPUT NAME
parser.add_argument("--out",
                    help="directs the output to a name of your choice (if not deleted, the directory with temporary files will be called in the same way)")


# END PARSER
args = parser.parse_args()

# SELECT ONLY COMMON ELEMENT BETWEEN FAM AND PHENO FILES


# PRINT AND READ


print("Bfile path:" + os.path.join(os.getcwd(), args.bfile))
fam = pd.read_csv(os.path.join(os.getcwd(), (args.bfile + '.fam')), header=None, sep=" ")
fam = fam.dropna()
fam.rename(columns={0: 'FID', 1: 'subject', 5: 'status'}, inplace=True)
print(str(fam.shape[0]) + " subject found in fam file")
# print(fam)

bim = pd.read_csv(os.path.join(os.getcwd(), (args.bfile + '.bim')), header=None, sep="\t")
bim = bim.dropna()
bim.rename(columns={0: 'chr', 1: 'idsnp', 3: 'bp'}, inplace=True)
bim[["chr", "bp"]] = bim[["chr", "bp"]].apply(pd.to_numeric)
print(str(bim.shape[0]) + " SNPs found in bim file")
# print(bim)


if args.pheno is not None:

    print('Select common subject between pheno and fam files.')
    print(str(args.pheno.shape[0]) + " subject found in pheno file")
    famsubj = list(set(fam['subject']))
    indiv = args.pheno[args.pheno['subject'].isin(famsubj)].dropna()
    subj = pd.merge(indiv, fam, on='subject')
    subj = subj[['FID', 'subject', 'status_x']]
    subj.rename(columns={'status_x': 'status'}, inplace=True)

    print(str(indiv.shape[0]) + ' between pheno and fam files.')
    print(str(fam.shape[0] - indiv.shape[0]) + ' subject exluded from FAM file.')
    print(str(args.pheno.shape[0] - indiv.shape[0]) + ' subject exluded from Pheno file.')
else:
    print('All subject in FAM file are included.')
    subj = fam[['FID', 'subject', 'status']].copy()

del fam
# print(subj)

category = subj['status'].value_counts().index.tolist()
category_number = subj['status'].value_counts()

print('Phenotype found: ' + ' '.join([str(c) for c in category]))
print('P-value random iterations set to: ' + str(args.iteration))

if len(category) == 1:
    raise ValueError('Phenotype value shoould be greater then one.')

if args.out is not None:
    folder = str(args.out)
else:
    folder = 'digastmp'

if not os.path.exists(folder):
    os.mkdir(folder)

print('Create plink file.')
for c in category:  ##AA: controllare se è int se funziona il primo controllo
    #subj.loc[subj['status'] == c, ('FID', 'subject')].to_csv(folder+'\'' + str(c) + '_patients.txt', sep='\t',
    #                                                         index=False, header=False)

    subj.loc[subj['status'] == c, ('FID', 'subject')].to_csv(os.path.join(folder, str(c) + '_patients.txt'), sep='\t',
                                                             index=False, header=False)

    # AA rimuovo standard output per warning in adni3 no qc
    subprocess.run(
        ["plink.exe", "--bfile", args.bfile, "--keep", os.path.join(folder, str(c) + "_patients.txt"), "--make-bed", "--out",
         os.path.join(folder,  str(c))], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

print('Selecting region.')

region = args.region
snptogene = []
for i in range(0, region.shape[0]):
    tmp = []
    tmp.append(region.iloc[i, 0])
    chromosome = int(region.iloc[i, 1].split(':')[0])
    start = int(region.iloc[i, 1].split(str(chromosome) + ':')[1])
    end = int(region.iloc[i, 2].split(str(chromosome) + ':')[1])
    chrbim = bim.loc[bim['chr'] == chromosome, ('idsnp', 'bp')]
    chrbim = chrbim.loc[chrbim['bp'].between(start, end, inclusive=False)]

    chrbim['bp'] = str(chromosome) + ':' + chrbim['bp'].astype(str)

    count = chrbim.shape[0]
    pos = chrbim['bp'].str.cat(sep=',')
    name = chrbim['idsnp'].str.cat(sep=',')
    tmp.append(count)
    tmp.append(pos)
    tmp.append(name)
    snpttogene = snptogene.append(tmp)


df = pd.DataFrame(snptogene, columns=['Genes', 'SNPsCount', 'SNPsNames', 'SNPsNames_Converter'])
df = df[df['SNPsCount'] > 0]

#Exit if no SNPs are in region specified
if df.shape[0] == 0:
    raise SystemExit('ERROR: No SNPs in the region specified in regionFile.')

df.to_csv(os.path.join(folder,"SnpToGene_RSID.txt"), sep='\t', index=False)

snpstokeep = df['SNPsNames_Converter']
snpstokeep = snpstokeep.str.split(",").values
snpstokeep = np.concatenate(snpstokeep)
snpstokeep = np.char.lstrip(snpstokeep)

# AA: rifare
path = os.path.join(folder, "SNPS_to_keep.txt")
file = open(path, "w+")
for i in snpstokeep.tolist():
    file.write(i + "\n")
file.close()

del snpstokeep
del df
del bim

# create recoded with selected snps
for c in category:
    subprocess.run(
        ["plink.exe", "--bfile", os.path.join(folder,  str(c)), "--extract", os.path.join(folder, "SNPS_to_keep.txt"), "--recodeA", "--out",
         os.path.join(folder,  str(c) + "_recoded")], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# create merged file AA rifare
with open(os.path.join(folder, "merged.raw"), "w") as merge:
    for i in range(0, len(category)):
        iline = 0
        for line in open(os.path.join(folder,  str(category[i])) + '_recoded.raw'):

            tmp = [v for v in line.strip().split(' ')]
            tmp = tmp[6:]
            if iline == 0:
                if i == 0:  # write header only on first file, snps corrersponds
                    merge.write(' '.join(tmp))
            else:
                for j in range(len(tmp)):
                    if tmp[j] == '2':
                        tmp[j] = '1'
                    elif tmp[j] == 'NA':
                        tmp[j] = '0'
                merge.write('\n' + ' '.join(tmp))
            iline += 1

print("Start analisys.")

def get_hills_by_sum(cat_nof, sample2cat, sample2snp, snp_list, gene2snp):
    snp_count = dict()
    for snp in snp_list:
        snp_count[snp] = {c: 0 for c in cat_nof.keys()}
    for sample, snps in sample2snp.items():
        for snp in snps:
            snp_count[snp][sample2cat[sample]] += 1

    hills = {gene: dict() for gene in gene2snp.keys()}
    for gene in gene2snp.keys():
        for cat in cat_nof.keys():
            hills[gene][cat] = 0
        for snp in gene2snp[gene]:
            for cat in cat_nof.keys():
                hills[gene][cat] += snp_count[snp][cat] / (cat_nof[cat] * len(gene2snp[gene]))
    return hills


snp2gene = dict()
gene2snp = dict()
iline = 0
for line in open(os.path.join(folder, 'SnpToGene_RSID.txt')):
    if iline > 0:
        cc = line.strip().split('\t')
        gene = cc[0]
        if gene not in gene2snp:
            gene2snp[gene] = list()
        snps = cc[3].split(',')
        for snp in snps:
            snp = snp.strip()
            snp2gene[snp] = gene
            gene2snp[gene].append(snp)
    iline += 1

# cn= 226
# mci= 59
# ad= 17


# cat_nof = {'CN':cn,'MCI':mci,'AD':ad}
zip_iterator = zip(category, category_number)
cat_nof = dict(zip_iterator)

# sample2cat = (['CN']*cn) + (['MCI']*mci) + (['AD']*ad)
sample2cat = ([category[0]] * category_number[0])
for i in range(1, len(category)):
    sample2cat = sample2cat + ([category[i]] * category_number[i])

sample2snp = dict()
snp_list = list()

iline = 0
for line in open(os.path.join(folder, 'merged.raw')):
    if iline == 0:
        snp_list = [v[:-2] for v in line.strip().split(' ')]
    else:
        cc = line.strip().split(' ')
        sample2snp[iline - 1] = [snp_list[i] for i in range(len(cc)) if cc[i] == '1']

    iline += 1

nof_rands = int(args.iteration)
random.seed(0)

shuffle_perc = 0.5

file_name = 'DiGASCompleteResults.txt'
f = open(os.path.join(folder,file_name), 'a+')
w = ["Cat1", "Cat2", "Gene", "FC", "Pvalue", "HillCat1", "HillCat2"]
w = ' '.join(w)
f.write(w + '\n')

output_name = folder + '.txt'
out = open(output_name, 'a+')
w = ["Cat1", "Cat2", "Gene", "Pvalue"]
w = ' '.join(w)
out.write(w + '\n')


nof_shuffs = math.ceil(len(sample2cat) * shuffle_perc)
real_hills = get_hills_by_sum(cat_nof, sample2cat, sample2snp, snp_list, gene2snp)

rand_hills = {gene: dict() for gene in gene2snp.keys()}
for gene in gene2snp.keys():
    for cat in cat_nof.keys():
        rand_hills[gene][cat] = [0 for i in range(nof_rands)]


print("Iterations:", end=" ", flush=True)

for r in range(nof_rands):

    if ((r + 1) % 100) == 0:
        print(r + 1, end=" ", flush=True)

    rand_sample2cat = [x for x in sample2cat]
    for s in range(nof_shuffs):
        s1 = random.randint(0, len(rand_sample2cat) - 1)
        s2 = random.randint(0, len(rand_sample2cat) - 1)
        c = rand_sample2cat[s1]
        rand_sample2cat[s1] = rand_sample2cat[s2]
        rand_sample2cat[s2] = c

    hills = get_hills_by_sum(cat_nof, rand_sample2cat, sample2snp, snp_list, gene2snp)
    for gene in hills.keys():
        for cat in cat_nof.keys():
            rand_hills[gene][cat][r] = hills[gene][cat]

counts = dict()
count_nozero = dict()
print(" ")
print('Computing p-value:', end=" ")
for cat1_i in range(len(category)):
    for cat2_i in range(cat1_i + 1, len(category)):
        cat1 = category[cat1_i]
        cat2 = category[cat2_i]
        cca = cat1 + '_' + cat2
        print(cca + ' ', end='')
        counts[cca] = 0
        count_nozero[cca] = 0

        for gene in sorted(real_hills.keys()):
            a = rand_hills[gene][cat1]
            b = rand_hills[gene][cat2]
            diff_rand = [abs(math.log(a[i] + 1) - math.log(b[i] + 1)) for i in range(len(a))]
            diff_real = abs(math.log(real_hills[gene][cat1] + 1) - math.log(real_hills[gene][cat2] + 1))
            pvalue = 0.0
            for i in range(len(diff_rand)):
                if diff_rand[i] >= diff_real:
                    pvalue += 1.0
            pvalue = pvalue / float(len(diff_rand))

            if pvalue <= 0.05:
                counts[cca] += 1
                #Results complete
                w = [cat1, cat2, gene, str(diff_real), str(pvalue), str(real_hills[gene][cat1]),str(real_hills[gene][cat2])]
                w = ' '.join(w)
                f.write(w + '\n')
                #Final results
                w = [cat1, cat2, gene, str(pvalue)]
                w = ' '.join(w)
                out.write(w + '\n')
            if diff_real > 0.0:
                count_nozero[cca] += 1
f.close()
out.close()
if args.saveIntermediate:
    shutil.rmtree(folder)

print(" ")
print("Done!")
