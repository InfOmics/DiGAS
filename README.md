# DiGAS

A novel computational model for identifying genomic elements (from single exon to entire genomic regions)  that are linked to a given phenotypic condition (i.e. a disease) investigating for single nucleotide polymorphisms (SNPs).

Genomic elements are investigated by means of a newly introduced descriptor of genomic information, the generalized allele spectrum. In contrast to allele spectrum, the novel descriptor takes into account the complete set of SNPs of a region at once. Thus, frequency is computed at region level rather than at SNP level. Then, we introduce the Differential Generalized Allele Spectrum (DiGAS), which captures the differences in frequency allele spectra between healthy and ill sets (control and case respectively). Statistical significance of a gene is then evaluated by means of a differential analysis between the healthy and ill portions of the population.

DiGAS is distributed as a command-line tool developed in Python and it is available for both Windows and Unix systems. It uses PLINK to retrieve data about genotyping and takes as input the coordinates of genomic regions to be analyzed.


## Installation and usage

DiGAS can be easliy downloaded from GitHub `src` directory (Windows or Unix version). A tested PLINK tool (version 1.9) is also contained in both directories in zip format, you have to extract it before running DiGAS.

### System requirements

DiGAS can run on Windows or Unix systems where Python 3 (or higher) was previously installed.

DiGAS depends on PLINK. PLINK executable must be contained in the same folder as the DiGAS executable.
A compressed PLINK version is contained in the `src` directory for both Windows and Unix versions. 
Alternatively, you can download it from the official website: https://www.cog-genomics.org/plink/.

DiGAS depends on a number of Python packages. To install all the required packages:

```
pip3 install pandas
pip3 install numpy
pip3 install argparse
pip3 install csv
pip3 install subprocess
pip3 install shutil
pip3 install re
```

### Usage
The script `digas.py` runs the DiGAS workflow.
```
python digas.py [-h] --bfile BFILE --region REGION [--iterations [ITERATION]] [--pheno PHENO] [--saveIntermediate] [--out OUT]
```

BFILE and REGION files are mandatory. 

#### -h
or `--help` shows the help message. 

#### --bfile 
Path of the PLINK binary files ***without extensions*** contained data to be analyzed. `bim`, `bed` and `fam` files are required.

#### --region
Tab-separated (txt) or comma-separated (csv) input files containing genomic position of genomic regions to be analyzed in *chr:position* format.
The name of the genomic region is not mandatory, while the *start* and *end* genomic positions are mandatory.


```
region1	1:154432419	1:154457855
region2	2:160220164	2:160220337
region3	3:14379478	3:14379480
region4	4:3813929	4:3825268
region5	5:1231231	5:1345235
...
```

```
1:154432419,1:154457855
2:160220164,2:160220337
3:14379478,3:14379480
4:3813929,4:3825268
5:1231231,5:1345235
...
```

#### --iterations
Number of iterations to calculate the p-value. Default 1000.

#### --pheno
Optional tab-separated (txt) or comma-separated (csv) input file used to overwrite subjects phenotype from original fam. It accepts either numeric or categorical labels. 
Only subjects in common between fam and pheno files are included in the analysis. DiGAS command-line tool prints the number of subjects excluded from the analysis.

Examples of valid pheno files are as follow:

```
subject1  1  
subject2  3
subject3  1
subject4  3
subject5  2
...
```

```
subject1,healty  
subject2,ill
subject3,healty
subject4,ill
subject5,unknown
...
```

#### --saveIntermediate
Save all intermediate files in a directory called as `--out` argument (or `digastmp` as default): 
- list of subjects and raw data divided by phenotype;
- `SnpToGene_RSID.txt`: a tab-separated file where for each genomic region (specified in `--region` argument) is reported the number of SNPs located in that region (i.e. *2*), SNPs coordinates in *chr:position* format (i.e. *10:43156390,10:43206700* ) and SNPs names in rsID format (i.e. *rs7100293,rs10128322*);
- `SNPS_to_keep.txt`: a list of SNPs located in defined genomic regions;
- 'DiGASCompleteResults.txt': For each DiGAS significant genomic region reports the generalized allele spectrum with respect to each phenotype category, fold change with respect to two phenotype categories and the p-value.


#### --out
Name of the output file. Default `digastmp.txt`.
For each significant genomic region, it shows the name of the region, the phenotype categories and the respective p-values.

```
Cat1  Cat2  Gene  Pvalue
ill healthy region1 0.038
ill unkown region1 0.047
ill healthy region3 0.038
ill unkown region4 0.038
...
```

## Data
To test DiGAS we used genotyping and demographic data available at ADNI (http://adni.loni.usc.edu) data portal.
Any kind of PLINK binary format file can be used.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests as appropriate.

## Support
You can contact directly the authors by their e-mail addresses.

## Authors and acknowledgement
Antonino Aparo (antonino.aparo@univr.it) and Vincenzo Bonnici (vincenzo.bonnici@unipr.it) created the workflow. 

## Citation
Submitted.

## License
[MIT](https://choosealicense.com/licenses/mit/)

