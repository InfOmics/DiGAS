# DiGAS

A novel computational model for identifying genomic elements (from single exon to entire genomic regions)  that are linked to a given phenotypic condition (i.e. a disease) investigating for single nucleotide polymorphisms (SNP).

Genomic elements are investigated by means of a newly introduced descriptor of genomic information, the generalized allele spectrum. In contrast to allele spectrum, the novel descriptor takes into account the complete set of SNPs of a region at once. Thus, frequency is computed at region level rather than at SNP level. Then, we introduce the Differential Generalized Allele Spectrum (DiGAS), which captures the differences in frequency allele spectra between healthy and ill sets (control and case respectively). Statistical significance of a gene is then evaluated by means of a differential analysis between the healthy and ill portions of the population.

DiGAS is distributed as a command-line tool developed in Python and it is available for both Windows and Unix systems. It uses PLINK to retrieve data about genotyping and takes as input the coordinates of genomic regions to be analyzed.


## Installation and usage

DiGAS can be easly downloaded from GitHub `src` directory (Windows or Unix version). A tested PLINK tool (version 1.9) is also contained in both directory in zip format, you have to extract it before running DiGAS.

### System requirements

DiGAS can run on Windows or Unix system where Python 3 (or higher) was previously installed.

DiGAS depends on PLINK. PLINK executable must be contained in the same folder as the DiGAS executable.
A compressed PLINK version is contained in the `src` directory for both Windows and Unix versions. 
Alternatively you can download it from the official website: https://www.cog-genomics.org/plink/.

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
The script `digas.py` will run the DiGAS workflow.
```
python digas.py [-h] --bfile BFILE --region REGION [--iterations [ITERATION]] [--pheno PHENO] [--saveIntermediate] [--out OUT]
```

BFILE and REGION files are mandatory. 

#### -h
or `--help` show the help message. 

#### --bfile 
Path of the PLINK binary files ***witouth extensions*** contained data to be analyzed. `bim`, `bed` and `fam` files are required.

#### --region
TO DO

#### --iterations
Number of iterations to calculate p-value. Default 1000.

#### --pheno
Optional tab-separeted (txt) or comma-separated (csv) input file used to overwrite subjects phenotype from original fam. It accepts either numeric or categorical labels. 
Only subjects in common between fam and pheno files are included in the analysis. DiGAS command-line tool prints the number of subjects excluded from the analysis.

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
Save all intermediate file in a directory called as `--out` argument (or `digastmp` as default): 
- list of subjects and raw data divided by phenotype;
- `SnpToGene_RSID.txt`: a tab-separated file where for each genomic region (specified in `--region` argument) is reported the number of SNPs located in that region (i.e. *2*), SNPs coordinates in *chr:position* format (i.e. *10:43156390,10:43206700* ) and SNPs names in rsID format (i.e. *rs7100293,rs10128322*);
- `SNPS_to_keep.txt`: a list of SNPs located in defined genomic regions;
- 'DiGASCompleteResults.txt': For each DiGAS significant genomic region reports the generalized allele spectrum with respect to each phenotype category, fold change with respect to two phenotype  categies and the p-value.

#### --out
Name of the output file. Default `digastmp.txt`.



## Data

To test DiGAS we used just genotyping and demographic data avaiable at ADNI (http://adni.loni.usc.edu) data portal.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Support

You can contact directly the authors by their e-mail addresses.

## Authors and acknowledgment

Antonino Aparo (antonino.aparo@univr.it) and Vincenzo Bonnici (vincenzo.bonnici@univr.it) created the workflow. 

## License
[MIT](https://choosealicense.com/licenses/mit/)

